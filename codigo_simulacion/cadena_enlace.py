"""
cadena_enlace.py — Cadena de enlace reutilizable (PRBS->PAM-4->FFE->canal->CTLE->LPF->DFE).

Módulo compartido por el barrido de tasa (y el re-DoE). Porta los bloques del
notebook `simulacion_canal_rcEjecutadoActual.ipynb`, pero está optimizado para
las MILES de evaluaciones que exigen 30 realizaciones × barrido de tasa:

  - El canal y el CTLE se aplican en el DOMINIO DE LA FRECUENCIA (multiplicación
    por H(f)), equivalente a `signal.lsim` para un sistema LTI en estado nulo
    pero mucho más rápido.
  - La alineación de muestreo se hace con UNA correlación por FFT (no 32 barridos
    de fase con np.correlate).

Modelo de ruido fiel al artículo: AWGN sigma_in referido a la ENTRADA del RX
(antes del CTLE) + LPF de BW finito. BER semi-analítico (Q-función) con ajuste
Gaussiano del ojo.
"""

import numpy as np
import scipy.signal as signal
from scipy.special import erfc

# --- Constantes de la cadena (co-diseño y receptor del artículo) ---------- #
SPS = 32
SIGMA_IN = 0.10
CTLE_P1, CTLE_P2, CTLE_ADC = 2 * np.pi * 15e9, 2 * np.pi * 30e9, 2.5
NIVELES_PAM4 = np.array([-3.0, -1.0, 1.0, 3.0])

# Co-diseño "heredado del RC" (DoE fijado en checkpoint.md)
W_Z_RC = 2 * np.pi * 2.35e9
C_POST_RC = -0.301

# Rango de búsqueda del DoE (idéntico al notebook)
WZ_MIN, WZ_MAX = 1e9, 10e9
CPOST_MIN, CPOST_MAX = -0.40, -0.05


# =========================================================================== #
# Generación / modulación / FFE  (idénticos al notebook)
# =========================================================================== #
def generar_prbs(n_bits, orden=23, semilla=None):
    TAPS = {7: [7, 6], 9: [9, 5], 11: [11, 9], 15: [15, 14], 23: [23, 18], 31: [31, 28]}
    taps = TAPS[orden]
    reg = [1] * orden if semilla is None else list(semilla)
    bits = np.empty(n_bits, dtype=np.uint8)
    for i in range(n_bits):
        bits[i] = reg[-1]
        realim = 0
        for t in taps:
            realim ^= reg[t - 1]
        reg = [realim] + reg[:-1]
    return bits


def modular_pam4(bits):
    if len(bits) % 2 != 0:
        bits = np.append(bits, 0)
    dibits = bits.reshape(-1, 2)
    mapa = {0: -3.0, 1: -1.0, 3: 1.0, 2: 3.0}
    return np.array([mapa[d[0] * 2 + d[1]] for d in dibits])


def aplicar_ffe(simbolos, c_pre, c_main, c_post):
    taps = np.array([c_pre, c_main, c_post])
    taps = taps / np.sum(np.abs(taps))
    return signal.lfilter(taps, 1.0, simbolos)


def aplicar_rx_lpf(y, ts_val, fc_hz, orden=3):
    fs_val = 1.0 / ts_val
    wn = min(fc_hz / (fs_val / 2.0), 0.99)
    b, a = signal.butter(orden, wn, btype="low")
    return signal.filtfilt(b, a, y)


def Q(x):
    return 0.5 * erfc(np.asarray(x, dtype=float) / np.sqrt(2.0))


# =========================================================================== #
# Canal y CTLE en frecuencia
# =========================================================================== #
def aplicar_canal_freq(x, ts_val, H_func):
    """y = IFFT( H_canal(f)·FFT(x) )."""
    n = len(x)
    f = np.fft.rfftfreq(n, ts_val)
    return np.fft.irfft(np.fft.rfft(x) * H_func(f), n=n)


def ctle_response(f, w_z, w_p1=CTLE_P1, w_p2=CTLE_P2, A_dc=CTLE_ADC):
    """H_CTLE(jw) = A_dc (1 + jw/wz) / [(1+jw/wp1)(1+jw/wp2)]."""
    jw = 1j * 2.0 * np.pi * np.asarray(f, dtype=float)
    return A_dc * (1.0 + jw / w_z) / ((1.0 + jw / w_p1) * (1.0 + jw / w_p2))


def aplicar_ctle_freq(y, ts_val, w_z, a_dc=CTLE_ADC):
    """Aplica el CTLE como multiplicación en frecuencia (equivalente a lsim)."""
    n = len(y)
    f = np.fft.rfftfreq(n, ts_val)
    return np.fft.irfft(np.fft.rfft(y) * ctle_response(f, w_z, A_dc=a_dc), n=n)


# =========================================================================== #
# Alineación (1 correlación por FFT) y métricas de ojo / BER
# =========================================================================== #
def alinear_full(y, sym, sps_val, n_score=400):
    """Muestrea 1 valor/símbolo en la FASE y el retardo óptimos.

    1) Retardo entero (en símbolos) grueso por correlación FFT a resolución
       completa.  2) Fase de muestreo óptima (0..sps-1) por máxima correlación
       del flujo submuestreado con los símbolos reales — evita muestrear en el
       borde del símbolo (cruce del ojo)."""
    n = len(y)
    sym = np.asarray(sym, dtype=float)
    x_up = np.repeat(sym, sps_val)
    x_up = x_up[:n] if len(x_up) >= n else np.pad(x_up, (0, n - len(x_up)))
    cc = np.fft.irfft(np.fft.rfft(y - y.mean()) * np.conj(np.fft.rfft(x_up - x_up.mean())), n=n)
    off0 = int(np.argmax(cc[: n // 2]))     # retardo positivo (circular)
    lag0 = off0 // sps_val                   # retardo aproximado en símbolos

    ns = min(n_score, len(sym))
    b = sym[:ns] - sym[:ns].mean()
    nb = np.linalg.norm(b)
    mejor = None
    for ph in range(sps_val):
        ds = y[ph::sps_val]
        for lag in (lag0 - 1, lag0, lag0 + 1):
            if lag < 0 or lag + ns > len(ds):
                continue
            a = ds[lag:lag + ns]
            a = a - a.mean()
            na = np.linalg.norm(a)
            if na <= 0 or nb <= 0:
                continue
            score = abs(float(np.dot(a, b)) / (na * nb))
            if mejor is None or score > mejor[0]:
                mejor = (score, ph, lag)
    if mejor is None:
        return y[::sps_val][:len(sym)], sym[:len(y[::sps_val])]
    _, ph, lag = mejor
    ds = y[ph::sps_val]
    k = min(len(ds) - lag, len(sym))
    return ds[lag:lag + k], sym[:k]


def estadisticas_niveles(muestras, tx, niveles=NIVELES_PAM4):
    mu = np.zeros(4); s = np.zeros(4); cnt = np.zeros(4, dtype=int)
    for k, lv in enumerate(niveles):
        sel = muestras[np.isclose(tx, lv)]
        cnt[k] = len(sel)
        if len(sel) >= 2:
            mu[k] = sel.mean(); s[k] = sel.std()
        elif len(sel) == 1:
            mu[k] = sel[0]
    return mu, s, cnt


def _decidir_pam4(yk, escala):
    niveles = NIVELES_PAM4 * escala
    return NIVELES_PAM4[int(np.argmin(np.abs(yk - niveles)))]


def evaluar_enlace(y, tx_syms, sps_val, dfe_ntaps=0, n_descartar=100):
    """Alinea, aplica DFE opcional y estima BER semi-analítico + apertura de ojo.
    El ruido ya está en la forma de onda `y` (sigma_ruido implícito = 0).

    DFE DATA-AIDED: los taps se estiman por mínimos cuadrados directamente sobre
    las muestras alineadas ( m_k ≈ h0·a_k + Σ h_i·a_{k-i} ), de modo que coinciden
    con la ISI real en la FASE de muestreo de los datos. Evita el desajuste de
    fase diseño↔datos que hacía al DFE inyectar ISI fantasma y diverger."""
    m, tx = alinear_full(y, tx_syms, sps_val)
    m, tx = m[n_descartar:], tx[n_descartar:]
    if len(m) < 50:
        return {"ber": 1.0, "eye_min": -99.0}

    if dfe_ntaps and dfe_ntaps > 0 and len(m) > dfe_ntaps + 20:
        nt = dfe_ntaps
        N = len(m)
        # Regresor: columnas tx[k], tx[k-1], ..., tx[k-nt]
        Xreg = np.column_stack([tx[nt - i: N - i] for i in range(nt + 1)])
        h, *_ = np.linalg.lstsq(Xreg, m[nt:N], rcond=None)
        h0 = h[0] if abs(h[0]) > 1e-12 else 1.0     # ganancia del cursor (datos)
        post = h[1:]                                 # post-cursores en voltios
        hist = [0.0] * nt
        m_eq = np.empty_like(m)
        for kk in range(len(m)):
            fb = sum(post[i] * hist[i] for i in range(nt))
            yk = m[kk] - fb
            m_eq[kk] = yk
            hist = [_decidir_pam4(yk, h0)] + hist[:-1]
        m = m_eq

    mu, s_isi, cnt = estadisticas_niveles(m, tx)
    if np.any(cnt == 0):
        return {"ber": 1.0, "eye_min": -99.0}
    o = np.argsort(mu); mu, s_isi, cnt = mu[o], s_isi[o], cnt[o]
    thr = (mu[:-1] + mu[1:]) / 2.0
    eps = 1e-300
    Pe = np.zeros(4)
    Pe[0] = Q((thr[0] - mu[0]) / (s_isi[0] + eps))
    Pe[1] = Q((mu[1] - thr[0]) / (s_isi[1] + eps)) + Q((thr[1] - mu[1]) / (s_isi[1] + eps))
    Pe[2] = Q((mu[2] - thr[1]) / (s_isi[2] + eps)) + Q((thr[2] - mu[2]) / (s_isi[2] + eps))
    Pe[3] = Q((mu[3] - thr[2]) / (s_isi[3] + eps))
    pesos = cnt / cnt.sum()
    ser = float(np.sum(pesos * Pe)); ber = ser / 2.0
    ap = (mu[1:] - 3 * s_isi[1:]) - (mu[:-1] + 3 * s_isi[:-1])
    return {"ber": ber, "eye_min": float(np.min(ap))}


def disenar_dfe_zf(H_func, w_z, c_post, sps_val, ts_val, fc_rx, n_taps=2, n_sym=256, pre=16):
    """DFE-ZF desde la respuesta al pulso del enlace (FFE->canal->CTLE->LPF_RX).

    CLAVE: se aplica el MISMO LPF de RX que ve la señal de datos. Sin él, el
    sobreimpulso de alta frecuencia que el CTLE (mal adaptado al canal dispersivo)
    introduce queda en el pulso de diseño; `argmax` toma ese pico espurio como
    cursor y los taps salen disparatados -> el DFE diverge. Con el LPF, el pulso
    de diseño coincide con el de decisión y el cursor es el real.

    Los post-cursores se toman estrictamente tras el pico global, en resolución
    completa ( y[pico + i·sps] )."""
    x = np.zeros(n_sym); x[pre] = 1.0
    ffe = aplicar_ffe(x, 0.0, 1.0, c_post)
    xu = np.repeat(ffe, sps_val)
    yc = aplicar_canal_freq(xu, ts_val, H_func)
    y = aplicar_ctle_freq(yc, ts_val, w_z)
    y = aplicar_rx_lpf(y, ts_val, fc_rx)      # mismo front-end que la señal de datos

    # Cursor en la REJILLA de símbolo: fase que maximiza el pico (no el argmax
    # de resolución completa, que puede caer fuera de la rejilla de muestreo).
    mejor = None
    for ph in range(sps_val):
        p = y[ph::sps_val]
        ci = int(np.argmax(np.abs(p)))
        if mejor is None or abs(p[ci]) > abs(mejor[2]):
            mejor = (ph, ci, p[ci])
    ph, ci, cursor = mejor
    if abs(cursor) < 1e-12:
        return np.zeros(n_taps)
    p = y[ph::sps_val]
    # RATIOS de ISI normalizados (post-cursor / cursor), adimensionales.
    return np.array([(p[ci + i] / cursor) if ci + i < len(p) else 0.0
                     for i in range(1, n_taps + 1)])


# =========================================================================== #
# Evaluación de una configuración completa (A/B/C/D) a una tasa dada
# =========================================================================== #
def evaluar_config(H_func, sym, Rs, w_z, c_post, n_dfe, rng, sps=SPS, sigma_in=SIGMA_IN, a_dc=CTLE_ADC):
    ts = 1.0 / (Rs * sps)
    fc = Rs

    # A: canal sin EQ
    y_a = aplicar_canal_freq(np.repeat(sym, sps), ts, H_func)
    yA = aplicar_rx_lpf(y_a + rng.normal(0, sigma_in, len(y_a)), ts, fc)
    rA = evaluar_enlace(yA, sym, sps)

    # B: FFE -> canal
    ffe = aplicar_ffe(sym, 0.0, 1.0, c_post)
    y_chan = aplicar_canal_freq(np.repeat(ffe, sps), ts, H_func)
    yB = aplicar_rx_lpf(y_chan + rng.normal(0, sigma_in, len(y_chan)), ts, fc)
    rB = evaluar_enlace(yB, sym, sps)

    # C: FFE + CTLE (ruido ANTES del CTLE)
    nin = rng.normal(0, sigma_in, len(y_chan))
    yC = aplicar_rx_lpf(aplicar_ctle_freq(y_chan + nin, ts, w_z, a_dc), ts, fc)
    rC = evaluar_enlace(yC, sym, sps)

    # D: + DFE data-aided (taps estimados de los datos en su fase de muestreo)
    rD = evaluar_enlace(yC, sym, sps, dfe_ntaps=n_dfe)

    return {"A": rA, "B": rB, "C": rC, "D": rD}


# =========================================================================== #
# DoE-LHS propio (sin pyDOE) y utilidad de tasa viable interpolada
# =========================================================================== #
def lhs_center(d, n, rng):
    H = np.zeros((n, d))
    for j in range(d):
        H[:, j] = (rng.permutation(n) + 0.5) / n
    return H


def doe_optimizar(H_func, sym, Rs_op, n_dfe, n_muestras=30, seed=42):
    """Re-optimiza (w_z, c_post) minimizando el BER del escenario C (FFE+CTLE)."""
    rng = np.random.default_rng(seed)
    dis = lhs_center(2, n_muestras, rng)
    wz = 2 * np.pi * (WZ_MIN + dis[:, 0] * (WZ_MAX - WZ_MIN))
    cp = CPOST_MIN + dis[:, 1] * (CPOST_MAX - CPOST_MIN)
    mejor = None
    for i in range(n_muestras):
        r = evaluar_config(H_func, sym, Rs_op, wz[i], cp[i], n_dfe,
                           np.random.default_rng(2026))
        if mejor is None or r["C"]["ber"] < mejor["ber_c"]:
            mejor = {"w_z": float(wz[i]), "c_post": float(cp[i]),
                     "ber_c": r["C"]["ber"]}
    return mejor


def viable_interp(bers, gbx, umbral=1e-2):
    """Tasa (GBaud) donde el BER cruza el umbral, interpolando en log10(BER)."""
    b = np.clip(np.asarray(bers, float), 1e-300, 1.0)
    lb = np.log10(b); lt = np.log10(umbral); below = lb <= lt
    if not below.any():
        return 0.0
    if below.all():
        return float(gbx[-1])
    cruce = None
    for i in range(len(gbx) - 1):
        if below[i] and not below[i + 1]:
            frac = (lt - lb[i]) / (lb[i + 1] - lb[i])
            cruce = gbx[i] + frac * (gbx[i + 1] - gbx[i])
    return float(cruce) if cruce is not None else float(gbx[below][-1])
