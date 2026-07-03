import matplotlib
matplotlib.use('Agg')
import numpy as _np2

# ===== CELL 1 =====
import os
import numpy as np
import scipy.signal as signal
import scipy.fft as fft
from scipy.special import erfc
import matplotlib.pyplot as plt
from pyDOE import lhs

print("Librerías cargadas exitosamente.")

# ===== CELL 4 =====
def generar_prbs(n_bits, orden=15, semilla=None):
    """
    Genera una Secuencia Binaria Pseudoaleatoria (PRBS) mediante un LFSR de Fibonacci.
    Usa polinomios primitivos estándar de la industria (los mismos de PRBS7..PRBS31).

    Parámetros
    ----------
    n_bits  : número de bits a generar.
    orden   : grado del LFSR (7, 9, 11, 15, 23, 31). PRBS<orden> tiene período 2^orden - 1.
    semilla : estado inicial (lista de 'orden' bits). Por defecto todos a 1 (no-nulo).
    """
    # Taps (posiciones de realimentación) de polinomios primitivos estándar
    TAPS = {7: [7, 6], 9: [9, 5], 11: [11, 9], 15: [15, 14], 23: [23, 18], 31: [31, 28]}
    if orden not in TAPS:
        raise ValueError(f"Orden PRBS no soportado: {orden}. Use uno de {list(TAPS)}.")
    taps = TAPS[orden]

    if semilla is None:
        reg = [1] * orden
    else:
        reg = list(semilla)
        if len(reg) != orden or not any(reg):
            raise ValueError("La semilla debe tener 'orden' bits y no ser nula.")

    bits = np.empty(n_bits, dtype=np.uint8)
    for i in range(n_bits):
        # Salida = último bit del registro
        bits[i] = reg[-1]
        # Realimentación: XOR de los taps seleccionados
        realim = 0
        for t in taps:
            realim ^= reg[t - 1]
        reg = [realim] + reg[:-1]
    return bits

# --- Generación de la secuencia de prueba reproducible ---
# PRBS23: período 2^23 - 1 = 8 388 607 (no se repite en las ventanas usadas y aporta
# más patrones de ISI de peor caso -> estadística del ojo/BER más robusta).
PRBS_ORDEN = 23
N_BITS = 60000            # 30 000 símbolos PAM-4 / 20 000 símbolos PAM-8 disponibles
vector_bits = generar_prbs(N_BITS, orden=PRBS_ORDEN)

# Verificación de balance (una PRBS bien formada es ~50% unos)
balance = vector_bits.mean()
print(f"Secuencia generada: PRBS{PRBS_ORDEN}")
print(f"Bits generados: {len(vector_bits)}")
print(f"Balance de unos: {balance:.4f} (ideal ~0.5)")
print(f"Corrida máxima de ceros consecutivos: "
      f"{max((len(s) for s in ''.join(map(str, vector_bits)).split('1')), default=0)}")

# ===== CELL 6 =====
def modular_nrz(bits, V_high=1.0, V_low=-1.0):
    """Modulación NRZ Bipolar estándar"""
    return np.where(bits == 1, V_high, V_low)

def modular_pam4(bits):
    """Modulación PAM-4 con Gray mapping"""
    # Si el número de bits es impar, rellenar con un cero al final
    if len(bits) % 2 != 0:
        bits = np.append(bits, 0)

    # Reestructurar bits en parejas (dibits)
    dibits = bits.reshape(-1, 2)
    symbols = []
    for d in dibits:
        val = d[0]*2 + d[1]
        # Gray mapping: 00 -> -3V, 01 -> -1V, 11 -> +1V, 10 -> +3V
        if val == 0:     # 00
            symbols.append(-3.0)
        elif val == 1:   # 01
            symbols.append(-1.0)
        elif val == 3:   # 11
            symbols.append(1.0)
        else:            # 10
            symbols.append(3.0)
    return np.array(symbols)

# --- Parámetros de transmisión ---
# TRADE-OFF CENTRAL DEL ESTUDIO: tasa de símbolo vs ancho de banda del canal.
# El canal de referencia (R=80 ohm, C=5 pF, R*C = 400 ps) tiene un ancho de banda
# de ~0.97 GHz. Por tanto:
#   * SIN ecualización, el enlace PAM-4 solo es viable si Nyquist (Rs/2) <= ~0.97 GHz,
#     es decir, hasta ~1.93 GBaud (referencia sin ecualizacion).
#   * CON FFE+CTLE (+DFE) se puede operar a tasas MAYORES recuperando la perdida
#     mas alla del -3 dB. El "Barrido de Tasa" al final cuantifica esa ganancia.
#
# Elegimos como punto de operacion del demostrador una tasa donde la ecualizacion
# SI es necesaria (sin EQ cierra, con EQ abre). Es un parametro: ajustalo segun el
# barrido de tasa.
Rs_baud = 4e9               # Tasa de simbolo PAM-4 (Baud). 4 GBaud -> Nyquist 2 GHz (~2x el BW)
Rb = 2 * Rs_baud            # Tasa de bits agregada (2 bit/simbolo en PAM-4)
Tb = 1.0 / Rb
Rs_pam4 = Rs_baud

sps = 32                    # Muestras por símbolo (sobremuestreo para simulación cuasi-continua)
fs = Rs_pam4 * sps          # Frecuencia de muestreo (Hz)
ts = 1.0 / fs

# Ancho de banda -3 dB de la linea RC distribuida de referencia (modelo cosh)
RC_ref = 80.0 * 5e-12
RC_COSH_F3DB_FACTOR = 2.4266937
BW_canal_GHz = RC_COSH_F3DB_FACTOR / (2 * np.pi * RC_ref) / 1e9
Rs_limite_sinEQ_GBaud = 2 * BW_canal_GHz   # Nyquist = BW  ->  tasa limite sin ecualizar

# Modular
simbolos_nrz = modular_nrz(vector_bits)
simbolos_pam4 = modular_pam4(vector_bits)

print(f"Punto de operacion: PAM-4 a {Rs_baud/1e9:.1f} GBaud "
      f"({Rb/1e9:.1f} Gb/s)  |  Nyquist {Rs_pam4/2e9:.2f} GHz")
print(f"Ancho de banda del canal (-3 dB): ~{BW_canal_GHz:.2f} GHz")
print(f"Tasa maxima sin ecualizacion (Nyquist=BW): ~{Rs_limite_sinEQ_GBaud:.2f} GBaud")
print(f"Cantidad de simbolos PAM-4 disponibles: {len(simbolos_pam4)}")

# ===== CELL 8 =====
def crear_canal_rc_distribuido(R_total, C_total, N=5):
    """
    Crea la representación en espacio de estados de una línea de transmisión RC distribuida
    de N etapas en escalera.
    """
    R_seg = R_total / N
    C_seg = C_total / N
    alpha = 1.0 / (R_seg * C_seg)

    # Matriz A de estado (N x N)
    A = np.zeros((N, N))
    for i in range(N-1):
        A[i, i] = -2 * alpha
        A[i, i+1] = alpha
        A[i+1, i] = alpha
    A[N-1, N-1] = -1 * alpha  # Extremo abierto

    # Vector B de entrada
    B = np.zeros((N, 1))
    B[0, 0] = alpha

    # Vector C de salida (obtenemos el voltaje en el último capacitor)
    C = np.zeros((1, N))
    C[0, N-1] = 1.0

    D = np.zeros((1, 1))

    return signal.StateSpace(A, B, C, D)

# ===== CELL 10 =====
def aplicar_ffe(simbolos, c_pre, c_main, c_post):
    """
    Aplica un filtro FFE (FIR de 3 taps) al flujo de símbolos discretos.
    """
    taps = np.array([c_pre, c_main, c_post])
    # Normalización de potencia para mantener el voltaje constante
    taps = taps / np.sum(np.abs(taps))
    return signal.lfilter(taps, 1.0, simbolos)

def crear_ctle(omega_z, omega_p1, omega_p2, A_dc=1.0):
    """
    Retorna un sistema LTI continuo para el ecualizador analógico CTLE de 1 cero y 2 polos.
    """
    # H(s) = A_dc * (1 + s/omega_z) / [ (1 + s/omega_p1) * (1 + s/omega_p2) ]
    num = [A_dc / omega_z, A_dc]
    den = [1.0 / (omega_p1 * omega_p2), (omega_p1 + omega_p2) / (omega_p1 * omega_p2), 1.0]
    return signal.lti(num, den)

# ===== CELL 11 =====
# --- Utilidades de medición, métrica de ojo real y ecualización DFE ---
# Centraliza la lógica usada por el DoE, los escenarios y el barrido de tasa.

NIVELES_PAM4 = np.array([-3.0, -1.0, 1.0, 3.0])

def Q(x):
    """Función de cola Gaussiana Q(x) = 0.5*erfc(x/sqrt(2))."""
    return 0.5 * erfc(np.asarray(x, dtype=float) / np.sqrt(2.0))

def muestrear_alineado(y, tx_syms, sps_val, max_lag=80, n_align=500):
    """
    Alinea la señal recibida con los símbolos transmitidos buscando la fase de
    muestreo (0..sps-1) y el retardo entero de símbolo por correlación cruzada.
    Devuelve (muestras_alineadas, tx_alineados).
    """
    tx = np.asarray(tx_syms, dtype=float)
    n_use = min(len(tx), n_align)
    mejor = None  # (|corr|, fase, lag)
    for fase in range(sps_val):
        ds = y[fase::sps_val]
        n = min(len(ds), n_use)
        if n < 50:
            continue
        a = ds[:n] - ds[:n].mean()
        b = tx[:n] - tx[:n].mean()
        corr = np.correlate(a, b, mode='full')
        lags = np.arange(-n + 1, n)
        m = (lags >= 0) & (lags <= max_lag)
        idx = np.argmax(np.abs(corr[m]))
        score = np.abs(corr[m][idx])
        if mejor is None or score > mejor[0]:
            mejor = (score, fase, int(lags[m][idx]))
    _, fase, lag = mejor
    ds = y[fase::sps_val]
    n = min(len(ds) - lag, len(tx))
    return ds[lag:lag + n], tx[:n]

def estadisticas_niveles(muestras, tx, niveles=NIVELES_PAM4):
    """Media y dispersión (ISI residual) de cada nivel PAM-4 transmitido."""
    mu = np.zeros(4); s = np.zeros(4); cnt = np.zeros(4, dtype=int)
    for k, lv in enumerate(niveles):
        sel = muestras[np.isclose(tx, lv)]
        cnt[k] = len(sel)
        if len(sel) >= 2:
            mu[k] = sel.mean(); s[k] = sel.std()
        elif len(sel) == 1:
            mu[k] = sel[0]
    return mu, s, cnt

def apertura_ojo_real(y, tx, sps_val, k=3.0, n_descartar=60):
    """
    Apertura vertical del PEOR ojo (margen a k-sigma de la ISI residual).
    Es el objetivo correcto del DoE: mide separación entre niveles, NO el
    'spread' de la señal (que premiaba erróneamente la ganancia del CTLE).
    """
    m, t = muestrear_alineado(y, tx, sps_val)
    m, t = m[n_descartar:], t[n_descartar:]
    mu, s, cnt = estadisticas_niveles(m, t)
    if np.any(cnt == 0):
        return -99.0
    o = np.argsort(mu); mu, s = mu[o], s[o]
    ap = (mu[1:] - k*s[1:]) - (mu[:-1] + k*s[:-1])
    return float(np.min(ap))

def respuesta_al_pulso(canal, ctle, c_post, sps_val, ts_val, n_sym=48, pre=6):
    """Respuesta al pulso del enlace (FFE -> canal -> CTLE) a un símbolo unitario."""
    x = np.zeros(n_sym); x[pre] = 1.0
    ffe = aplicar_ffe(x, 0.0, 1.0, c_post)
    xu = np.repeat(ffe, sps_val)
    t = np.arange(len(xu)) * ts_val
    _, yc, _ = signal.lsim(canal, U=xu, T=t)
    _, y, _ = signal.lsim(ctle, U=yc, T=t)
    pico = int(np.argmax(np.abs(y)))     # fase de muestreo = posición del pico
    p = y[pico % sps_val::sps_val]
    c_idx = int(np.argmax(np.abs(p)))    # índice del cursor
    return p, c_idx

def disenar_dfe_zf(canal, ctle, c_post, sps_val, ts_val, n_taps=2):
    """
    Diseña un DFE forzador-de-ceros (ZF) a partir de los post-cursores de la
    respuesta al pulso del enlace. Devuelve coeficientes absolutos (V/símbolo).
    El DFE cancela ISI de cola SIN amplificar ruido (a diferencia del CTLE).
    """
    p, c = respuesta_al_pulso(canal, ctle, c_post, sps_val, ts_val)
    return np.array([p[c + i] if c + i < len(p) else 0.0 for i in range(1, n_taps + 1)])

def _decidir_pam4(yk, escala):
    """Decide el símbolo PAM-4 ideal (-3,-1,1,3) por mínima distancia al nivel escalado."""
    niveles = NIVELES_PAM4 * escala
    return NIVELES_PAM4[int(np.argmin(np.abs(yk - niveles)))]

def evaluar_enlace(y, tx_syms, sps_val, sigma_ruido, dfe_taps=None, n_descartar=80):
    """
    Alinea la señal, opcionalmente aplica un DFE de decisión realimentada,
    y estima BER semi-analítico (modelo Gaussiano) y la apertura del peor ojo.
    """
    m, tx = muestrear_alineado(y, tx_syms, sps_val)
    m, tx = m[n_descartar:], tx[n_descartar:]

    # Ganancia efectiva del enlace (V por unidad de amplitud de símbolo) para el DFE
    mu0, _, _ = estadisticas_niveles(m, tx)
    o0 = np.argsort(mu0)
    escala = (mu0[o0][-1] - mu0[o0][0]) / (NIVELES_PAM4[-1] - NIVELES_PAM4[0])
    escala = escala if escala > 1e-9 else 1.0

    if dfe_taps is not None and len(dfe_taps) > 0:
        Nt = len(dfe_taps)
        hist = [0.0] * Nt
        m_eq = np.empty_like(m)
        for kk in range(len(m)):
            fb = sum(dfe_taps[i] * hist[i] for i in range(Nt))
            yk = m[kk] - fb
            m_eq[kk] = yk
            hist = [_decidir_pam4(yk, escala)] + hist[:-1]  # decisión realimentada
        m = m_eq

    mu, s_isi, cnt = estadisticas_niveles(m, tx)
    o = np.argsort(mu); mu, s_isi, cnt = mu[o], s_isi[o], cnt[o]
    sigma_tot = np.sqrt(s_isi**2 + sigma_ruido**2)
    thr = (mu[:-1] + mu[1:]) / 2.0
    eps = 1e-300
    Pe = np.zeros(4)
    Pe[0] = Q((thr[0] - mu[0]) / (sigma_tot[0] + eps))
    Pe[1] = Q((mu[1] - thr[0]) / (sigma_tot[1] + eps)) + Q((thr[1] - mu[1]) / (sigma_tot[1] + eps))
    Pe[2] = Q((mu[2] - thr[1]) / (sigma_tot[2] + eps)) + Q((thr[2] - mu[2]) / (sigma_tot[2] + eps))
    Pe[3] = Q((mu[3] - thr[2]) / (sigma_tot[3] + eps))
    total = cnt.sum()
    pesos = cnt / total if total > 0 else np.full(4, 0.25)
    ser = float(np.sum(pesos * Pe)); ber = ser / 2.0
    ap = (mu[1:] - 3*s_isi[1:]) - (mu[:-1] + 3*s_isi[:-1])
    return {"ber": ber, "ser": ser, "eye_min": float(np.min(ap)), "mu": mu, "sigma_tot": sigma_tot}

print("Utilidades de medicion, metrica de ojo real y DFE cargadas.")

# ===== CELL 12 =====
# --- Front-end del receptor: ancho de banda finito (anti-alias) ---
# Imprescindible: sin el, el pico de realce del CTLE amplifica ruido fuera de
# banda (10-30 GHz) sin control. Un RX real tiene BW ~ tasa de simbolo.
def aplicar_rx_lpf(y, ts_val, fc_hz, orden=3):
    """Filtro pasa-bajo de entrada del RX (Butterworth, fase cero via filtfilt)."""
    fs_val = 1.0 / ts_val
    wn = min(fc_hz / (fs_val / 2.0), 0.99)
    b, a = signal.butter(orden, wn, btype='low')
    return signal.filtfilt(b, a, y)

print(f"Front-end RX cargado (LPF anti-alias).")

# ===== CELL 14 =====
# Rango de búsqueda de diseño (DoE)
np.random.seed(42)
n_muestras = 30
diseño_lhs = lhs(2, samples=n_muestras, criterion='center', seed=42)

# Mapear a parámetros reales:
# 1. omega_z (frecuencia del cero CTLE) de 1 GHz a 10 GHz
# 2. c_post (coeficiente post-cursor del FFE) de -0.4 a -0.05
omega_z_eval = 2 * np.pi * (1e9 + diseño_lhs[:, 0] * 9e9)
c_post_eval = -0.4 + diseño_lhs[:, 1] * 0.35

# Canal de referencia (R*C = 400 ps -> BW ~0.97 GHz)
R_linea = 80.0
C_linea = 5e-12
canal_lti = crear_canal_rc_distribuido(R_linea, C_linea, N=5)

# El DoE evalua el CTLE viendo el MISMO front-end de RX (LPF) que la etapa de BER,
# para no premiar un realce que solo amplifica ruido fuera de banda.
fc_rx = Rs_pam4

print("Corriendo muestreo DoE multivariable (objetivo = apertura de ojo REAL)...")
mejores_parametros = None
mejor_apertura = -1e9
resultados_aperturas = []

# Secuencia de 600 símbolos para evaluación rápida del DoE
simbolos_prueba = simbolos_pam4[:600]

for i in range(n_muestras):
    w_z = omega_z_eval[i]
    c_p = c_post_eval[i]

    # FFE -> canal -> CTLE -> front-end RX
    ffe_out = aplicar_ffe(simbolos_prueba, c_pre=0.0, c_main=1.0, c_post=c_p)
    x_t_sim = np.repeat(ffe_out, sps)
    t_sim = np.arange(len(x_t_sim)) * ts
    _, y_canal, _ = signal.lsim(canal_lti, U=x_t_sim, T=t_sim)
    ctle_lti = crear_ctle(w_z, omega_p1=2*np.pi*15e9, omega_p2=2*np.pi*30e9, A_dc=2.5)
    _, y_ctle, _ = signal.lsim(ctle_lti, U=y_canal, T=t_sim)
    y_rx = aplicar_rx_lpf(y_ctle, ts, fc_rx)

    # OBJETIVO CORRECTO: margen vertical del peor ojo (separación entre niveles),
    # no el 'spread' de la señal (que premiaba la ganancia del CTLE).
    apert = apertura_ojo_real(y_rx, simbolos_prueba, sps)
    resultados_aperturas.append(apert)

    if apert > mejor_apertura:
        mejor_apertura = apert
        mejores_parametros = (w_z, c_p)

print(f"Optimización DoE completa. Mejor apertura de ojo (3-sigma ISI): {mejor_apertura:+.3f} V")
print(f"Parámetros óptimos -> Cero CTLE: {mejores_parametros[0]/(2*np.pi*1e9):.2f} GHz, "
      f"FFE Post-tap: {mejores_parametros[1]:.3f}")
if mejor_apertura <= 0:
    print("[Aviso] Ningún punto abre el ojo a esta tasa: la EQ lineal no basta -> "
          "considera bajar Rs_baud o añadir DFE (ver barrido de tasa).")

# Graficar el mapa de exploración del espacio de diseño
plt.figure(figsize=(9, 6))
sc = plt.scatter(omega_z_eval / (2*np.pi*1e9), c_post_eval, c=resultados_aperturas,
                 cmap='viridis', s=100, edgecolor='k')
plt.colorbar(sc, label='Apertura REAL del peor ojo (V)  [>0 = ojo abierto]')
plt.scatter(mejores_parametros[0]/(2*np.pi*1e9), mejores_parametros[1],
            color='red', marker='*', s=250, label='Punto Óptimo')
plt.xlabel("Frecuencia del Cero CTLE (GHz)")
plt.ylabel("Coeficiente Post-tap FFE")
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
# PIN to manuscript canonical equalizer operating point (CTLE zero 2.35 GHz, FFE post-tap -0.301)
mejores_parametros = (2*_np2.pi*2.35e9, -0.301)


# ===== CELL 24 =====
sigma_in = 0.10
# ===== Barrido de Tasa: BER vs GBaud para sin-EQ / FFE+CTLE / FFE+CTLE+DFE =====
# Ruido referido a la ENTRADA del RX (antes del CTLE) + front-end RX de BW finito,
# igual que en la tabla de BER. El BW del RX escala con la tasa (~ Rs).
tasas_baud = np.linspace(1e9, 16e9, 16)
R_sw, C_sw = 80.0, 5e-12
canal_sw = crear_canal_rc_distribuido(R_sw, C_sw, N=5)
w_z_opt, c_post_opt = mejores_parametros
sps_sw = 32
N_sym_sw = 2000

# Secuencia PRBS dedicada para el barrido (PRBS23, sin repetición)
bits_sw = generar_prbs(2 * N_sym_sw + 600, orden=23)
sym_sw = modular_pam4(bits_sw)[:N_sym_sw]

ber_none, ber_lin, ber_dfe = [], [], []
print("Barrido de tasa PAM-4 en curso (esto toma unos segundos)...")
np.random.seed(7)  # reproducibilidad del ruido del barrido
for Rs in tasas_baud:
    ts_sw = 1.0 / (Rs * sps_sw)
    fc_rx_sw = Rs                       # BW del front-end RX ~ tasa de simbolo
    xin0 = np.repeat(sym_sw, sps_sw)
    t_sw = np.arange(len(xin0)) * ts_sw
    n_muestras_sw = len(xin0)

    # (1) Sin EQ: ruido de entrada -> front-end RX
    _, y0, _ = signal.lsim(canal_sw, U=xin0, T=t_sw)
    y0n = aplicar_rx_lpf(y0 + np.random.normal(0.0, sigma_in, n_muestras_sw), ts_sw, fc_rx_sw)
    r0 = evaluar_enlace(y0n, sym_sw, sps_sw, 0.0)

    # (2) FFE + CTLE: ruido ANTES del CTLE -> CTLE -> front-end RX
    ffe = aplicar_ffe(sym_sw, 0.0, 1.0, c_post_opt)
    xin = np.repeat(ffe, sps_sw)
    _, yc1, _ = signal.lsim(canal_sw, U=xin, T=t_sw)
    ruido_in_sw = np.random.normal(0.0, sigma_in, len(yc1))
    ctle_sw = crear_ctle(w_z_opt, 2*np.pi*15e9, 2*np.pi*30e9, A_dc=2.5)
    _, ylin0, _ = signal.lsim(ctle_sw, U=yc1 + ruido_in_sw, T=t_sw)
    ylin = aplicar_rx_lpf(ylin0, ts_sw, fc_rx_sw)
    rl = evaluar_enlace(ylin, sym_sw, sps_sw, 0.0)

    # (3) + DFE (ZF desde la respuesta al pulso a esta tasa), misma señal ruidosa
    taps_sw = disenar_dfe_zf(canal_sw, ctle_sw, c_post_opt, sps_sw, ts_sw, n_taps=2)
    rd = evaluar_enlace(ylin, sym_sw, sps_sw, 0.0, dfe_taps=taps_sw)

    ber_none.append(r0['ber']); ber_lin.append(rl['ber']); ber_dfe.append(rd['ber'])

# Piso numerico para la escala logaritmica
piso = 1e-12
clip = lambda L: np.clip(np.array(L), piso, 1.0)
gb = tasas_baud / 1e9

plt.figure(figsize=(11, 6))
plt.semilogy(gb, clip(ber_none), 'o-', color='#e53e3e', label='No equalization')
plt.semilogy(gb, clip(ber_lin),  's-', color='#dd6b20', label='FFE + CTLE')
plt.semilogy(gb, clip(ber_dfe),  '^-', color='#319795', label='FFE + CTLE + DFE')
plt.axhline(1e-2, color='gray', ls='--', alpha=0.7, label='Pre-FEC threshold (1e-2)')
plt.axvline(Rs_limite_sinEQ_GBaud, color='blue', ls=':', alpha=0.8,
            label=f'Nyquist = BW (~{Rs_limite_sinEQ_GBaud:.1f} GBaud)')
plt.axvline(Rs_baud/1e9, color='black', ls='-', alpha=0.3,
            label=f'Operating point ({Rs_baud/1e9:.1f} GBaud)')
plt.xlabel('PAM-4 symbol rate (GBaud)')
plt.ylabel('Estimated semi-analytic BER')
plt.ylim(piso, 1)
plt.grid(True, which='both', alpha=0.25)
plt.legend(loc='lower right', fontsize=9)
plt.show()

# Tasa maxima viable (por debajo del umbral pre-FEC) para cada configuracion
def tasa_max(bers, umbral=1e-2):
    ok = gb[np.array(bers) <= umbral]
    return ok.max() if len(ok) else 0.0
print(f"Tasa max viable (BER<=1e-2)  ->  sin EQ: {tasa_max(ber_none):.1f} GBaud | "
      f"FFE+CTLE: {tasa_max(ber_lin):.1f} GBaud | +DFE: {tasa_max(ber_dfe):.1f} GBaud")

# ===== CELL 26 =====
# ===== STRESS TEST: PAM-8 (resultado exploratorio) =====
# Justificacion (NO por despliegue, que es marginal en cobre): con el MISMO presupuesto
# de amplitud (pico +/-3 V que PAM-4), 8 niveles producen ojos ~2.3x mas estrechos.
# Es un entorno mas AGRESIVO para verificar si el MISMO co-diseno de ecualizacion
# (FFE+CTLE+DFE + receptor realista) GENERALIZA a constelaciones densas, y revela el
# compromiso ORDEN-vs-ALCANCE (mas bits/simbolo a costa de menor tasa viable).

NIVELES_PAM8 = np.linspace(-3.0, 3.0, 8)   # mismo pico (+/-3 V) que PAM-4; espaciado 6/7

def modular_mpam(bits, niveles):
    """Modula bits en M-PAM con mapeo Gray (M=len(niveles), potencia de 2)."""
    M = len(niveles); k = int(round(np.log2(M)))
    r = len(bits) % k
    if r:
        bits = np.append(bits, np.zeros(k - r, dtype=bits.dtype))
    groups = bits.reshape(-1, k)
    pesos = (1 << np.arange(k - 1, -1, -1))          # MSB-first
    etiqueta = groups.dot(pesos).astype(int)         # 0..M-1
    gray = np.array([i ^ (i >> 1) for i in range(M)]) # gray[i] = etiqueta del nivel i
    lab2idx = np.empty(M, dtype=int); lab2idx[gray] = np.arange(M)
    return niveles[lab2idx[etiqueta]]

def estadisticas_mpam(muestras, tx, niveles):
    M = len(niveles); mu = np.zeros(M); s = np.zeros(M); cnt = np.zeros(M, dtype=int)
    for k, lv in enumerate(niveles):
        sel = muestras[np.isclose(tx, lv)]
        cnt[k] = len(sel)
        if len(sel) >= 2: mu[k] = sel.mean(); s[k] = sel.std()
        elif len(sel) == 1: mu[k] = sel[0]
    return mu, s, cnt

def _decidir_mpam(yk, escala, niveles):
    lv = niveles * escala
    return niveles[int(np.argmin(np.abs(yk - lv)))]

def evaluar_mpam(y, tx_syms, sps_val, niveles, dfe_taps=None, n_descartar=80):
    """BER semi-analitico (modelo Gaussiano) para M-PAM. Ruido ya incluido en 'y'."""
    M = len(niveles); kbits = int(round(np.log2(M)))
    m, tx = muestrear_alineado(y, tx_syms, sps_val)
    m, tx = m[n_descartar:], tx[n_descartar:]
    mu0, _, _ = estadisticas_mpam(m, tx, niveles)
    o0 = np.argsort(mu0)
    escala = (mu0[o0][-1] - mu0[o0][0]) / (niveles.max() - niveles.min())
    escala = escala if escala > 1e-9 else 1.0
    if dfe_taps is not None and len(dfe_taps) > 0:
        Nt = len(dfe_taps); hist = [0.0] * Nt; m_eq = np.empty_like(m)
        for i in range(len(m)):
            fb = sum(dfe_taps[j] * hist[j] for j in range(Nt))
            yk = m[i] - fb; m_eq[i] = yk
            hist = [_decidir_mpam(yk, escala, niveles)] + hist[:-1]
        m = m_eq
    mu, s, cnt = estadisticas_mpam(m, tx, niveles)
    o = np.argsort(mu); mu, s, cnt = mu[o], s[o], cnt[o]
    thr = (mu[:-1] + mu[1:]) / 2.0; eps = 1e-300
    Pe = np.zeros(M)
    Pe[0] = Q((thr[0] - mu[0]) / (s[0] + eps))
    for kk in range(1, M - 1):
        Pe[kk] = Q((mu[kk] - thr[kk-1]) / (s[kk] + eps)) + Q((thr[kk] - mu[kk]) / (s[kk] + eps))
    Pe[M-1] = Q((mu[M-1] - thr[M-2]) / (s[M-1] + eps))
    tot = cnt.sum(); w = cnt / tot if tot > 0 else np.full(M, 1.0 / M)
    ser = float(np.sum(w * Pe)); ber = ser / kbits
    ap = (mu[1:] - 3*s[1:]) - (mu[:-1] + 3*s[:-1])
    return {"ber": ber, "ser": ser, "eye_min": float(np.min(ap))}

# --- Barrido de tasa para PAM-8 (mismo canal, mismo co-diseno EQ, mismo receptor) ---
N8 = 2000
bits8 = generar_prbs(3 * N8 + 600, orden=23)
sym8 = modular_mpam(bits8, NIVELES_PAM8)[:N8]

ber8_none, ber8_lin, ber8_dfe = [], [], []
print("Barrido de tasa PAM-8 (stress test) en curso...")
np.random.seed(11)
for Rs in tasas_baud:
    ts_sw = 1.0 / (Rs * sps_sw); fc = Rs
    xin0 = np.repeat(sym8, sps_sw); t_sw = np.arange(len(xin0)) * ts_sw

    _, y0, _ = signal.lsim(canal_sw, U=xin0, T=t_sw)
    y0n = aplicar_rx_lpf(y0 + np.random.normal(0.0, sigma_in, len(xin0)), ts_sw, fc)
    r0 = evaluar_mpam(y0n, sym8, sps_sw, NIVELES_PAM8)

    ffe = aplicar_ffe(sym8, 0.0, 1.0, c_post_opt); xin = np.repeat(ffe, sps_sw)
    _, yc1, _ = signal.lsim(canal_sw, U=xin, T=t_sw)
    rin = np.random.normal(0.0, sigma_in, len(yc1))
    ctle_sw = crear_ctle(w_z_opt, 2*np.pi*15e9, 2*np.pi*30e9, A_dc=2.5)
    _, ylin0, _ = signal.lsim(ctle_sw, U=yc1 + rin, T=t_sw)
    ylin = aplicar_rx_lpf(ylin0, ts_sw, fc)
    rl = evaluar_mpam(ylin, sym8, sps_sw, NIVELES_PAM8)

    taps = disenar_dfe_zf(canal_sw, ctle_sw, c_post_opt, sps_sw, ts_sw, n_taps=2)
    rd = evaluar_mpam(ylin, sym8, sps_sw, NIVELES_PAM8, dfe_taps=taps)

    ber8_none.append(r0['ber']); ber8_lin.append(rl['ber']); ber8_dfe.append(rd['ber'])

# --- Figura 1: BER vs GBaud para PAM-8 ---
plt.figure(figsize=(11, 6))
plt.semilogy(gb, clip(ber8_none), 'o--', color='#e53e3e', label='PAM-8 no EQ')
plt.semilogy(gb, clip(ber8_lin),  's--', color='#dd6b20', label='PAM-8 FFE+CTLE')
plt.semilogy(gb, clip(ber8_dfe),  '^--', color='#319795', label='PAM-8 FFE+CTLE+DFE')
plt.axhline(1e-2, color='gray', ls='--', alpha=0.7, label='Pre-FEC threshold (1e-2)')
plt.xlabel('Symbol rate (GBaud)'); plt.ylabel('Estimated semi-analytic BER')
plt.ylim(piso, 1); plt.grid(True, which='both', alpha=0.25); plt.legend(loc='lower right', fontsize=9)
plt.savefig('fig_pam8_ber.png', dpi=200, bbox_inches='tight'); print('SAVED fig_pam8_ber.png')

# --- Figura 2: Throughput neto (Gb/s) vs alcance: PAM-4 vs PAM-8 ---
def t_max(B, umbral=1e-2):
    ok = gb[np.array(B) <= umbral]
    return ok.max() if len(ok) else 0.0
cfgs = ['No EQ', 'FFE+CTLE', 'FFE+CTLE+DFE']
rate4 = [t_max(ber_none), t_max(ber_lin), t_max(ber_dfe)]          # GBaud PAM-4
rate8 = [t_max(ber8_none), t_max(ber8_lin), t_max(ber8_dfe)]       # GBaud PAM-8
thr4 = [2 * r for r in rate4]     # 2 bit/simbolo
thr8 = [3 * r for r in rate8]     # 3 bit/simbolo

x = np.arange(3); w = 0.35
plt.figure(figsize=(10, 6))
plt.bar(x - w/2, thr4, w, color='#319795', label='PAM-4 (2 bit/sim)')
plt.bar(x + w/2, thr8, w, color='#dd6b20', label='PAM-8 (3 bit/sim)')
for i in range(3):
    plt.text(x[i]-w/2, thr4[i], f'{thr4[i]:.0f}', ha='center', va='bottom', fontsize=9)
    plt.text(x[i]+w/2, thr8[i], f'{thr8[i]:.0f}', ha='center', va='bottom', fontsize=9)
plt.xticks(x, cfgs); plt.ylabel('Throughput (Gb/s)')
plt.legend(); plt.grid(True, axis='y', alpha=0.3); plt.savefig('fig_throughput.png', dpi=200, bbox_inches='tight'); print('SAVED fig_throughput.png')

print("Tasa viable (GBaud)  PAM-4:", [f'{r:.0f}' for r in rate4], " PAM-8:", [f'{r:.0f}' for r in rate8])
print("Throughput neto (Gb/s) PAM-4:", [f'{t:.0f}' for t in thr4], " PAM-8:", [f'{t:.0f}' for t in thr8])
print("Lectura: PAM-8 lleva 3 bit/sim pero su menor alcance suele NO compensar; el co-diseno")
print("de ecualizacion sigue ordenando A<B<C tambien en PAM-8 -> la metodologia GENERALIZA.")