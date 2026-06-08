"""
simulacion_ojo.py — Ojo y señal transmitida sobre el canal FÍSICO dispersivo.

Replica la cadena de `simulacion_canal_rcEjecutadoActual.ipynb`
(PRBS23 -> PAM-4 -> FFE -> canal -> CTLE -> LPF RX -> DFE -> ojo / BER), pero
sustituye el canal RC en espacio de estados por el canal de línea de transmisión
con pérdidas dependientes de frecuencia (efecto pelicular + dieléctrico) de
`modelo_rlgc.py`. El canal se aplica en el dominio de la frecuencia: y = IFFT( H(f)·FFT(x) ).

Mantiene IDÉNTICOS el resto de bloques y el co-diseño DoE fijado en checkpoint.md:
    omega_z = 2π·2.35 GHz,  c_post = -0.301,  DFE ZF de 2 taps,
    CTLE: 1 cero + 2 polos (15 / 30 GHz), A_dc = 2.5,
    ruido sigma_in = 0.10 V referido a la ENTRADA del RX + LPF de BW finito.

Entregables (carpeta ./figuras/):
  - fig_ojo_senal_transmitida.png : señal transmitida vs recibida (transitorio)
  - fig_ojo_escenarios.png        : ojos A / B / C + niveles muestreados C vs D(+DFE)
  - fig_ojo_comparativa_canal.png : ojo del escenario C  -> RC vs FR-4 vs Megtron 6
  - tabla_ojo_ber.csv             : BER y apertura de ojo por canal y escenario

Uso:
    python simulacion_ojo.py
"""

import csv
import os

import numpy as np
import scipy.signal as signal
from scipy.special import erfc
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from parametros import SUBSTRATOS, FR4, MEGTRON6, OPERACION, INCH
from modelo_rlgc import LineaTransmision
from modelo_rc import CanalRC

DIR_FIG = os.path.join(os.path.dirname(__file__), "figuras")
DPI = 300

# --------------------------------------------------------------------------- #
# Punto de operación y co-diseño (tomados del notebook / checkpoint)
# --------------------------------------------------------------------------- #
RS_BAUD = OPERACION.baud_op          # 4 GBaud PAM-4
SPS = 32                             # muestras por símbolo
FS = RS_BAUD * SPS
TS = 1.0 / FS

SIGMA_IN = 0.10                      # V RMS, referido a la entrada del RX
FC_RX = RS_BAUD                      # BW del front-end RX (~ tasa de símbolo)

W_Z = 2 * np.pi * 2.35e9             # cero del CTLE (DoE fijado)
C_POST = -0.301                      # post-cursor FFE (DoE fijado)
CTLE_P1, CTLE_P2, CTLE_ADC = 2 * np.pi * 15e9, 2 * np.pi * 30e9, 2.5
DFE_NTAPS = 2

CANT_SIMBOLOS = 4000                 # ventana de símbolos para ojo/BER
NIVELES_PAM4 = np.array([-3.0, -1.0, 1.0, 3.0])


# =========================================================================== #
# Bloques portados del notebook (idénticos)
# =========================================================================== #
def generar_prbs(n_bits, orden=23, semilla=None):
    """LFSR de Fibonacci con polinomios primitivos estándar (PRBS7..31)."""
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
    """PAM-4 con Gray mapping: 00->-3, 01->-1, 11->+1, 10->+3."""
    if len(bits) % 2 != 0:
        bits = np.append(bits, 0)
    dibits = bits.reshape(-1, 2)
    out = []
    for d in dibits:
        val = d[0] * 2 + d[1]
        out.append({0: -3.0, 1: -1.0, 3: 1.0, 2: 3.0}[val])
    return np.array(out)


def aplicar_ffe(simbolos, c_pre, c_main, c_post):
    """FFE FIR de 3 taps con normalización de potencia."""
    taps = np.array([c_pre, c_main, c_post])
    taps = taps / np.sum(np.abs(taps))
    return signal.lfilter(taps, 1.0, simbolos)


def crear_ctle(omega_z, omega_p1, omega_p2, A_dc=1.0):
    """CTLE analógico: 1 cero, 2 polos."""
    num = [A_dc / omega_z, A_dc]
    den = [1.0 / (omega_p1 * omega_p2), (omega_p1 + omega_p2) / (omega_p1 * omega_p2), 1.0]
    return signal.lti(num, den)


def aplicar_rx_lpf(y, ts_val, fc_hz, orden=3):
    """LPF anti-alias del front-end del RX (Butterworth, fase cero)."""
    fs_val = 1.0 / ts_val
    wn = min(fc_hz / (fs_val / 2.0), 0.99)
    b, a = signal.butter(orden, wn, btype="low")
    return signal.filtfilt(b, a, y)


def Q(x):
    return 0.5 * erfc(np.asarray(x, dtype=float) / np.sqrt(2.0))


def _find_align(y, tx_syms, sps_val, max_lag=80, n_align=500):
    """Devuelve (fase, lag) óptimos de muestreo por correlación cruzada."""
    tx = np.asarray(tx_syms, dtype=float)
    n_use = min(len(tx), n_align)
    mejor = None
    for fase in range(sps_val):
        ds = y[fase::sps_val]
        n = min(len(ds), n_use)
        if n < 50:
            continue
        a = ds[:n] - ds[:n].mean()
        b = tx[:n] - tx[:n].mean()
        corr = np.correlate(a, b, mode="full")
        lags = np.arange(-n + 1, n)
        m = (lags >= 0) & (lags <= max_lag)
        idx = np.argmax(np.abs(corr[m]))
        score = np.abs(corr[m][idx])
        if mejor is None or score > mejor[0]:
            mejor = (score, fase, int(lags[m][idx]))
    return mejor[1], mejor[2]


def muestrear_alineado(y, tx_syms, sps_val):
    """Muestra 1 valor por símbolo en la fase/lag óptimos."""
    fase, lag = _find_align(y, tx_syms, sps_val)
    tx = np.asarray(tx_syms, dtype=float)
    ds = y[fase::sps_val]
    n = min(len(ds) - lag, len(tx))
    return ds[lag:lag + n], tx[:n]


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


def evaluar_enlace(y, tx_syms, sps_val, sigma_ruido, dfe_taps=None,
                   n_descartar=80, devolver_muestras=False):
    """Alinea, aplica DFE opcional, estima BER semi-analítico y apertura de ojo."""
    m, tx = muestrear_alineado(y, tx_syms, sps_val)
    m, tx = m[n_descartar:], tx[n_descartar:]

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
            hist = [_decidir_pam4(yk, escala)] + hist[:-1]
        m = m_eq

    mu, s_isi, cnt = estadisticas_niveles(m, tx)
    o = np.argsort(mu); mu, s_isi, cnt = mu[o], s_isi[o], cnt[o]
    sigma_tot = np.sqrt(s_isi ** 2 + sigma_ruido ** 2)
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
    ap = (mu[1:] - 3 * s_isi[1:]) - (mu[:-1] + 3 * s_isi[:-1])
    res = {"ber": ber, "ser": ser, "eye_min": float(np.min(ap)), "mu": mu, "sigma": sigma_tot}
    if devolver_muestras:
        res["muestras"] = m
        res["tx"] = tx
    return res


# =========================================================================== #
# Bloque NUEVO: canal aplicado en frecuencia (RC cosh o línea dispersiva)
# =========================================================================== #
def aplicar_canal_freq(x, ts_val, H_func):
    """y = IFFT( H(f)·FFT(x) ).  H_func(f) -> respuesta compleja del canal."""
    n = len(x)
    f = np.fft.rfftfreq(n, ts_val)
    X = np.fft.rfft(x)
    return np.fft.irfft(X * H_func(f), n=n)


def disenar_dfe_zf(H_func, ctle, c_post, sps_val, ts_val, n_taps=2, n_sym=48, pre=6):
    """DFE forzador-de-ceros desde la respuesta al pulso (FFE -> canal -> CTLE)."""
    x = np.zeros(n_sym); x[pre] = 1.0
    ffe = aplicar_ffe(x, 0.0, 1.0, c_post)
    xu = np.repeat(ffe, sps_val)
    yc = aplicar_canal_freq(xu, ts_val, H_func)
    t = np.arange(len(xu)) * ts_val
    _, y, _ = signal.lsim(ctle, U=yc, T=t)
    pico = int(np.argmax(np.abs(y)))
    p = y[pico % sps_val::sps_val]
    c_idx = int(np.argmax(np.abs(p)))
    return np.array([p[c_idx + i] if c_idx + i < len(p) else 0.0
                     for i in range(1, n_taps + 1)])


# =========================================================================== #
# Generación de formas de onda por escenario (A/B/C) para un canal dado
# =========================================================================== #
def construir_escenarios(H_func, sym_ventana, rng):
    """Devuelve formas de onda limpias (para ojo) y ruidosas (para BER)."""
    t = np.arange(len(sym_ventana) * SPS) * TS

    # --- A: canal sin EQ ---
    x_a = np.repeat(sym_ventana, SPS)
    y_a = aplicar_canal_freq(x_a, TS, H_func)

    # --- B/C: FFE -> canal (mismo canal out); C añade CTLE ---
    ffe = aplicar_ffe(sym_ventana, 0.0, 1.0, C_POST)
    x_ffe = np.repeat(ffe, SPS)
    y_chan_ffe = aplicar_canal_freq(x_ffe, TS, H_func)          # = escenario B (limpio)
    ctle = crear_ctle(W_Z, CTLE_P1, CTLE_P2, CTLE_ADC)
    _, y_c, _ = signal.lsim(ctle, U=y_chan_ffe, T=t)            # escenario C (limpio)

    # --- Versiones ruidosas (ruido referido a la ENTRADA del RX) ---
    yA_n = aplicar_rx_lpf(y_a + rng.normal(0, SIGMA_IN, len(y_a)), TS, FC_RX)
    yB_n = aplicar_rx_lpf(y_chan_ffe + rng.normal(0, SIGMA_IN, len(y_chan_ffe)), TS, FC_RX)
    ruido_in = rng.normal(0, SIGMA_IN, len(y_chan_ffe))         # ANTES del CTLE
    _, yC_ctle, _ = signal.lsim(ctle, U=y_chan_ffe + ruido_in, T=t)
    yC_n = aplicar_rx_lpf(yC_ctle, TS, FC_RX)

    dfe_taps = disenar_dfe_zf(H_func, ctle, C_POST, SPS, TS, n_taps=DFE_NTAPS)

    return {
        "t": t, "x_a": x_a, "x_ffe": x_ffe,
        "limpio": {"A": y_a, "B": y_chan_ffe, "C": y_c},
        "ruidoso": {"A": yA_n, "B": yB_n, "C": yC_n},
        "dfe_taps": dfe_taps,
    }


# =========================================================================== #
# Figuras
# =========================================================================== #
def graficar_ojo(senal, sps_val, titulo, color, ax):
    seg = 2 * sps_val
    total = (len(senal) - 100 * sps_val) // seg
    t_ojo = np.arange(seg) * (1e12 * TS)  # ps
    for k in range(min(total, 400)):
        idx = 100 * sps_val + k * seg
        ax.plot(t_ojo, senal[idx:idx + seg], color=color, alpha=0.08, lw=0.8)
    ax.set_title(titulo, fontsize=10, fontweight="bold")
    ax.set_xlabel("Time (ps)", fontsize=8)
    ax.set_ylabel("Voltage (V)", fontsize=8)
    ax.grid(True, alpha=0.2)


def fig_senal_transmitida(esc, sym_ventana, nombre_canal):
    """Señal transmitida x(t) vs recibida y(t) (transitorio), escenarios A y C."""
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)

    n_sym_ventana = 16
    for ax, escn, col, tit in [
        (axes[0], "A", "#e53e3e", "Scenario A — channel output (no equalization)"),
        (axes[1], "C", "#319795", "Scenario C — FFE + CTLE (co-design)"),
    ]:
        y = esc["limpio"][escn]
        fase, lag = _find_align(y, sym_ventana, SPS)
        off = lag * SPS + fase
        i0 = 150 * SPS
        sl = slice(i0, i0 + n_sym_ventana * SPS)
        t_ns = (np.arange(n_sym_ventana * SPS)) * TS * 1e9
        x_step = esc["x_a"][slice(i0 - off, i0 - off + n_sym_ventana * SPS)] \
            if (i0 - off) >= 0 else esc["x_a"][sl]
        ax.step(t_ns, x_step, where="post", color="#718096", ls="--", alpha=0.8,
                label="Transmitted x(t) [PAM-4 symbols]")
        ax.plot(t_ns, y[sl], color=col, lw=2, label="Received y(t)")
        ax.set_title(tit, fontsize=10, fontweight="bold")
        ax.set_ylabel("Voltage (V)", fontsize=9)
        ax.grid(True, alpha=0.25)
        ax.legend(loc="upper right", fontsize=8)
    axes[1].set_xlabel("Time (ns)", fontsize=9)
    fig.tight_layout()
    ruta = os.path.join(DIR_FIG, "fig_ojo_senal_transmitida.png")
    fig.savefig(ruta, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return ruta


def fig_ojos_escenarios(esc, sym_ventana, nombre_canal):
    """Ojos A/B/C (limpios) + niveles muestreados C vs C+DFE (ruidosos)."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    graficar_ojo(esc["limpio"]["A"], SPS, "A — No equalization", "#e53e3e", axes[0, 0])
    graficar_ojo(esc["limpio"]["B"], SPS, "B — FFE only (TX)", "#dd6b20", axes[0, 1])
    graficar_ojo(esc["limpio"]["C"], SPS, "C — FFE + CTLE", "#319795", axes[1, 0])

    # Panel D: niveles muestreados en el instante de decisión, C vs C+DFE
    ax = axes[1, 1]
    rC = evaluar_enlace(esc["ruidoso"]["C"], sym_ventana, SPS, 0.0, devolver_muestras=True)
    rD = evaluar_enlace(esc["ruidoso"]["C"], sym_ventana, SPS, 0.0,
                        dfe_taps=esc["dfe_taps"], devolver_muestras=True)
    ax.scatter(np.arange(len(rC["muestras"])), rC["muestras"], s=4, alpha=0.25,
               color="#319795", label="C: FFE+CTLE")
    ax.scatter(np.arange(len(rD["muestras"])), rD["muestras"], s=4, alpha=0.25,
               color="#6b46c1", label="D: +DFE")
    for lv in NIVELES_PAM4:
        ax.axhline(lv * (rC["mu"][-1] / NIVELES_PAM4[-1]), color="0.6", lw=0.6, ls=":")
    ax.set_title("D — Sampled levels: DFE effect", fontsize=10, fontweight="bold")
    ax.set_xlabel("Symbol index", fontsize=8)
    ax.set_ylabel("Sampled value (V)", fontsize=8)
    ax.grid(True, alpha=0.2)
    ax.legend(fontsize=8, loc="upper right")

    fig.tight_layout()
    ruta = os.path.join(DIR_FIG, "fig_ojo_escenarios.png")
    fig.savefig(ruta, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return ruta


def fig_comparativa_canal(escenarios_por_canal):
    """Ojo del escenario C: RC vs FR-4 vs Megtron 6 (canal físico vs difusivo)."""
    nombres = list(escenarios_por_canal.keys())
    fig, axes = plt.subplots(1, len(nombres), figsize=(5 * len(nombres), 4.3), sharey=True)
    colores = {"RC (cosh)": "k", "FR-4": "#2b6cb0", "Megtron 6": "#dd6b20"}
    for ax, nb in zip(np.atleast_1d(axes), nombres):
        graficar_ojo(escenarios_por_canal[nb]["limpio"]["C"], SPS, nb,
                     colores.get(nb, "#319795"), ax)
    fig.tight_layout()
    ruta = os.path.join(DIR_FIG, "fig_ojo_comparativa_canal.png")
    fig.savefig(ruta, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return ruta


# =========================================================================== #
def main():
    os.makedirs(DIR_FIG, exist_ok=True)
    print("== Recomendación 01 — Ojo y señal transmitida sobre canal físico ==")
    print(f"PAM-4 {RS_BAUD/1e9:.0f} GBaud | sps={SPS} | sigma_in={SIGMA_IN} V | "
          f"co-diseño: w_z=2π·2.35GHz, c_post={C_POST}")

    # Secuencia PRBS23 -> PAM-4
    bits = generar_prbs(2 * CANT_SIMBOLOS + 600, orden=23)
    sym = modular_pam4(bits)[:CANT_SIMBOLOS]

    # Canal RC de referencia y su pérdida en Nyquist (para igualar longitudes)
    rc = CanalRC()
    il_obj = float(-20 * np.log10(np.abs(rc.H(np.array([OPERACION.f_nyquist])))[0]))

    # Definir los canales como proveedores de H(f)
    canales = {"RC (cosh)": (lambda f: rc.H(f), None)}
    for sub in SUBSTRATOS:
        ln = LineaTransmision(sub)
        long_m = ln.longitud_para_il(il_obj, OPERACION.f_nyquist)
        canales[sub.nombre] = ((lambda f, ln=ln, L=long_m: ln.H(f, L)), long_m)

    # Construir escenarios por canal y evaluar BER
    escenarios = {}
    filas = []
    for nb, (Hf, long_m) in canales.items():
        rng = np.random.default_rng(2026)
        esc = construir_escenarios(Hf, sym, rng)
        escenarios[nb] = esc

        rA = evaluar_enlace(esc["ruidoso"]["A"], sym, SPS, 0.0)
        rB = evaluar_enlace(esc["ruidoso"]["B"], sym, SPS, 0.0)
        rC = evaluar_enlace(esc["ruidoso"]["C"], sym, SPS, 0.0)
        rD = evaluar_enlace(esc["ruidoso"]["C"], sym, SPS, 0.0, dfe_taps=esc["dfe_taps"])
        long_txt = "-" if long_m is None else f"{long_m/INCH:.1f}"
        for escn, r in [("A sinEQ", rA), ("B FFE", rB), ("C FFE+CTLE", rC), ("D +DFE", rD)]:
            filas.append({
                "canal": nb, "long_pulg": long_txt, "escenario": escn,
                "BER": f"{r['ber']:.2e}",
                "apertura_ojo_V": round(r["eye_min"], 3),
                "estado": "ABIERTO" if r["eye_min"] > 0 else "cerrado",
            })

    # Figuras: usar FR-4 como canal físico de demostración
    nb_demo = FR4.nombre
    r1 = fig_senal_transmitida(escenarios[nb_demo], sym, nb_demo)
    r2 = fig_ojos_escenarios(escenarios[nb_demo], sym, nb_demo)
    r3 = fig_comparativa_canal(escenarios)

    # Tabla
    cols = ["canal", "long_pulg", "escenario", "BER", "apertura_ojo_V", "estado"]
    anchos = {c: max(len(c), *(len(str(f[c])) for f in filas)) for c in cols}
    linea = " | ".join(c.ljust(anchos[c]) for c in cols)
    print("\n" + linea); print("-" * len(linea))
    for f in filas:
        print(" | ".join(str(f[c]).ljust(anchos[c]) for c in cols))

    ruta_csv = os.path.join(DIR_FIG, "tabla_ojo_ber.csv")
    with open(ruta_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader(); w.writerows(filas)

    print("\nEntregables generados:")
    for r in (r1, r2, r3, ruta_csv):
        print(f"  - {r}")


if __name__ == "__main__":
    main()
