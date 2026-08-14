"""
r1c4_detector_pdf_niveles.py — Revisor 1, comentario 4.

    "The symbol detection mechanism at the receiver to generate BER is not clear."

El mecanismo ya existía en `cadena_enlace.evaluar_enlace()`, pero el manuscrito
lo describía solo en prosa. Este script lo hace explícito y verificable: expone
las tres etapas de la decisión y produce la figura y la tabla que se insertan en
la respuesta a los revisores.

Mecanismo de detección (lo que hace el receptor, paso a paso)
-------------------------------------------------------------
1. INSTANTE DE MUESTREO. Se toma una muestra por símbolo. La fase
   phi* in {0..sps-1} y el retardo entero l se eligen maximizando la correlación
   normalizada entre el flujo submuestreado y los símbolos transmitidos
   (`alinear_full`). Es una sincronización asistida por datos: acota el
   rendimiento por arriba y evita muestrear en el cruce del ojo.
2. REALIMENTACIÓN DE DECISIONES (opcional). Con DFE de nt taps, los coeficientes
   se ajustan por mínimos cuadrados sobre las muestras alineadas y la muestra
   corregida es  y_k = m_k - sum_i h_i * a_hat_{k-i} , donde a_hat es la decisión
   dura previa (nivel PAM-4 más próximo, escalado por la ganancia del cursor h0).
3. UMBRALES Y BER. Se estiman la media mu_k y la desviación sigma_k de cada uno
   de los cuatro niveles a partir de las muestras cuyo símbolo transmitido es
   ese nivel. Los umbrales son los puntos medios entre medias contiguas,
   tau_k = (mu_k + mu_{k+1})/2 — es decir, umbrales adaptados a la señal
   recibida, no los niveles nominales. La probabilidad de error por nivel se
   evalúa con la función Q sobre las colas que cruzan el umbral (los niveles
   interiores aportan dos colas), se pondera por la frecuencia de cada nivel
   para obtener la SER y, por el mapeo Gray de PAM-4 —donde un error de símbolo
   entre niveles contiguos cambia un solo bit de los dos—, BER = SER/2.

Verificación de identidad
-------------------------
El script recalcula el BER con las ecuaciones de arriba y lo compara contra el
valor que devuelve `evaluar_enlace()`, que es la función que produce TODOS los
BER del artículo. La coincidencia exacta demuestra que la figura no ilustra un
detector "de ejemplo", sino el que genera cada número publicado.

Entregables en ./figuras/ y ./tablas/:
  - fig_detector_niveles.png  : (a) elección del instante de muestreo,
                             (b) niveles, umbrales y ajuste Gaussiano (FFE+CTLE),
                             (c) lo mismo con DFE de 2 taps.
  - tabla_r1c4_detector.csv: mu_k, sigma_k, cuentas, Pe_k, umbrales, SER y BER.

Uso:
    python revision01/r1c4_detector_pdf_niveles.py
"""

import csv
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from cadena_enlace import (                                  # noqa: E402
    NIVELES_PAM4, SIGMA_IN, SPS, W_Z_RC, C_POST_RC, Q,
    generar_prbs, modular_pam4, aplicar_ffe, aplicar_canal_freq,
    aplicar_ctle_freq, aplicar_rx_lpf, alinear_full, estadisticas_niveles,
    _decidir_pam4, evaluar_enlace,
)
from canales import construir_canales                        # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
DIR_FIG = os.path.join(AQUI, "figuras")
DIR_TAB = os.path.join(AQUI, "tablas")
DPI = 300

# Punto de operación del artículo, canal y realización 0 del conjunto publicado.
CANAL = "FR-4"
CRITERIO = "equal_il"
RS = 4e9                # 4 GBaud (Nyquist 2 GHz)
N_SYM = 1000
N_DESCARTAR = 100       # transitorio descartado, igual que en evaluar_enlace()
N_DFE = 2
SEMILLA_RUIDO = 1000    # realización r=0  ->  semilla 1000+r

COLORES = ("#e53e3e", "#dd6b20", "#319795", "#2b6cb0")


# --------------------------------------------------------------------------- #
# Cadena hasta la entrada del detector (escenario C del artículo: FFE+CTLE)
# --------------------------------------------------------------------------- #
def forma_de_onda_rx():
    """Devuelve (y_rx, sym, ts): la señal en la entrada del detector."""
    canales, _ = construir_canales(CRITERIO)
    H = canales[CANAL]["H"]

    bits = generar_prbs(2 * N_SYM + 200, orden=23)
    sym = modular_pam4(bits)[:N_SYM]

    ts = 1.0 / (RS * SPS)
    rng = np.random.default_rng(SEMILLA_RUIDO)

    ffe = aplicar_ffe(sym, 0.0, 1.0, C_POST_RC)
    y_chan = aplicar_canal_freq(np.repeat(ffe, SPS), ts, H)
    # Ruido REFERIDO A LA ENTRADA: se suma antes del CTLE, de modo que el
    # ecualizador amplifica señal y ruido por igual.
    nin = rng.normal(0, SIGMA_IN, len(y_chan))
    y_rx = aplicar_rx_lpf(aplicar_ctle_freq(y_chan + nin, ts, W_Z_RC), ts, RS)
    return y_rx, sym, ts


# --------------------------------------------------------------------------- #
# Etapa 1: instante de muestreo
# --------------------------------------------------------------------------- #
def perfil_de_fase(y, sym, n_score=400):
    """Correlación normalizada |rho| contra la fase de muestreo phi.

    Reproduce el criterio interno de `alinear_full`, para poder dibujarlo.
    """
    n = len(y)
    x_up = np.repeat(sym, SPS)
    x_up = x_up[:n] if len(x_up) >= n else np.pad(x_up, (0, n - len(x_up)))
    cc = np.fft.irfft(np.fft.rfft(y - y.mean())
                      * np.conj(np.fft.rfft(x_up - x_up.mean())), n=n)
    lag0 = int(np.argmax(cc[: n // 2])) // SPS

    ns = min(n_score, len(sym))
    b = sym[:ns] - sym[:ns].mean()
    nb = np.linalg.norm(b)
    rho = np.zeros(SPS)
    for ph in range(SPS):
        ds = y[ph::SPS]
        mejor = 0.0
        for lag in (lag0 - 1, lag0, lag0 + 1):
            if lag < 0 or lag + ns > len(ds):
                continue
            a = ds[lag:lag + ns]
            a = a - a.mean()
            na = np.linalg.norm(a)
            if na > 0 and nb > 0:
                mejor = max(mejor, abs(float(np.dot(a, b)) / (na * nb)))
        rho[ph] = mejor
    return rho, int(np.argmax(rho))


# --------------------------------------------------------------------------- #
# Etapas 2 y 3: realimentación de decisiones, umbrales y BER
# --------------------------------------------------------------------------- #
def detectar(y, sym, n_dfe=0, sps=SPS, n_descartar=N_DESCARTAR, n_score=400):
    """Aplica el detector y devuelve las magnitudes que lo definen.

    `sps`, `n_descartar` y `n_score` son parámetros para que otros scripts de la
    revisión —en particular `r1c3_constelacion.py`, que procesa formas de onda
    de LTspice con otro sobremuestreo y secuencias más cortas— reutilicen ESTE
    detector en lugar de reimplementarlo.

    Cuidado con `n_score`: `alinear_full` solo evalúa retardos que dejen espacio
    para su ventana de correlación, de modo que con una secuencia de 400
    símbolos y la ventana por defecto de 400 el único retardo admisible sería 0
    y la alineación fallaría en silencio. La ventana debe ser holgadamente menor
    que la secuencia.
    """
    m, tx = alinear_full(y, sym, sps, n_score=n_score)
    m, tx = m[n_descartar:], tx[n_descartar:]

    if n_dfe > 0:
        N = len(m)
        Xreg = np.column_stack([tx[n_dfe - i: N - i] for i in range(n_dfe + 1)])
        h, *_ = np.linalg.lstsq(Xreg, m[n_dfe:N], rcond=None)
        h0 = h[0] if abs(h[0]) > 1e-12 else 1.0
        post = h[1:]
        hist = [0.0] * n_dfe
        m_eq = np.empty_like(m)
        for kk in range(len(m)):
            fb = sum(post[i] * hist[i] for i in range(n_dfe))
            yk = m[kk] - fb
            m_eq[kk] = yk
            hist = [_decidir_pam4(yk, h0)] + hist[:-1]
        m = m_eq

    mu, sig, cnt = estadisticas_niveles(m, tx)
    o = np.argsort(mu)
    mu, sig, cnt = mu[o], sig[o], cnt[o]
    thr = (mu[:-1] + mu[1:]) / 2.0          # umbrales = puntos medios

    eps = 1e-300
    Pe = np.zeros(4)
    Pe[0] = Q((thr[0] - mu[0]) / (sig[0] + eps))
    Pe[1] = Q((mu[1] - thr[0]) / (sig[1] + eps)) + Q((thr[1] - mu[1]) / (sig[1] + eps))
    Pe[2] = Q((mu[2] - thr[1]) / (sig[2] + eps)) + Q((thr[2] - mu[2]) / (sig[2] + eps))
    Pe[3] = Q((mu[3] - thr[2]) / (sig[3] + eps))

    pesos = cnt / cnt.sum()
    ser = float(np.sum(pesos * Pe))
    ber = ser / 2.0                          # Gray: 1 bit erróneo por error entre vecinos
    apertura = (mu[1:] - 3 * sig[1:]) - (mu[:-1] + 3 * sig[:-1])
    return {"m": m, "tx": tx, "mu": mu, "sigma": sig, "cnt": cnt, "thr": thr,
            "Pe": Pe, "ser": ser, "ber": ber, "eye_min": float(np.min(apertura))}


# --------------------------------------------------------------------------- #
# Figura
# --------------------------------------------------------------------------- #
def graficar(rho, ph_opt, det_c, det_d, ruta):
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.6))

    ax = axes[0]
    ui = np.arange(SPS) / SPS
    ax.plot(ui, rho, color="#2b6cb0", lw=2)
    ax.axvline(ph_opt / SPS, color="#c05621", ls="--", lw=1.5)
    ax.plot([ph_opt / SPS], [rho[ph_opt]], "o", color="#c05621", ms=7)
    ax.annotate(r"$\varphi^{*}$", xy=(ph_opt / SPS, rho[ph_opt]),
                xytext=(6, -14), textcoords="offset points",
                color="#c05621", fontsize=11)
    ax.set_xlabel("(a) Sampling phase (UI)")
    ax.set_ylabel(r"Normalized correlation $|\rho|$")
    ax.grid(alpha=0.3)

    for ax, det, etiqueta in ((axes[1], det_c, "(b) FFE+CTLE"),
                              (axes[2], det_d, "(c) FFE+CTLE+DFE")):
        for k, lv in enumerate(NIVELES_PAM4):
            sel = det["m"][np.isclose(det["tx"], lv)]
            ax.hist(sel, bins=40, density=True, alpha=0.45, color=COLORES[k])
            rejilla = np.linspace(det["mu"][k] - 4 * det["sigma"][k],
                                  det["mu"][k] + 4 * det["sigma"][k], 200)
            pdf = np.exp(-0.5 * ((rejilla - det["mu"][k]) / det["sigma"][k]) ** 2) \
                / (det["sigma"][k] * np.sqrt(2 * np.pi))
            ax.plot(rejilla, pdf, color=COLORES[k], lw=1.8,
                    label=r"$a_k=%+d$" % lv)
        for t in det["thr"]:
            ax.axvline(t, color="0.35", ls=":", lw=1.4)
        ax.set_xlabel("%s\nSampled voltage (V)   [dotted: decision thresholds $\\tau_k$]"
                      % etiqueta)
        ax.set_ylabel("Density")
        ax.grid(alpha=0.3)
        ax.text(0.02, 0.95, "BER = %.1e" % det["ber"], transform=ax.transAxes,
                va="top", fontsize=9)

    axes[2].legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(ruta, dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def escribir_tabla(det_c, det_d, ruta):
    with open(ruta, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["escenario", "nivel", "mu_V", "sigma_V", "cuenta", "Pe"])
        for nombre, det in (("FFE+CTLE", det_c), ("FFE+CTLE+DFE", det_d)):
            for k, lv in enumerate(NIVELES_PAM4):
                w.writerow([nombre, "%+d" % lv, "%.4f" % det["mu"][k],
                            "%.4f" % det["sigma"][k], int(det["cnt"][k]),
                            "%.3e" % det["Pe"][k]])
        w.writerow([])
        w.writerow(["escenario", "umbral_1_V", "umbral_2_V", "umbral_3_V",
                    "SER", "BER", "apertura_ojo_V"])
        for nombre, det in (("FFE+CTLE", det_c), ("FFE+CTLE+DFE", det_d)):
            w.writerow([nombre] + ["%.4f" % t for t in det["thr"]]
                       + ["%.3e" % det["ser"], "%.3e" % det["ber"],
                          "%.4f" % det["eye_min"]])


def main():
    os.makedirs(DIR_FIG, exist_ok=True)
    os.makedirs(DIR_TAB, exist_ok=True)

    y, sym, _ = forma_de_onda_rx()
    rho, ph_opt = perfil_de_fase(y, sym)
    det_c = detectar(y, sym, n_dfe=0)
    det_d = detectar(y, sym, n_dfe=N_DFE)

    # --- Identidad con la función que genera todos los BER del artículo ---- #
    ref_c = evaluar_enlace(y, sym, SPS)
    ref_d = evaluar_enlace(y, sym, SPS, dfe_ntaps=N_DFE)
    ok_c = np.isclose(det_c["ber"], ref_c["ber"], rtol=1e-12, atol=0.0)
    ok_d = np.isclose(det_d["ber"], ref_d["ber"], rtol=1e-12, atol=0.0)

    print("== R1.4 — mecanismo de detección (canal %s, %s, %.0f GBaud) =="
          % (CANAL, CRITERIO, RS / 1e9))
    print("Instante de muestreo: fase óptima phi* = %d/%d = %.3f UI  (|rho| = %.4f)"
          % (ph_opt, SPS, ph_opt / SPS, rho[ph_opt]))
    print("Peor fase: |rho| = %.4f  ->  la elección del instante no es cosmética"
          % rho.min())
    for nombre, det in (("FFE+CTLE", det_c), ("FFE+CTLE+DFE", det_d)):
        print("\n-- %s" % nombre)
        for k, lv in enumerate(NIVELES_PAM4):
            print("   nivel %+d : mu = %+7.4f V   sigma = %.4f V   n = %4d   Pe = %.3e"
                  % (lv, det["mu"][k], det["sigma"][k], det["cnt"][k], det["Pe"][k]))
        print("   umbrales tau = %s V"
              % np.array2string(det["thr"], precision=4, floatmode="fixed"))
        print("   SER = %.3e   BER = SER/2 = %.3e   apertura de ojo = %.4f V"
              % (det["ser"], det["ber"], det["eye_min"]))

    print("\nIdentidad con evaluar_enlace() (la función que produce los BER del artículo):")
    print("   FFE+CTLE     : %.6e vs %.6e  -> %s"
          % (det_c["ber"], ref_c["ber"], "COINCIDE" if ok_c else "DIFIERE"))
    print("   FFE+CTLE+DFE : %.6e vs %.6e  -> %s"
          % (det_d["ber"], ref_d["ber"], "COINCIDE" if ok_d else "DIFIERE"))
    if not (ok_c and ok_d):
        raise SystemExit("ERROR: el detector documentado no reproduce el del artículo.")

    ruta_fig = os.path.join(DIR_FIG, "fig_detector_niveles.png")
    ruta_tab = os.path.join(DIR_TAB, "tabla_r1c4_detector.csv")
    graficar(rho, ph_opt, det_c, det_d, ruta_fig)
    escribir_tabla(det_c, det_d, ruta_tab)
    print("\nEntregables:\n  - %s\n  - %s" % (ruta_fig, ruta_tab))


if __name__ == "__main__":
    main()
