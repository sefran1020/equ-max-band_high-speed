"""
r1c3_constelacion.py — Revisor 1, comentario 3.

    "How hardware and software implementation for generating the graphical
     results including constellation diagram and BER are correlated?"

Aclaración previa: en este trabajo **no hay implementación en hardware**. Todos
los resultados provienen de simulación, y así se declara. Lo que sí existe es
una segunda implementación *independiente* de la misma cadena en un simulador de
circuitos (LTspice), que sirve para comprobar que las figuras no dependen del
programa que las produce. Este script responde midiendo esa correlación
exactamente sobre las dos magnitudes que el revisor menciona —el diagrama de
constelación y el BER— y no solo sobre la apertura del ojo, que era lo único
que el manuscrito comparaba.

Procedimiento
-------------
Se parte de los dos conjuntos ya generados por el flujo de validación:

  - `ref_python.npz`      : salida de la cadena evaluada en frecuencia (Python),
                            400 símbolos PAM-4 a 4 GBaud sobre el canal RC.
  - `cadena_completa.raw` : la MISMA cadena resuelta por LTspice, excitada con la
                            misma secuencia y la misma realización de ruido.

Ambas formas de onda se remuestrean a una rejilla común, se procesan con el
MISMO detector del artículo —se importa `detectar()` de
`r1c4_detector_pdf_niveles.py`, no se reimplementa— y se comparan:

  - constelación: valor muestreado por símbolo, con sus umbrales de decisión;
  - estadística por nivel: mu_k y sigma_k en cada motor;
  - BER semi-analítico de cada motor;
  - correlación muestra a muestra entre ambos (coeficiente de Pearson y RMS).

Entregables en ./figuras/ y ./tablas/:
  - fig_constelacion_motores.png : (a) constelación del motor en frecuencia,
                                   (b) constelación de LTspice, (c) correlación
                                   muestra a muestra entre ambos.
  - tabla_r1c3_constelacion.csv  : mu_k, sigma_k y BER por motor, y su desviación.

Uso:
    python revision01/r1c3_constelacion.py
"""

import csv
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(AQUI)
sys.path.insert(0, RAIZ)
sys.path.insert(0, AQUI)
sys.path.insert(0, os.path.join(RAIZ, "validacion_ltspice"))

from cadena_enlace import NIVELES_PAM4                      # noqa: E402
from lt_raw import leer_raw                                 # noqa: E402
from r1c4_detector_pdf_niveles import detectar              # noqa: E402

DIR_LT = os.path.join(RAIZ, "validacion_ltspice")
DIR_FIG = os.path.join(AQUI, "figuras")
DIR_TAB = os.path.join(AQUI, "tablas")
DPI = 300

SPU = 32                 # muestras por símbolo de la rejilla común
N_DESCARTAR = 20         # transitorio (la secuencia tiene solo 400 símbolos)
N_SCORE = 300            # ventana de correlación: debe dejar sitio al retardo
COLORES = ("#e53e3e", "#dd6b20", "#319795", "#2b6cb0")


PISO_BER = 1e-12         # piso de reporte del artículo


def fmt_ber(b, tex=False):
    """El artículo no reporta BER por debajo de 1e-12: ahí la cola gaussiana
    deja de ser fiable, así que se declara el piso en lugar de un número."""
    if b < PISO_BER:
        return r"BER $< 10^{-12}$" if tex else "<1e-12 (bajo el piso de reporte)"
    return (r"BER $= %.1e$" % b) if tex else "%.3e" % b


def remuestrear(t, v, dt):
    """Lleva una traza de paso variable a una rejilla uniforme."""
    tu = np.arange(t[0], t[-1], dt)
    return np.interp(tu, t, v)


def cargar():
    ref = np.load(os.path.join(DIR_LT, "ref_python.npz"))
    raw = leer_raw(os.path.join(DIR_LT, "cadena_completa.raw"))
    dt = float(ref["Tui"]) / SPU
    y_py = remuestrear(ref["t"], ref["vout"], dt)
    y_lt = remuestrear(raw["x"], raw["V(out)"], dt)
    n = min(len(y_py), len(y_lt))
    return y_py[:n], y_lt[:n], np.asarray(ref["syms"], dtype=float), \
        float(ref["baud"]), str(ref["canal"])


def graficar(det_py, det_lt, r_pearson, rms, ruta):
    """Figura a UNA columna: la constelación superpuesta demuestra el acuerdo
    mejor que dos paneles lado a lado, y ahorra media página en el artículo."""
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9))

    ax = axes[0]
    k = np.arange(len(det_py["m"]))
    for j, lv in enumerate(NIVELES_PAM4):
        sel = np.isclose(det_py["tx"], lv)
        ax.plot(k[sel], det_py["m"][sel], ".", color=COLORES[j], ms=3.0,
                alpha=0.85)
        sel_lt = np.isclose(det_lt["tx"], lv)
        ax.plot(np.arange(len(det_lt["m"]))[sel_lt], det_lt["m"][sel_lt], "x",
                color="0.25", ms=2.6, mew=0.5, alpha=0.55)
    for t in det_py["thr"]:
        ax.axhline(t, color="0.35", ls=":", lw=1.1)
    ax.plot([], [], ".", color="0.5", ms=6, label="frequency-domain chain")
    ax.plot([], [], "x", color="0.25", ms=5, mew=0.8, label="LTspice")
    ax.set_xlabel("(a) Symbol index")
    ax.set_ylabel("Sampled voltage (V)")
    ax.legend(fontsize=7, loc="center left", framealpha=0.9)
    ax.grid(alpha=0.3)

    ax = axes[1]
    n = min(len(det_py["m"]), len(det_lt["m"]))
    a, b = det_py["m"][:n], det_lt["m"][:n]
    lim = [min(a.min(), b.min()) - 0.3, max(a.max(), b.max()) + 0.3]
    ax.plot(lim, lim, color="0.6", lw=1.1, ls="--")
    ax.plot(a, b, ".", color="#805ad5", ms=3.0, alpha=0.7)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("(b) Sample, frequency-domain chain (V)")
    ax.set_ylabel("Sample, LTspice (V)")
    ax.grid(alpha=0.3)
    ax.text(0.04, 0.95, "$r = %.5f$\nRMS $= %.1f$ mV" % (r_pearson, rms * 1e3),
            transform=ax.transAxes, va="top", fontsize=8)

    fig.tight_layout()
    fig.savefig(ruta, dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def main():
    os.makedirs(DIR_FIG, exist_ok=True)
    os.makedirs(DIR_TAB, exist_ok=True)

    y_py, y_lt, syms, baud, canal = cargar()
    print("== R1.3 — correlación entre implementaciones (canal %s, %.0f GBaud) =="
          % (canal.upper(), baud / 1e9))
    print("Secuencia común: %d símbolos PAM-4, misma realización de ruido."
          % len(syms))

    det_py = detectar(y_py, syms, sps=SPU, n_descartar=N_DESCARTAR,
                      n_score=N_SCORE)
    det_lt = detectar(y_lt, syms, sps=SPU, n_descartar=N_DESCARTAR,
                      n_score=N_SCORE)

    n = min(len(det_py["m"]), len(det_lt["m"]))
    a, b = det_py["m"][:n], det_lt["m"][:n]
    r_pearson = float(np.corrcoef(a, b)[0, 1])
    rms = float(np.sqrt(np.mean((a - b) ** 2)))
    escala = float(np.max(np.abs(a)))

    print("\n-- Estadística por nivel (V)")
    print("   nivel      mu (freq)   mu (LTspice)   dif      "
          "sigma (freq)  sigma (LTspice)")
    for k, lv in enumerate(NIVELES_PAM4):
        print("   %+d       %+8.4f     %+8.4f    %+7.4f    %8.4f       %8.4f"
              % (lv, det_py["mu"][k], det_lt["mu"][k],
                 det_lt["mu"][k] - det_py["mu"][k],
                 det_py["sigma"][k], det_lt["sigma"][k]))

    print("\n-- Umbrales de decisión (V)")
    print("   frecuencia : %s" % np.array2string(det_py["thr"], precision=4,
                                                 floatmode="fixed"))
    print("   LTspice    : %s" % np.array2string(det_lt["thr"], precision=4,
                                                 floatmode="fixed"))

    print("\n-- Correlación entre motores")
    print("   Pearson r        : %.6f" % r_pearson)
    print("   RMS de diferencia: %.4f V  (%.2f %% del fondo de escala)"
          % (rms, 100 * rms / escala))
    print("   BER              : %s (frecuencia)  vs  %s (LTspice)"
          % (fmt_ber(det_py["ber"]), fmt_ber(det_lt["ber"])))
    print("   Apertura de ojo  : %.4f V  vs  %.4f V  (dif %.4f V)"
          % (det_py["eye_min"], det_lt["eye_min"],
             det_lt["eye_min"] - det_py["eye_min"]))

    ruta_fig = os.path.join(DIR_FIG, "fig_constelacion_motores.png")
    ruta_tab = os.path.join(DIR_TAB, "tabla_r1c3_constelacion.csv")
    graficar(det_py, det_lt, r_pearson, rms, ruta_fig)

    with open(ruta_tab, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["magnitud", "nivel", "cadena_frecuencia", "ltspice", "diferencia"])
        for k, lv in enumerate(NIVELES_PAM4):
            w.writerow(["mu_V", "%+d" % lv, "%.4f" % det_py["mu"][k],
                        "%.4f" % det_lt["mu"][k],
                        "%+.4f" % (det_lt["mu"][k] - det_py["mu"][k])])
        for k, lv in enumerate(NIVELES_PAM4):
            w.writerow(["sigma_V", "%+d" % lv, "%.4f" % det_py["sigma"][k],
                        "%.4f" % det_lt["sigma"][k],
                        "%+.4f" % (det_lt["sigma"][k] - det_py["sigma"][k])])
        w.writerow(["BER", "-", fmt_ber(det_py["ber"]), fmt_ber(det_lt["ber"]), ""])
        w.writerow(["apertura_ojo_V", "-", "%.4f" % det_py["eye_min"],
                    "%.4f" % det_lt["eye_min"],
                    "%+.4f" % (det_lt["eye_min"] - det_py["eye_min"])])
        w.writerow(["pearson_r", "-", "%.6f" % r_pearson, "", ""])
        w.writerow(["rms_diferencia_V", "-", "%.4f" % rms, "", ""])

    print("\nEntregables:\n  - %s\n  - %s" % (ruta_fig, ruta_tab))


if __name__ == "__main__":
    main()
