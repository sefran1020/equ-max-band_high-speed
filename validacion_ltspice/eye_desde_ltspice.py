"""
eye_desde_ltspice.py - Diagrama de ojo desde la salida de cadena_completa.cir
y CRUCE cuantitativo LTspice <-> Python (gemelo de ref_python.npz).

Hace dos cosas:
  1) Solapamiento de formas de onda V(out): LTspice vs el gemelo de Python,
     alineadas por correlacion, con el error RMS en banda (prueba de que ambos
     motores implementan la misma cadena LTI).
  2) Diagrama de ojo (2 UI) de la salida de LTspice, con la apertura estimada.

Salida: fig_eye_ltspice.png  (+ metricas por consola)

Uso:  python eye_desde_ltspice.py [cadena_completa.raw]
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from lt_raw import leer_raw

AQUI = os.path.dirname(__file__)
SPU = 100          # muestras por UI en la rejilla del ojo


def cargar():
    ref = np.load(os.path.join(AQUI, "ref_python.npz"))
    raw = sys.argv[1] if len(sys.argv) > 1 else os.path.join(AQUI, "cadena_completa.raw")
    d = leer_raw(raw)
    return ref, d


def remuestrear(t, v, dt):
    tu = np.arange(t[0], t[-1], dt)
    return tu, np.interp(tu, t, v)


def main():
    ref, d = cargar()
    Tui = float(ref["Tui"]); baud = float(ref["baud"])
    t_py, v_py = ref["t"], ref["vout"]
    t_lt, v_lt = d["x"], d["V(out)"]

    dt = Tui / SPU
    tu_py, vu_py = remuestrear(t_py, v_py, dt)
    tu_lt, vu_lt = remuestrear(t_lt, v_lt, dt)

    # --- alinear por correlacion (mismo TX, retardos de grupo casi iguales) -- #
    n = min(len(vu_py), len(vu_lt))
    a = vu_py[:n] - vu_py[:n].mean()
    b = vu_lt[:n] - vu_lt[:n].mean()
    cc = np.correlate(b, a, mode="full")
    lag = np.argmax(cc) - (n - 1)
    if lag >= 0:
        p, l = vu_py[: n - lag], vu_lt[lag:n]
    else:
        p, l = vu_py[-lag:n], vu_lt[: n + lag]
    m = min(len(p), len(l))
    # banda interior: descartar 10 UI de calentamiento y cola
    g = 10 * SPU
    pi_, li_ = p[g:m - g], l[g:m - g]
    rms = float(np.sqrt(np.mean((pi_ - li_) ** 2)))
    span = float(np.ptp(pi_))
    print("Cruce LTspice<->Python:  lag=%d muestras  RMS=%.4f  (%.2f%% del span %.3f)"
          % (lag, rms, 100 * rms / span, span))

    # --- diagrama de ojo de LTspice (2 UI) --------------------------------- #
    g0 = 10 * SPU
    v = vu_lt[g0: len(vu_lt) - SPU]
    nfil = (len(v) - 2 * SPU) // SPU
    eje = np.linspace(0, 2.0, 2 * SPU)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.3))

        # (1) solapamiento de formas de onda (primeras ~30 UI utiles)
        nshow = 30 * SPU
        tt = np.arange(min(nshow, len(pi_))) / SPU
        ax1.plot(tt, pi_[:len(tt)], label="Python (gemelo)", lw=1.4)
        ax1.plot(tt, li_[:len(tt)], "--", label="LTspice", lw=1.0)
        ax1.set_xlabel("Tiempo (UI)"); ax1.set_ylabel("V(out)")
        ax1.set_title("Cruce de formas de onda  (RMS=%.3f)" % rms)
        ax1.legend(fontsize=8); ax1.grid(alpha=0.3)

        # (2) ojo
        for k in range(nfil):
            seg = v[k * SPU: k * SPU + 2 * SPU]
            if len(seg) == 2 * SPU:
                ax2.plot(eje, seg, color="C0", alpha=0.05, lw=0.6)
        ax2.set_xlabel("Tiempo (UI)"); ax2.set_ylabel("V(out)")
        ax2.set_title("Diagrama de ojo PAM-4 (LTspice, %.0f GBaud)" % (baud / 1e9))
        ax2.grid(alpha=0.3)
        plt.tight_layout()
        ruta = os.path.join(AQUI, "fig_eye_ltspice.png")
        plt.savefig(ruta, dpi=200, bbox_inches="tight")
        print("Figura:", ruta)
    except Exception as e:
        print("(figura omitida: %s)" % e)


if __name__ == "__main__":
    main()
