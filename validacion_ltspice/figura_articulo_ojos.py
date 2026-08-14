"""
figura_articulo_ojos.py - Figura sobria (lista para el manuscrito) del cruce de
ojos LTspice vs pipeline de Python sobre el canal RC.

Usa ref_python.npz (gemelo del pipeline) y cadena_completa.raw (LTspice), ya
generados por el flujo de validacion. Guarda la figura directamente en la carpeta
de figuras del articulo, con estilo neutro (sin titulo promocional).

Uso:  python figura_articulo_ojos.py
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from lt_raw import leer_raw

AQUI = os.path.dirname(__file__)
DEST = os.path.normpath(os.path.join(AQUI, "fig_ltspice_ojos.png"))
SPU = 100


def remuestrear(t, v, dt):
    tu = np.arange(t[0], t[-1], dt)
    return np.interp(tu, t, v)


def ojo(ax, vu, color, panel=None):
    g0 = 10 * SPU
    v = vu[g0: len(vu) - SPU]
    nfil = (len(v) - 2 * SPU) // SPU
    eje = np.linspace(0, 2.0, 2 * SPU)
    for k in range(nfil):
        seg = v[k * SPU: k * SPU + 2 * SPU]
        if len(seg) == 2 * SPU:
            ax.plot(eje, seg, color=color, alpha=0.05, lw=0.6)
    if panel:
        ax.set_xlabel("%s\nTime (UI)" % panel)
    else:
        ax.set_xlabel("Time (UI)")
    ax.grid(alpha=0.3)


def main():
    ref = np.load(os.path.join(AQUI, "ref_python.npz"))
    Tui = float(ref["Tui"])
    d = leer_raw(os.path.join(AQUI, "cadena_completa.raw"))
    dt = Tui / SPU
    vu_py = remuestrear(ref["t"], ref["vout"], dt)
    vu_lt = remuestrear(d["x"], d["V(out)"], dt)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8), sharey=True)
    ojo(ax[0], vu_py, "#2b6cb0", "(a) Frequency-domain chain (Python)")
    ojo(ax[1], vu_lt, "#c05621", "(b) Circuit simulator (LTspice)")
    ojo(ax[2], vu_py, "#2b6cb0", "(c) Overlay")
    ojo(ax[2], vu_lt, "#c05621", "(c) Overlay")
    ax[0].set_ylabel("Voltage (V)")
    plt.tight_layout()
    plt.savefig(DEST, dpi=300, bbox_inches="tight")
    print("Figura del articulo:", DEST)


if __name__ == "__main__":
    main()
