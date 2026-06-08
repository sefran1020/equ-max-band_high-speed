"""
comparar_bode.py - Cruza el Bode (.ac) de LTspice contra la H(f) de Python.

Lee rc_bode.raw y ctle_bode.raw (genéralos con el runner o LTspice -b) y los
compara con modelo_rc.CanalRCEscalera.H y cadena_enlace.ctle_response.
Reporta el error maximo en dB en banda y guarda fig_bode_cruce.png.

Uso:  python comparar_bode.py
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))
from lt_raw import leer_raw
from modelo_rc import CanalRCEscalera
from cadena_enlace import ctle_response, W_Z_RC, CTLE_ADC

AQUI = os.path.dirname(__file__)
DIR_S2P = os.path.join(AQUI, "..", "figuras")


def s21_de(nombre_archivo):
    d = np.loadtxt(os.path.join(DIR_S2P, nombre_archivo), comments=["!", "#"])
    return d[:, 0], d[:, 3] + 1j * d[:, 4]


def db(x):
    return 20 * np.log10(np.abs(x) + 1e-300)


def main():
    casos = []
    rc_raw = os.path.join(AQUI, "rc_bode.raw")
    ctle_raw = os.path.join(AQUI, "ctle_bode.raw")
    if os.path.exists(rc_raw):
        d = leer_raw(rc_raw); f = d["x"]
        H_py = CanalRCEscalera(80.0, 5e-12, 5).H(f)
        casos.append(("Canal RC N=5", f, d["V(out)"], H_py))
    if os.path.exists(ctle_raw):
        d = leer_raw(ctle_raw); f = d["x"]
        H_py = ctle_response(f, W_Z_RC, A_dc=CTLE_ADC)
        casos.append(("CTLE", f, d["V(out)"], H_py))

    disp_raw = os.path.join(AQUI, "disp_bode.raw")
    if os.path.exists(disp_raw):
        d = leer_raw(disp_raw); f = d["x"]
        for nodo, s2p, nb in [("V(o_fr4)", "FR4_10in.s2p", "Disp. FR-4"),
                              ("V(o_meg)", "Megtron6_10in.s2p", "Disp. Megtron 6")]:
            if nodo in d:
                fo, s21 = s21_de(s2p)
                H_py = np.interp(f, fo, s21.real) + 1j * np.interp(f, fo, s21.imag)
                casos.append((nb, f, d[nodo], H_py))

    if not casos:
        print("No hay .raw de Bode. Corre:  LTspice -b rc_bode.cir  y  ctle_bode.cir")
        return

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 4.5))
    except Exception:
        ax = None

    for nombre, f, H_lt, H_py in casos:
        banda = f <= 15e9
        err = np.max(np.abs(db(H_lt[banda]) - db(H_py[banda])))
        print("[%s]  err_max(|H|) en <=15 GHz = %.4f dB" % (nombre, err))
        if ax is not None:
            ax.semilogx(f, db(H_lt), label="%s LTspice" % nombre)
            ax.semilogx(f, db(H_py), "--", label="%s Python" % nombre)

    if ax is not None:
        ax.set_xlabel("Frecuencia (Hz)"); ax.set_ylabel("Magnitud (dB)")
        ax.set_title("Bode: LTspice vs Python"); ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=8)
        import matplotlib.pyplot as plt
        plt.tight_layout()
        ruta = os.path.join(AQUI, "fig_bode_cruce.png")
        plt.savefig(ruta, dpi=200, bbox_inches="tight")
        print("Figura:", ruta)


if __name__ == "__main__":
    main()
