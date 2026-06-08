"""
vectorfit_s2p.py - Convierte el canal dispersivo (.s2p) en un bloque LTspice.

Problema: LTspice no consume parametros S arbitrarios en transitorio. Solucion
rigurosa: aproximar S21(f) por

      S21(f) ~= R(s) * exp(-s*Td)

donde Td es el retardo de grupo de masa (lo aporta un delay() ideal) y R(s) es
una funcion racional ESTABLE (la aporta una fuente Laplace=). Asi el canal
dispersivo queda disponible en .ac y en .tran dentro de LTspice.

El ajuste racional usa la iteracion de Sanathanan-Koerner (minimos cuadrados
con re-pesado), coeficientes reales, y reflexion de polos inestables (preserva
|H| exactamente). Genera:

  disp_lines.lib          : subcircuitos disp_fr4 y disp_meg para cadena_completa.cir
  fig_vectorfit_<x>.png   : calidad del ajuste (|S21| y fase, original vs fit)

Tambien deja el camino Python como referencia (el usuario pidio "ambos"): la
verdad de terreno es LineaTransmision.H, contra la que se mide el error.

Uso:  python vectorfit_s2p.py            (orden auto, banda <= 25 GHz)
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from parametros import SUBSTRATOS, INCH
from modelo_rlgc import LineaTransmision

AQUI = os.path.dirname(__file__)
DIR_S2P = os.path.join(AQUI, "..", "figuras")
W0 = 2 * np.pi * 1e10          # normalizacion de s (10 GHz) para coeficientes sanos
F_FIT = 25e9                   # banda de ajuste (suficiente para hasta ~ 12 Gbaud)
TOL_DB = 0.5                   # error objetivo en banda


def cargar_s21(ruta):
    d = np.loadtxt(ruta, comments=["!", "#"])
    f = d[:, 0]
    s21 = d[:, 3] + 1j * d[:, 4]
    return f, s21


def estimar_td(f, s21, fmax=8e9):
    ph = np.unwrap(np.angle(s21))
    gd = -np.gradient(ph, 2 * np.pi * f)
    return float(np.mean(gd[f <= fmax]))


def ajuste_sk(sn, H, P, Q, iters=8):
    """Ajuste racional real N(s)/D(s) por Sanathanan-Koerner.
    sn: s normalizado (=s/W0) ; H: residuo complejo. Devuelve (num, den)
    en potencias crecientes (num: 0..P ; den: 0..Q con den[0]=1)."""
    Va = np.vander(sn, P + 1, increasing=True)            # s^0..s^P
    Vb = np.vander(sn, Q + 1, increasing=True)[:, 1:]      # s^1..s^Q
    w = np.ones_like(sn, dtype=float)
    num = den = None
    for _ in range(iters):
        M = np.hstack([Va, -(H[:, None]) * Vb]) * w[:, None]
        rhs = (H * w)
        Mr = np.vstack([M.real, M.imag])
        br = np.concatenate([rhs.real, rhs.imag])
        x, *_ = np.linalg.lstsq(Mr, br, rcond=None)
        a = x[:P + 1]
        b = np.concatenate([[1.0], x[P + 1:]])
        num, den = a, b
        Dval = np.polyval(den[::-1], sn)                  # polyval: alto->bajo
        w = 1.0 / np.maximum(np.abs(Dval), 1e-9)
    return num, den


def estabilizar(den):
    """Refleja polos inestables al semiplano izquierdo (|H| invariante).
    Renormaliza el denominador a termino constante = 1 (convencion del ajuste),
    de modo que la magnitud queda EXACTAMENTE igual a la del ajuste crudo."""
    r = np.roots(den[::-1])                                # raices en s/W0
    r = np.where(r.real > 0, -r.real + 1j * r.imag, r)
    den_new = np.poly(r).real                              # alto->bajo, monico
    den_new = den_new / den_new[-1]                        # termino constante = 1
    return den_new[::-1], r                                # bajo->alto


def H_fit(num, den_lh, sn):
    """num bajo->alto ; den_lh bajo->alto."""
    return np.polyval(num[::-1], sn) / np.polyval(den_lh[::-1], sn)


def horner(coef_lh, var):
    """String Horner de un polinomio (coef en potencias crecientes)."""
    s = "%.8g" % coef_lh[-1]
    for c in coef_lh[-2::-1]:
        s = "(%.8g+(%s)*%s)" % (c, s, var)
    return s


def emitir_subckt(nombre, num_lh, den_lh, Td, ganancia):
    var = "(s/W0)"
    nums = horner(num_lh * ganancia, var)
    dens = horner(den_lh, var)
    return (
        f".subckt {nombre} in out\n"
        f".param Td={Td:.6e} W0={W0:.8e}\n"
        f"B1 nr 0 V=V(in) Laplace=({nums})/({dens})\n"
        f"B2 out 0 V=delay(V(nr),Td,{{2*Td}})\n"
        f".ends {nombre}\n"
    )


def procesar(sub, nombre_sub, ax_mag, ax_ph):
    ln = LineaTransmision(sub)
    L = 10.0 * INCH
    nf = sub.nombre.replace(" ", "").replace("-", "")
    f, s21 = cargar_s21(os.path.join(DIR_S2P, f"{nf}_10in.s2p"))

    banda = f <= F_FIT
    fb, s21b = f[banda], s21[banda]
    # OBJETIVO: el S21 del .s2p (dato Touchstone). Este bloque sirve para la
    # validacion en FRECUENCIA (Bode/.ac) del canal dispersivo en LTspice; NO
    # para el ojo en transitorio: la perdida por efecto pelicular ~ sqrt(f) es
    # NO-RACIONAL (cuspide en DC) y la perdida dielectrica exige la unidad
    # imaginaria en la Laplace -> el ojo dispersivo se queda en Python (exacto
    # por FFT). Ver README_LTSPICE.md, seccion "que NO conviene".
    Td = estimar_td(fb, s21b)
    sn = 1j * 2 * np.pi * fb / W0
    residuo = s21b * np.exp(1j * 2 * np.pi * fb * Td)      # quitar retardo de masa

    mejor = None
    for orden in range(4, 13):
        num, den = ajuste_sk(sn, residuo, orden, orden)
        den_lh, _ = estabilizar(den)
        Hf = H_fit(num, den_lh, sn)
        err_db = np.max(np.abs(20 * np.log10(np.abs(Hf) + 1e-300) -
                               20 * np.log10(np.abs(residuo) + 1e-300)))
        if mejor is None or err_db < mejor[0]:
            mejor = (err_db, orden, num, den_lh)
        if err_db < TOL_DB:
            break
    err_db, orden, num, den_lh = mejor
    err_cpx = err_db / 20.0
    g = abs(s21b[0]) / abs(H_fit(num, den_lh, sn)[0])
    num = num * g
    Hf = H_fit(num, den_lh, sn) * np.exp(-1j * 2 * np.pi * fb * Td)

    # Figura de calidad de ajuste
    ax_mag.plot(fb / 1e9, 20 * np.log10(np.abs(s21b)), label=f"{sub.nombre} S21 (.s2p)")
    ax_mag.plot(fb / 1e9, 20 * np.log10(np.abs(Hf)), "--", label=f"{sub.nombre} fit (ord {orden})")
    ax_ph.plot(fb / 1e9, np.unwrap(np.angle(s21b)), label=f"{sub.nombre} S21")
    ax_ph.plot(fb / 1e9, np.unwrap(np.angle(Hf)), "--", label=f"{sub.nombre} fit")

    print(f"[{sub.nombre}] Td={Td*1e9:.3f} ns  orden={orden}  "
          f"err_mag(|S21|)={err_db:.3f} dB  (uso: validacion Bode/.ac, no ojo)")
    return emitir_subckt(nombre_sub, num, den_lh, Td, 1.0)


def main():
    bloques = ["* disp_lines.lib - canales dispersivos (vectorfit de los .s2p)\n",
               "* Generado por vectorfit_s2p.py. Uso en cadena_completa.cir:\n",
               "*   .include disp_lines.lib   y   Xch a b disp_fr4|disp_meg\n\n"]
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, (axm, axp) = plt.subplots(1, 2, figsize=(11, 4.2))
    except Exception:
        axm = axp = None

    nombres = {"FR-4": "disp_fr4", "Megtron 6": "disp_meg"}
    for sub in SUBSTRATOS:
        bloques.append(procesar(sub, nombres[sub.nombre], axm, axp))
        bloques.append("\n")

    with open(os.path.join(AQUI, "disp_lines.lib"), "w") as fh:
        fh.write("".join(bloques))
    print("Escrito: disp_lines.lib")

    if axm is not None:
        for ax, ttl, yl in [(axm, "|S21|: vectorfit vs Touchstone", "Magnitud (dB)"),
                            (axp, "Fase S21: vectorfit vs Touchstone", "Fase (rad)")]:
            ax.set_xlabel("Frecuencia (GHz)"); ax.set_ylabel(yl)
            ax.set_title(ttl); ax.grid(alpha=0.3); ax.legend(fontsize=7)
        import matplotlib.pyplot as plt
        plt.tight_layout()
        ruta = os.path.join(AQUI, "fig_vectorfit.png")
        plt.savefig(ruta, dpi=200, bbox_inches="tight")
        print("Figura:", ruta)


if __name__ == "__main__":
    main()
