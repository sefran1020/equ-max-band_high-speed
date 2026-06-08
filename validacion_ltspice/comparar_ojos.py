"""
comparar_ojos.py - Ojo de LTspice  vs  ojo del pipeline de Python (lado a lado).

Compara, bajo CONDICIONES IDENTICAS (mismo PRBS, mismo canal, mismo ruido
muestra a muestra), el diagrama de ojo que produce LTspice (solver de circuitos
causal) con el que produce el pipeline de Python del articulo.

La referencia de Python es ref_python.npz, generado por generar_tx_pwl.py con
LOS MISMOS bloques del pipeline (generar_prbs / modular_pam4 / aplicar_ffe /
H del canal / ctle_response / Butterworth de RX) y el ruido IDENTICO inyectado
en LTspice. Asi, cualquier diferencia entre ojos es atribuible solo al motor.

Salida:
  fig_comparacion_ojos.png   (Python | LTspice | superpuestos)
  metricas por consola: apertura de cada sub-ojo (mu +- 3 sigma) en ambos motores

Uso:  python comparar_ojos.py [cadena_completa.raw]
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from lt_raw import leer_raw

AQUI = os.path.dirname(__file__)
SPU = 100
NIVELES = np.array([-3.0, -1.0, 1.0, 3.0])


def remuestrear(t, v, dt):
    tu = np.arange(t[0], t[-1], dt)
    return np.interp(tu, t, v)


def medir_ojo(vu, syms, spu, descartar=20):
    """Mejor instante de decision: apertura de los 3 sub-ojos (mu +- 3 sigma)."""
    nsym = len(syms)
    mejor = None
    for ph in range(spu):
        ds = vu[ph::spu]
        n = min(len(ds), nsym)
        if n < 60:
            continue
        a = ds[:n] - ds[:n].mean()
        b = syms[:n] - syms[:n].mean()
        cc = np.correlate(a, b, mode="full")
        lag = int(np.argmax(cc)) - (n - 1)
        if lag < 0 or lag + 60 > len(ds):
            continue
        k = min(len(ds) - lag, nsym)
        m = ds[lag:lag + k][descartar:]
        tx = syms[:k][descartar:]
        mu, sd, ok = [], [], True
        for lv in NIVELES:
            sel = m[np.isclose(tx, lv)]
            if len(sel) < 2:
                ok = False
                break
            mu.append(sel.mean()); sd.append(sel.std())
        if not ok:
            continue
        mu, sd = np.array(mu), np.array(sd)
        op = (mu[1:] - 3 * sd[1:]) - (mu[:-1] + 3 * sd[:-1])
        sc = float(np.min(op))
        if mejor is None or sc > mejor["score"]:
            mejor = {"score": sc, "op": op, "mu": mu, "sd": sd, "ph": ph}
    return mejor


def dibujar_ojo(ax, vu, spu, color, titulo):
    g0 = 10 * spu
    v = vu[g0: len(vu) - spu]
    nfil = (len(v) - 2 * spu) // spu
    eje = np.linspace(0, 2.0, 2 * spu)
    for k in range(nfil):
        seg = v[k * spu: k * spu + 2 * spu]
        if len(seg) == 2 * spu:
            ax.plot(eje, seg, color=color, alpha=0.05, lw=0.6)
    ax.set_xlabel("Tiempo (UI)"); ax.set_ylabel("V(out)")
    ax.set_title(titulo); ax.grid(alpha=0.3)


def main():
    ref = np.load(os.path.join(AQUI, "ref_python.npz"))
    Tui = float(ref["Tui"]); baud = float(ref["baud"]); syms = ref["syms"]
    canal = str(ref["canal"]) if "canal" in ref else "rc"
    raw = sys.argv[1] if len(sys.argv) > 1 else os.path.join(AQUI, "cadena_completa.raw")
    d = leer_raw(raw)

    dt = Tui / SPU
    vu_py = remuestrear(ref["t"], ref["vout"], dt)
    vu_lt = remuestrear(d["x"], d["V(out)"], dt)

    mp = medir_ojo(vu_py, syms, SPU)
    ml = medir_ojo(vu_lt, syms, SPU)

    print("== Comparacion de ojos (canal=%s, %.0f GBaud) ==" % (canal, baud / 1e9))
    print("  sub-ojo:           inferior   central   superior")
    print("  Python  apertura:  %8.3f %9.3f %9.3f" % tuple(mp["op"]))
    print("  LTspice apertura:  %8.3f %9.3f %9.3f" % tuple(ml["op"]))
    print("  |dif| apertura  :  %8.3f %9.3f %9.3f" % tuple(np.abs(mp["op"] - ml["op"])))
    print("  niveles mu Python : %s" % np.round(mp["mu"], 3))
    print("  niveles mu LTspice: %s" % np.round(ml["mu"], 3))
    dmu = float(np.max(np.abs(mp["mu"] - ml["mu"])))
    dop = float(np.max(np.abs(mp["op"] - ml["op"])))
    print("  max|d mu|=%.4f V   max|d apertura|=%.4f V" % (dmu, dop))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.4), sharey=True)
        dibujar_ojo(axes[0], vu_py, SPU, "#2b6cb0",
                    "Pipeline Python (%.0f GBaud)" % (baud / 1e9))
        dibujar_ojo(axes[1], vu_lt, SPU, "#dd6b20", "LTspice")
        # superpuestos
        dibujar_ojo(axes[2], vu_py, SPU, "#2b6cb0", "Superpuestos (Py azul / LT naranja)")
        dibujar_ojo(axes[2], vu_lt, SPU, "#dd6b20", "Superpuestos (Py azul / LT naranja)")
        # marcar instante de decision e niveles
        for ax in axes:
            for lv in mp["mu"]:
                ax.axhline(lv, color="0.6", lw=0.5, ls=":")
        fig.suptitle("Diagrama de ojo: pipeline Python vs LTspice  "
                     "(mismo PRBS / canal / ruido)  -  max|d apertura|=%.3f V" % dop,
                     fontsize=12, fontweight="bold")
        fig.tight_layout()
        ruta = os.path.join(AQUI, "fig_comparacion_ojos.png")
        fig.savefig(ruta, dpi=200, bbox_inches="tight")
        print("Figura:", ruta)
    except Exception as e:
        print("(figura omitida: %s)" % e)


if __name__ == "__main__":
    main()
