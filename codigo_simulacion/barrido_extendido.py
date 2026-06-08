"""
barrido_extendido.py — Barrido de tasa extendido + sensibilidad de umbral
(Bloques C3 y C4 del plan).

C4: extiende el barrido más allá de 32 GBaud para RESOLVER las saturaciones
    (Megtron 6 'fixed' sin EQ / +DFE) que en la Tabla I aparecían como ">=32".
C3: calcula la tasa viable a TRES umbrales pre-FEC (1e-2, 1e-3, 1e-4) para
    mostrar que las conclusiones no dependen críticamente del umbral.

Misma metodología que barrido_tasa_fisico.py (30 realizaciones, EQ heredada del
RC), solo cambia la rejilla de tasas y el post-proceso multi-umbral.

Uso:  python barrido_extendido.py
"""

import csv
import os
import numpy as np

from cadena_enlace import (
    generar_prbs, modular_pam4, evaluar_config, viable_interp,
    W_Z_RC, C_POST_RC,
)
from canales import construir_canales

DIR_FIG = os.path.join(os.path.dirname(__file__), "figuras")

N_REAL = 30
N_SYM = 1000
N_DFE = 2
TASAS = np.linspace(1e9, 64e9, 43)        # hasta 64 GBaud (resuelve saturaciones)
GB = TASAS / 1e9
UMBRALES = [1e-2, 1e-3, 1e-4]
CONFIGS = [("none", "Sin EQ", "A"), ("lin", "FFE+CTLE", "C"), ("dfe", "FFE+CTLE+DFE", "D")]


def correr(criterio):
    canales, _ = construir_canales(criterio)
    bits_pool = generar_prbs(N_REAL * 2 * N_SYM + 5000, orden=23)
    # curvas[nb][key] = lista de (n_real) curvas de BER
    curvas = {nb: {k: [] for k, _, _ in CONFIGS} for nb in canales}
    print(f"\n== Criterio '{criterio}' ({N_REAL} realizaciones, hasta {GB[-1]:.0f} GBaud) ==")
    for r in range(N_REAL):
        bsl = bits_pool[r * 2 * N_SYM: r * 2 * N_SYM + 2 * N_SYM + 200]
        sym = modular_pam4(bsl)[:N_SYM]
        for nb, info in canales.items():
            rng = np.random.default_rng(1000 + r)
            bn, bl, bd = [], [], []
            for Rs in TASAS:
                res = evaluar_config(info["H"], sym, Rs, W_Z_RC, C_POST_RC, N_DFE, rng)
                bn.append(res["A"]["ber"]); bl.append(res["C"]["ber"]); bd.append(res["D"]["ber"])
            curvas[nb]["none"].append(bn); curvas[nb]["lin"].append(bl); curvas[nb]["dfe"].append(bd)
        print(f"  realización {r + 1:2d}/{N_REAL}")
    return canales, curvas


def main():
    os.makedirs(DIR_FIG, exist_ok=True)
    print("== C3+C4 — Barrido extendido (multi-umbral) ==")
    filas = []
    for criterio in ("equal_il", "fixed"):
        canales, curvas = correr(criterio)
        long_in = {nb: ("-" if canales[nb]["long_in"] is None
                        else f"{canales[nb]['long_in']:.1f}") for nb in canales}
        for nb in canales:
            for key, lab, _ in CONFIGS:
                B = np.array(curvas[nb][key])          # (n_real, n_rates)
                fila = {"criterio": criterio, "canal": nb, "long_pulg": long_in[nb],
                        "config": lab}
                for u in UMBRALES:
                    vs = np.array([viable_interp(B[r], GB, umbral=u) for r in range(N_REAL)])
                    sat = np.mean(vs >= GB[-1] - 1e-6) > 0.5      # mayoría saturada
                    med = float(vs.mean()); sd = float(vs.std())
                    tag = f"{med:.2f}±{sd:.2f}" + ("(sat)" if sat else "")
                    fila[f"viable_{u:.0e}"] = tag
                filas.append(fila)
                print(f"  {criterio:<9}{nb:<11}{lab:<14} " +
                      "  ".join(f"{u:.0e}:{fila[f'viable_{u:.0e}']}" for u in UMBRALES))

    cols = ["criterio", "canal", "long_pulg", "config"] + [f"viable_{u:.0e}" for u in UMBRALES]
    ruta_csv = os.path.join(DIR_FIG, "tabla_extendido_umbral.csv")
    with open(ruta_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader(); w.writerows(filas)
    print(f"\nCSV: {ruta_csv}")

    # Figura: sensibilidad de umbral en 'fixed' (la que tenía saturaciones)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        sub = [f for f in filas if f["criterio"] == "fixed"]
        canales_fixed = sorted(set(f["canal"] for f in sub),
                               key=lambda n: (not n.startswith("RC"), n))
        fig, ax = plt.subplots(figsize=(10, 5))
        x = np.arange(len(canales_fixed))
        w = 0.12
        for j, u in enumerate(UMBRALES):
            for c, (key, lab, _) in enumerate(CONFIGS):
                vals = []
                for nb in canales_fixed:
                    rec = next(f for f in sub if f["canal"] == nb and f["config"] == lab)
                    vals.append(float(rec[f"viable_{u:.0e}"].split("±")[0]))
                off = (j * len(CONFIGS) + c - 4.5) * w
                ax.bar(x + off, vals, w, label=f"{lab} @ {u:.0e}", alpha=0.85)
        ax.set_xticks(x); ax.set_xticklabels(canales_fixed)
        ax.set_ylabel("Tasa viable (GBaud)")
        ax.set_title("Sensibilidad de la tasa viable al umbral pre-FEC (criterio fixed)",
                     fontsize=11, fontweight="bold", color="#1a365d")
        ax.legend(fontsize=6, ncol=3); ax.grid(axis="y", alpha=0.3)
        fig.tight_layout()
        ruta_fig = os.path.join(DIR_FIG, "fig_sensibilidad_umbral.png")
        fig.savefig(ruta_fig, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Figura: {ruta_fig}")
    except Exception as e:
        print(f"(figura omitida: {e})")


if __name__ == "__main__":
    main()
