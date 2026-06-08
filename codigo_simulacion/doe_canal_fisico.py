"""
doe_canal_fisico.py — Re-optimización del co-diseño (FFE+CTLE) sobre el canal FÍSICO.

Motivación: el co-diseño fijado en el checkpoint (ω_z=2π·2.35 GHz, c_post=−0.301)
se optimizó sobre el canal RC. El barrido mostró que NO es óptimo para los canales
dispersivos (FR-4 / Megtron 6), cuyo roll-off (skin+dieléctrico) es distinto. Aquí
se re-optimiza (ω_z, c_post) POR SUBSTRATO y se compara la tasa viable resultante
contra la EQ heredada del RC, con 30 realizaciones (patrón PRBS + ruido).

Procedimiento:
  1) DoE-LHS sobre (ω_z, c_post): para cada muestra, mini-barrido de tasa y se mide
     la tasa viable de FFE+CTLE (BER≤1e-2). Se elige la que la maximiza.
  2) Con la EQ re-optimizada, barrido completo (32 GBaud, 30 realizaciones) para
     FFE+CTLE y +DFE; se compara contra la EQ del RC.

Entregables en ./figuras/:
  - fig_redoe_{equal_il,fixed}.png : tasa viable EQ-RC vs EQ-reoptimizada (barras)
  - tabla_redoe.csv                 : resumen numérico
  - doe_resultados.json             : (ω_z, c_post) re-optimizados por substrato/criterio

Uso:  python doe_canal_fisico.py
"""

import csv
import json
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cadena_enlace import (
    generar_prbs, modular_pam4, evaluar_config, viable_interp, lhs_center,
    W_Z_RC, C_POST_RC, WZ_MIN, WZ_MAX, CPOST_MIN, CPOST_MAX,
)
from canales import construir_canales

DIR_FIG = os.path.join(os.path.dirname(__file__), "figuras")
DPI = 300

N_DFE = 2
# DoE (objetivo): rejilla de tasa más gruesa y menos símbolos -> rápido
N_LHS = 30
N_SYM_DOE = 800
TASAS_DOE = np.linspace(2e9, 26e9, 13)
GB_DOE = TASAS_DOE / 1e9
# Barrido final (reporte): fino, 30 realizaciones
N_SYM_FIN = 1000
N_REAL_FIN = 30
TASAS = np.linspace(1e9, 32e9, 32)
GB = TASAS / 1e9

SUBSTR = ("FR-4", "Megtron 6")     # se re-optimiza solo el canal dispersivo


# --------------------------------------------------------------------------- #
# 1) Re-DoE: maximizar la tasa viable de FFE+CTLE
# --------------------------------------------------------------------------- #
def viable_de_params(H, sym, w_z, c_post, key, tasas, gbx, seed=2026):
    rng = np.random.default_rng(seed)
    bers = [evaluar_config(H, sym, Rs, w_z, c_post, N_DFE, rng)[key]["ber"] for Rs in tasas]
    return viable_interp(bers, gbx)


def re_doe(H):
    """LHS sobre (w_z, c_post); maximiza la tasa viable de FFE+CTLE (config 'C')."""
    rng_lhs = np.random.default_rng(42)
    dis = lhs_center(2, N_LHS, rng_lhs)
    wz = 2 * np.pi * (WZ_MIN + dis[:, 0] * (WZ_MAX - WZ_MIN))
    cp = CPOST_MIN + dis[:, 1] * (CPOST_MAX - CPOST_MIN)
    sym = modular_pam4(generar_prbs(2 * N_SYM_DOE + 400, orden=23))[:N_SYM_DOE]

    mejor = None
    for i in range(N_LHS):
        v = viable_de_params(H, sym, wz[i], cp[i], "C", TASAS_DOE, GB_DOE)
        if mejor is None or v > mejor["viable"]:
            mejor = {"w_z": float(wz[i]), "c_post": float(cp[i]), "viable": float(v)}
    return mejor


# --------------------------------------------------------------------------- #
# 2) Barrido final con 30 realizaciones (devuelve viable de C y D)
# --------------------------------------------------------------------------- #
def viable_30(H, w_z, c_post):
    bits_pool = generar_prbs(N_REAL_FIN * 2 * N_SYM_FIN + 5000, orden=23)
    vC, vD = [], []
    for r in range(N_REAL_FIN):
        bsl = bits_pool[r * 2 * N_SYM_FIN: r * 2 * N_SYM_FIN + 2 * N_SYM_FIN + 200]
        sym = modular_pam4(bsl)[:N_SYM_FIN]
        rng = np.random.default_rng(1000 + r)
        bC, bD = [], []
        for Rs in TASAS:
            res = evaluar_config(H, sym, Rs, w_z, c_post, N_DFE, rng)
            bC.append(res["C"]["ber"]); bD.append(res["D"]["ber"])
        vC.append(viable_interp(bC, GB)); vD.append(viable_interp(bD, GB))
    return (float(np.mean(vC)), float(np.std(vC)),
            float(np.mean(vD)), float(np.std(vD)))


# --------------------------------------------------------------------------- #
def fig_comparacion(criterio, datos):
    """datos[nb] = dict con medias/desv de C y D para EQ-RC y EQ-reopt."""
    nombres = list(datos.keys())
    x = np.arange(len(nombres))
    w = 0.2
    series = [
        ("C_rc", "FFE+CTLE (EQ-RC)", "#fbd38d"),
        ("C_re", "FFE+CTLE (re-opt.)", "#dd6b20"),
        ("D_rc", "+DFE (EQ-RC)", "#9decd9"),
        ("D_re", "+DFE (re-opt.)", "#319795"),
    ]
    fig, ax = plt.subplots(figsize=(2.8 * len(nombres) + 3, 5.4))
    for j, (k, lab, col) in enumerate(series):
        med = [datos[nb][k][0] for nb in nombres]
        sd = [datos[nb][k][1] for nb in nombres]
        ax.bar(x + (j - 1.5) * w, med, w, yerr=sd, capsize=3, color=col,
               edgecolor="k", alpha=0.9, label=lab)
        for i in range(len(nombres)):
            ax.text(x[i] + (j - 1.5) * w, med[i] + sd[i] + 0.15, f"{med[i]:.1f}",
                    ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x); ax.set_xticklabels(nombres)
    ax.set_ylabel("Viable rate (GBaud)  [BER ≤ 1e-2]")
    ax.grid(True, axis="y", alpha=0.3); ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    ruta = os.path.join(DIR_FIG, f"fig_redoe_{criterio}.png")
    fig.savefig(ruta, dpi=DPI, bbox_inches="tight"); plt.close(fig)
    return ruta


def main():
    os.makedirs(DIR_FIG, exist_ok=True)
    print("== Re-DoE del co-diseño (FFE+CTLE) sobre el canal físico ==")
    print(f"EQ heredada del RC: w_z=2π·2.35GHz, c_post={C_POST_RC}")

    filas, rutas = [], []
    params_out = {}
    for criterio in ("equal_il", "fixed"):
        print(f"\n--- Criterio '{criterio}' ---")
        canales, _ = construir_canales(criterio)
        datos = {}
        params_out[criterio] = {}
        for nb in SUBSTR:
            H = canales[nb]["H"]
            print(f"  [{nb}] re-optimizando (LHS {N_LHS})...")
            best = re_doe(H)
            wz_ghz = best["w_z"] / (2 * np.pi) / 1e9
            print(f"        re-opt: w_z=2π·{wz_ghz:.2f}GHz, c_post={best['c_post']:.3f}")

            print(f"  [{nb}] barrido final EQ-RC ({N_REAL_FIN} realiz.)...")
            cC_rc, sC_rc, dD_rc, sD_rc = viable_30(H, W_Z_RC, C_POST_RC)
            print(f"  [{nb}] barrido final EQ-reopt ({N_REAL_FIN} realiz.)...")
            cC_re, sC_re, dD_re, sD_re = viable_30(H, best["w_z"], best["c_post"])

            datos[nb] = {
                "C_rc": (cC_rc, sC_rc), "C_re": (cC_re, sC_re),
                "D_rc": (dD_rc, sD_rc), "D_re": (dD_re, sD_re),
            }
            params_out[criterio][nb] = {
                "w_z": best["w_z"], "w_z_GHz": wz_ghz, "c_post": best["c_post"],
            }
            long_in = canales[nb]["long_in"]
            filas.append({
                "criterio": criterio, "canal": nb,
                "long_pulg": round(long_in, 1) if long_in else "-",
                "wz_GHz_reopt": round(wz_ghz, 2), "c_post_reopt": round(best["c_post"], 3),
                "FFE+CTLE_EQ-RC": f"{cC_rc:.1f}±{sC_rc:.1f}",
                "FFE+CTLE_reopt": f"{cC_re:.1f}±{sC_re:.1f}",
                "+DFE_EQ-RC": f"{dD_rc:.1f}±{sD_rc:.1f}",
                "+DFE_reopt": f"{dD_re:.1f}±{sD_re:.1f}",
            })
        rutas.append(fig_comparacion(criterio, datos))

    # Tabla
    cols = ["criterio", "canal", "long_pulg", "wz_GHz_reopt", "c_post_reopt",
            "FFE+CTLE_EQ-RC", "FFE+CTLE_reopt", "+DFE_EQ-RC", "+DFE_reopt"]
    anchos = {c: max(len(c), *(len(str(f[c])) for f in filas)) for c in cols}
    linea = " | ".join(c.ljust(anchos[c]) for c in cols)
    print("\n" + linea); print("-" * len(linea))
    for f in filas:
        print(" | ".join(str(f[c]).ljust(anchos[c]) for c in cols))

    ruta_csv = os.path.join(DIR_FIG, "tabla_redoe.csv")
    with open(ruta_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(filas)
    ruta_json = os.path.join(DIR_FIG, "doe_resultados.json")
    with open(ruta_json, "w", encoding="utf-8") as fh:
        json.dump(params_out, fh, indent=2, ensure_ascii=False)

    print("\nEntregables generados:")
    for r in rutas + [ruta_csv, ruta_json]:
        print(f"  - {r}")


if __name__ == "__main__":
    main()
