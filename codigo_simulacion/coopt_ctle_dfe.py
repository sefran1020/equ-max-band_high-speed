r"""
coopt_ctle_dfe.py — Co-optimización CONJUNTA CTLE+DFE (Fase β).

Resuelve el matiz del re-DoE (`bitacora` §6.1.4): el óptimo de FFE+CTLE (config C)
no es el óptimo de la cadena completa con DFE (config D). Aquí el objetivo del DoE
es la **tasa viable de D** (+DFE), y el espacio se amplía a $(\omega_z, c_{post}, A_{DC})$.

Compara, por sustrato y criterio, la tasa viable de D bajo tres ecualizaciones:
  - EQ-RC      : co-diseño heredado del RC (w_z=2.35 GHz, c_post=-0.301, A_DC=2.5).
  - C-opt      : optimiza C (FFE+CTLE), A_DC fijo=2.5 (re-DoE previo).
  - conjunta   : optimiza D (cadena completa) sobre (w_z, c_post, A_DC).

Hipótesis: la optimización conjunta usa un CTLE más suave (menor A_DC/w_z), deja la
ISI de cola al DFE, y da tasa viable de D mayor o igual, eliminando los casos D<C.

Entregables en ./figuras/:
  - fig_coopt_{equal_il,fixed}.png : tasa viable de D (EQ-RC / C-opt / conjunta)
  - tabla_coopt.csv                : resumen con parámetros óptimos
  - coopt_resultados.json          : (w_z, c_post, A_DC) conjuntos

Uso:  python coopt_ctle_dfe.py
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
    W_Z_RC, C_POST_RC, CTLE_ADC, WZ_MIN, WZ_MAX, CPOST_MIN, CPOST_MAX,
)
from canales import construir_canales

DIR_FIG = os.path.join(os.path.dirname(__file__), "figuras")
DPI = 300

N_DFE = 2
ADC_MIN, ADC_MAX = 1.0, 3.0          # rango de ganancia DC del CTLE (nuevo eje)
# DoE (objetivo): rejilla gruesa + pocos símbolos
N_LHS = 24
N_SYM_DOE = 700
TASAS_DOE = np.linspace(2e9, 30e9, 12)
GB_DOE = TASAS_DOE / 1e9
# Barrido final (reporte)
N_SYM_FIN = 1000
N_REAL_FIN = 15
TASAS = np.linspace(1e9, 32e9, 32)
GB = TASAS / 1e9

SUBSTR = ("FR-4", "Megtron 6")


def viable_de(H, sym, w_z, c_post, a_dc, clave, tasas, gbx, seed=2026):
    rng = np.random.default_rng(seed)
    bers = [evaluar_config(H, sym, Rs, w_z, c_post, N_DFE, rng, a_dc=a_dc)[clave]["ber"]
            for Rs in tasas]
    return viable_interp(bers, gbx)


def optimizar(H, objetivo):
    """Maximiza la tasa viable de `objetivo` ('C' o 'D').
    Para 'C', A_DC fijo=2.5 (2 variables); para 'D', A_DC libre (3 variables)."""
    rng_lhs = np.random.default_rng(7 if objetivo == "D" else 42)
    if objetivo == "D":
        dis = lhs_center(3, N_LHS, rng_lhs)
        wz = 2 * np.pi * (WZ_MIN + dis[:, 0] * (WZ_MAX - WZ_MIN))
        cp = CPOST_MIN + dis[:, 1] * (CPOST_MAX - CPOST_MIN)
        adc = ADC_MIN + dis[:, 2] * (ADC_MAX - ADC_MIN)
    else:
        dis = lhs_center(2, N_LHS, rng_lhs)
        wz = 2 * np.pi * (WZ_MIN + dis[:, 0] * (WZ_MAX - WZ_MIN))
        cp = CPOST_MIN + dis[:, 1] * (CPOST_MAX - CPOST_MIN)
        adc = np.full(N_LHS, CTLE_ADC)
    sym = modular_pam4(generar_prbs(2 * N_SYM_DOE + 400, orden=23))[:N_SYM_DOE]
    mejor = None
    for i in range(N_LHS):
        v = viable_de(H, sym, wz[i], cp[i], adc[i], objetivo, TASAS_DOE, GB_DOE)
        if mejor is None or v > mejor["v"]:
            mejor = {"w_z": float(wz[i]), "c_post": float(cp[i]),
                     "a_dc": float(adc[i]), "v": float(v)}
    return mejor


def viable30_D(H, w_z, c_post, a_dc):
    pool = generar_prbs(N_REAL_FIN * 2 * N_SYM_FIN + 5000, orden=23)
    vs = []
    for r in range(N_REAL_FIN):
        bsl = pool[r * 2 * N_SYM_FIN: r * 2 * N_SYM_FIN + 2 * N_SYM_FIN + 200]
        sym = modular_pam4(bsl)[:N_SYM_FIN]
        rng = np.random.default_rng(1000 + r)
        bers = [evaluar_config(H, sym, Rs, w_z, c_post, N_DFE, rng, a_dc=a_dc)["D"]["ber"]
                for Rs in TASAS]
        vs.append(viable_interp(bers, GB))
    return float(np.mean(vs)), float(np.std(vs))


def fig_coopt(criterio, datos):
    nombres = list(datos.keys())
    x = np.arange(len(nombres))
    w = 0.26
    series = [("EQ-RC", "EQ-RC", "#e53e3e"), ("C-opt", "C-opt", "#dd6b20"),
              ("conjunta", "joint", "#319795")]
    fig, ax = plt.subplots(figsize=(2.8 * len(nombres) + 2, 5.2))
    for j, (k, lab, col) in enumerate(series):
        med = [datos[nb][k][0] for nb in nombres]
        sd = [datos[nb][k][1] for nb in nombres]
        ax.bar(x + (j - 1) * w, med, w, yerr=sd, capsize=4, color=col, edgecolor="k",
               alpha=0.9, label=lab)
        for i in range(len(nombres)):
            ax.text(x[i] + (j - 1) * w, med[i] + sd[i] + 0.1, f"{med[i]:.1f}",
                    ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x); ax.set_xticklabels(nombres)
    ax.set_ylabel("Viable rate of FFE+CTLE+DFE (GBaud)  [BER ≤ 1e-2]")
    ax.grid(True, axis="y", alpha=0.3); ax.legend(fontsize=9)
    fig.tight_layout()
    ruta = os.path.join(DIR_FIG, f"fig_coopt_{criterio}.png")
    fig.savefig(ruta, dpi=DPI, bbox_inches="tight"); plt.close(fig)
    return ruta


def main():
    os.makedirs(DIR_FIG, exist_ok=True)
    print("== Fase β — co-optimización conjunta CTLE+DFE (objetivo = cadena con DFE) ==")
    print(f"Espacio: w_z [1,10] GHz, c_post [-0.40,-0.05], A_DC [{ADC_MIN},{ADC_MAX}]")

    filas, rutas, params = [], [], {}
    for criterio in ("equal_il", "fixed"):
        print(f"\n--- Criterio '{criterio}' ---")
        canales, _ = construir_canales(criterio)
        datos = {}
        params[criterio] = {}
        for nb in SUBSTR:
            H = canales[nb]["H"]
            print(f"  [{nb}] optimizando C (re-DoE) y D (conjunta)...")
            c_opt = optimizar(H, "C")
            d_opt = optimizar(H, "D")
            wzC = c_opt["w_z"] / (2 * np.pi) / 1e9
            wzD = d_opt["w_z"] / (2 * np.pi) / 1e9
            print(f"        C-opt:    w_z={wzC:.2f}GHz, c_post={c_opt['c_post']:.3f}, A_DC=2.50")
            print(f"        conjunta: w_z={wzD:.2f}GHz, c_post={d_opt['c_post']:.3f}, A_DC={d_opt['a_dc']:.2f}")

            print(f"  [{nb}] barrido final D ({N_REAL_FIN} realiz.) x3 configs...")
            mRC = viable30_D(H, W_Z_RC, C_POST_RC, CTLE_ADC)
            mC = viable30_D(H, c_opt["w_z"], c_opt["c_post"], CTLE_ADC)
            mD = viable30_D(H, d_opt["w_z"], d_opt["c_post"], d_opt["a_dc"])
            datos[nb] = {"EQ-RC": mRC, "C-opt": mC, "conjunta": mD}
            params[criterio][nb] = {
                "w_z_GHz": round(wzD, 2), "c_post": round(d_opt["c_post"], 3),
                "a_dc": round(d_opt["a_dc"], 2),
            }
            long_in = canales[nb]["long_in"]
            filas.append({
                "criterio": criterio, "canal": nb,
                "long_pulg": round(long_in, 1) if long_in else "-",
                "D_EQ-RC": f"{mRC[0]:.1f}±{mRC[1]:.1f}",
                "D_C-opt": f"{mC[0]:.1f}±{mC[1]:.1f}",
                "D_conjunta": f"{mD[0]:.1f}±{mD[1]:.1f}",
                "wz_GHz_conj": round(wzD, 2), "A_DC_conj": round(d_opt["a_dc"], 2),
                "c_post_conj": round(d_opt["c_post"], 3),
            })
        rutas.append(fig_coopt(criterio, datos))

    cols = ["criterio", "canal", "long_pulg", "D_EQ-RC", "D_C-opt", "D_conjunta",
            "wz_GHz_conj", "A_DC_conj", "c_post_conj"]
    anchos = {c: max(len(c), *(len(str(f[c])) for f in filas)) for c in cols}
    linea = " | ".join(c.ljust(anchos[c]) for c in cols)
    print("\n" + linea); print("-" * len(linea))
    for f in filas:
        print(" | ".join(str(f[c]).ljust(anchos[c]) for c in cols))

    ruta_csv = os.path.join(DIR_FIG, "tabla_coopt.csv")
    with open(ruta_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(filas)
    ruta_json = os.path.join(DIR_FIG, "coopt_resultados.json")
    with open(ruta_json, "w", encoding="utf-8") as fh:
        json.dump(params, fh, indent=2, ensure_ascii=False)

    print("\nEntregables generados:")
    for r in rutas + [ruta_csv, ruta_json]:
        print(f"  - {r}")


if __name__ == "__main__":
    main()
