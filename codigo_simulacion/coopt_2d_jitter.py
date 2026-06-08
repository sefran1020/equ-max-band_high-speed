r"""
coopt_2d_jitter.py — Margen 2D (vertical x horizontal) bajo jitter (Fase γ).

Integra la Fase α (jitter / ojo horizontal) y la Fase β (co-optimización CTLE+DFE):
el objetivo ya no es la tasa viable vertical, sino la TASA VIABLE 2D, definida como
la tasa a la que el OJO HORIZONTAL CON JITTER cae por debajo de un mínimo (w_min UI).
Esta métrica es gain-invariante (eye width en UI) y penaliza el peaking del CTLE que
abre el ojo vertical pero cierra el horizontal y amplifica la sensibilidad al jitter.

Compara, por sustrato (criterio fixed), la tasa viable 2D bajo tres ecualizaciones:
  - EQ-RC   : heredada del RC.
  - β-opt   : óptimo vertical de la cadena con DFE (Fase β, coopt_resultados.json).
  - γ-opt   : óptimo del margen 2D (este script).

Entregables en ./figuras/:
  - fig_coopt2d_fixed.png   : tasa viable 2D (EQ-RC / β-opt / γ-opt)
  - fig_contorno_2d.png     : contorno vertical del ojo vs fase (EQ-RC vs γ-opt)
  - tabla_coopt2d.csv       : resumen con parámetros y márgenes
  - coopt2d_resultados.json : (w_z, c_post, A_DC) del óptimo 2D

Uso:  python coopt_2d_jitter.py
"""

import csv
import json
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cadena_enlace import (
    SPS, SIGMA_IN, W_Z_RC, C_POST_RC, CTLE_ADC,
    WZ_MIN, WZ_MAX, CPOST_MIN, CPOST_MAX, lhs_center,
    generar_prbs, modular_pam4, aplicar_ffe, aplicar_canal_freq,
    aplicar_ctle_freq, aplicar_rx_lpf,
)
from canales import construir_canales
from parametros import JITTER
from jitter import (
    banera_dfe, refinar, aplicar_jitter, eye_width,
)

DIR_FIG = os.path.join(os.path.dirname(__file__), "figuras")
DPI = 300

N_DFE = 2
ADC_MIN, ADC_MAX = 1.0, 3.0
W_MIN = 0.10                 # apertura horizontal mínima exigida (UI) bajo jitter
# Objetivo (DoE): rejilla gruesa (bañera DFE-aware -> N_SYM moderado por costo)
N_LHS = 16
N_SYM_OBJ = 550
TASAS_OBJ = np.linspace(6e9, 32e9, 10)
GB_OBJ = TASAS_OBJ / 1e9
# Reporte final
N_SYM_FIN = 600
N_REAL_FIN = 5
TASAS_FIN = np.linspace(6e9, 32e9, 14)
GB_FIN = TASAS_FIN / 1e9
SUBSTR = ("FR-4", "Megtron 6")
CRITERIO = "fixed"           # vista de diseño realista
RS_CONTORNO = 16e9           # tasa para el contorno 2D ilustrativo


def construir_C(H, sym, Rs, w_z, c_post, a_dc, rng):
    sps = SPS
    ts = 1.0 / (Rs * sps)
    fc = Rs
    ffe = aplicar_ffe(sym, 0.0, 1.0, c_post)
    ych = aplicar_canal_freq(np.repeat(ffe, sps), ts, H)
    nin = rng.normal(0, SIGMA_IN, len(ych))
    return aplicar_rx_lpf(aplicar_ctle_freq(ych + nin, ts, w_z, a_dc), ts, fc)


def ew_jit_rate(H, sym, Rs, w_z, c_post, a_dc, rng):
    """Apertura horizontal con jitter (UI) a la tasa Rs, sobre el ojo DFE-aware."""
    yC = construir_C(H, sym, Rs, w_z, c_post, a_dc, rng)
    fase, ber = banera_dfe(yC, sym, SPS, N_DFE)
    xf, bc = refinar(fase, ber)
    bj = aplicar_jitter(xf, bc, JITTER.sigma_rj_ui, JITTER.dj_pp_ui)
    return eye_width(xf, bj)


def cruce_decreciente(gbx, vals, umbral):
    """Tasa (GBaud) donde `vals` (decreciente) cruza `umbral`, interpolando."""
    vals = np.asarray(vals, dtype=float)
    above = vals >= umbral
    if not above.any():
        return 0.0
    if above.all():
        return float(gbx[-1])
    for i in range(len(gbx) - 1):
        if above[i] and not above[i + 1]:
            frac = (vals[i] - umbral) / (vals[i] - vals[i + 1] + 1e-12)
            return float(gbx[i] + frac * (gbx[i + 1] - gbx[i]))
    return float(gbx[above][-1])


def viable2d(H, sym, w_z, c_post, a_dc, tasas, gbx, seed=2026):
    rng = np.random.default_rng(seed)
    ews = [ew_jit_rate(H, sym, Rs, w_z, c_post, a_dc, rng) for Rs in tasas]
    return cruce_decreciente(gbx, ews, W_MIN)


def viable2d_reps(H, w_z, c_post, a_dc, bits_pool):
    vs = []
    for r in range(N_REAL_FIN):
        bsl = bits_pool[r * 2 * N_SYM_FIN: r * 2 * N_SYM_FIN + 2 * N_SYM_FIN + 200]
        sym = modular_pam4(bsl)[:N_SYM_FIN]
        rng = np.random.default_rng(1000 + r)
        ews = [ew_jit_rate(H, sym, Rs, w_z, c_post, a_dc, rng) for Rs in TASAS_FIN]
        vs.append(cruce_decreciente(GB_FIN, ews, W_MIN))
    return float(np.mean(vs)), float(np.std(vs))


def optimizar_2d(H):
    rng_lhs = np.random.default_rng(11)
    dis = lhs_center(3, N_LHS, rng_lhs)
    wz = 2 * np.pi * (WZ_MIN + dis[:, 0] * (WZ_MAX - WZ_MIN))
    cp = CPOST_MIN + dis[:, 1] * (CPOST_MAX - CPOST_MIN)
    adc = ADC_MIN + dis[:, 2] * (ADC_MAX - ADC_MIN)
    sym = modular_pam4(generar_prbs(2 * N_SYM_OBJ + 400, orden=23))[:N_SYM_OBJ]
    mejor = None
    for i in range(N_LHS):
        v = viable2d(H, sym, wz[i], cp[i], adc[i], TASAS_OBJ, GB_OBJ)
        if mejor is None or v > mejor["v"]:
            mejor = {"w_z": float(wz[i]), "c_post": float(cp[i]),
                     "a_dc": float(adc[i]), "v": float(v)}
    return mejor


def fig_barras(datos):
    nombres = list(datos.keys())
    x = np.arange(len(nombres))
    w = 0.26
    series = [("EQ-RC", "#e53e3e"), ("β-opt", "#dd6b20"), ("γ-opt (2D)", "#319795")]
    fig, ax = plt.subplots(figsize=(2.8 * len(nombres) + 2, 5.2))
    for j, (k, col) in enumerate(series):
        med = [datos[nb][k][0] for nb in nombres]
        sd = [datos[nb][k][1] for nb in nombres]
        ax.bar(x + (j - 1) * w, med, w, yerr=sd, capsize=4, color=col, edgecolor="k",
               alpha=0.9, label=k)
        for i in range(len(nombres)):
            ax.text(x[i] + (j - 1) * w, med[i] + sd[i] + 0.1, f"{med[i]:.1f}",
                    ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x); ax.set_xticklabels(nombres)
    ax.set_ylabel(f"Tasa viable 2D (GBaud)  [ojo horiz. ≥ {W_MIN:.2f} UI con jitter]")
    ax.set_title(f"Margen 2D bajo jitter — criterio '{CRITERIO}' ({N_REAL_FIN} realizaciones)",
                 fontsize=11, fontweight="bold", color="#1a365d")
    ax.grid(True, axis="y", alpha=0.3); ax.legend(fontsize=9)
    fig.tight_layout()
    ruta = os.path.join(DIR_FIG, "fig_coopt2d_fixed.png")
    fig.savefig(ruta, dpi=DPI, bbox_inches="tight"); plt.close(fig)
    return ruta


def fig_contorno(H, eqs):
    """Bañera DFE-aware (BER vs fase) a RS_CONTORNO, EQ-RC vs γ-opt, sin/con jitter."""
    sym = modular_pam4(generar_prbs(2 * 1500 + 200, orden=23))[:1500]
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    for lab, (wz, cp, adc), col in eqs:
        rng = np.random.default_rng(0)
        yC = construir_C(H, sym, RS_CONTORNO, wz, cp, adc, rng)
        fase, ber = banera_dfe(yC, sym, SPS, N_DFE)
        xf, bc = refinar(fase, ber)
        bj = aplicar_jitter(xf, bc, JITTER.sigma_rj_ui, JITTER.dj_pp_ui)
        ax.semilogy(xf, bc, color=col, lw=2, label=f"{lab} sin jitter")
        ax.semilogy(xf, bj, color=col, lw=2, ls="--", label=f"{lab} con jitter")
    ax.axhline(1e-2, color="gray", ls=":", lw=1)
    ax.set_title(f"Bañera DFE-aware (FR-4, {RS_CONTORNO/1e9:.0f} GBaud) — EQ-RC vs γ-opt",
                 fontsize=11, fontweight="bold", color="#1a365d")
    ax.set_xlabel("Fase de muestreo (UI)"); ax.set_ylabel("BER")
    ax.set_xlim(-0.5, 0.5); ax.set_ylim(1e-12, 1)
    ax.grid(True, which="both", alpha=0.25); ax.legend(fontsize=8)
    fig.tight_layout()
    ruta = os.path.join(DIR_FIG, "fig_contorno_2d.png")
    fig.savefig(ruta, dpi=DPI, bbox_inches="tight"); plt.close(fig)
    return ruta


def main():
    os.makedirs(DIR_FIG, exist_ok=True)
    print("== Fase γ — margen 2D (vertical x horizontal) bajo jitter ==")
    print(f"Objetivo: maximizar la tasa donde el ojo horizontal con jitter >= {W_MIN} UI "
          f"(RJ={JITTER.sigma_rj_ui:.3f}, DJ={JITTER.dj_pp_ui:.2f} UI)")

    canales, _ = construir_canales(CRITERIO)
    # Cargar β-opt (Fase β)
    ruta_beta = os.path.join(DIR_FIG, "coopt_resultados.json")
    beta = {}
    if os.path.exists(ruta_beta):
        with open(ruta_beta, encoding="utf-8") as fh:
            beta = json.load(fh).get(CRITERIO, {})

    bits_pool = generar_prbs(N_REAL_FIN * 2 * N_SYM_FIN + 5000, orden=23)
    datos, filas, params = {}, [], {}
    contorno_eqs = None
    for nb in SUBSTR:
        H = canales[nb]["H"]
        print(f"\n[{nb}] optimizando margen 2D (LHS {N_LHS})...")
        g = optimizar_2d(H)
        wzG = g["w_z"] / (2 * np.pi) / 1e9
        print(f"   γ-opt: w_z={wzG:.2f}GHz, c_post={g['c_post']:.3f}, A_DC={g['a_dc']:.2f}")

        # parámetros β (si existen); si no, usar EQ-RC como sustituto
        bo = beta.get(nb)
        if bo:
            wzB = 2 * np.pi * bo["w_z_GHz"] * 1e9
            cpB, adcB = bo["c_post"], bo["a_dc"]
        else:
            wzB, cpB, adcB = W_Z_RC, C_POST_RC, CTLE_ADC

        print(f"[{nb}] tasa viable 2D ({N_REAL_FIN} realiz.) x3 configs...")
        mRC = viable2d_reps(H, W_Z_RC, C_POST_RC, CTLE_ADC, bits_pool)
        mB = viable2d_reps(H, wzB, cpB, adcB, bits_pool)
        mG = viable2d_reps(H, g["w_z"], g["c_post"], g["a_dc"], bits_pool)
        datos[nb] = {"EQ-RC": mRC, "β-opt": mB, "γ-opt (2D)": mG}
        params[nb] = {"w_z_GHz": round(wzG, 2), "c_post": round(g["c_post"], 3),
                      "a_dc": round(g["a_dc"], 2)}
        filas.append({
            "canal": nb,
            "viable2d_EQ-RC": f"{mRC[0]:.1f}±{mRC[1]:.1f}",
            "viable2d_beta": f"{mB[0]:.1f}±{mB[1]:.1f}",
            "viable2d_gamma": f"{mG[0]:.1f}±{mG[1]:.1f}",
            "wz_GHz_gamma": round(wzG, 2), "A_DC_gamma": round(g["a_dc"], 2),
            "c_post_gamma": round(g["c_post"], 3),
        })
        if nb == "FR-4":
            contorno_eqs = [
                ("EQ-RC", (W_Z_RC, C_POST_RC, CTLE_ADC), "#e53e3e"),
                ("γ-opt", (g["w_z"], g["c_post"], g["a_dc"]), "#319795"),
            ]

    r1 = fig_barras(datos)
    r2 = fig_contorno(canales["FR-4"]["H"], contorno_eqs)

    cols = ["canal", "viable2d_EQ-RC", "viable2d_beta", "viable2d_gamma",
            "wz_GHz_gamma", "A_DC_gamma", "c_post_gamma"]
    anchos = {c: max(len(c), *(len(str(f[c])) for f in filas)) for c in cols}
    linea = " | ".join(c.ljust(anchos[c]) for c in cols)
    print("\n" + linea); print("-" * len(linea))
    for f in filas:
        print(" | ".join(str(f[c]).ljust(anchos[c]) for c in cols))

    ruta_csv = os.path.join(DIR_FIG, "tabla_coopt2d.csv")
    with open(ruta_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(filas)
    ruta_json = os.path.join(DIR_FIG, "coopt2d_resultados.json")
    with open(ruta_json, "w", encoding="utf-8") as fh:
        json.dump(params, fh, indent=2, ensure_ascii=False)

    print("\nEntregables generados:")
    for r in (r1, r2, ruta_csv, ruta_json):
        print(f"  - {r}")


if __name__ == "__main__":
    main()
