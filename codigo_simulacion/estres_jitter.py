r"""
estres_jitter.py — Prueba de estrés de jitter (¿cuándo liga el eje horizontal?).

Barre el nivel de jitter (RJ y DJ crecientes, dual-Dirac) y mide la TASA VIABLE 2D
(ojo horizontal DFE-aware con jitter >= W_MIN UI) frente a la tasa viable vertical
(sin jitter). Cuando el jitter es suficiente, la tasa viable 2D cae por debajo de la
vertical: el eje horizontal (temporización) pasa a ligar. También comprueba si, en
ese régimen, la elección de ecualización (EQ-RC vs γ-opt) vuelve a importar.

Optimización: las bañeras DFE-aware (caras) se calculan UNA vez por (EQ, tasa,
realización); cada nivel de jitter solo re-aplica la convolución (barata).

Entregables en ./figuras/:
  - fig_estres_jitter.png : tasa viable 2D vs nivel de jitter (EQ-RC vs γ-opt)
  - tabla_estres_jitter.csv

Uso:  python estres_jitter.py
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
    generar_prbs, modular_pam4, aplicar_ffe, aplicar_canal_freq,
    aplicar_ctle_freq, aplicar_rx_lpf,
)
from canales import construir_canales
from jitter import banera_dfe, refinar, aplicar_jitter, eye_width

DIR_FIG = os.path.join(os.path.dirname(__file__), "figuras")
DPI = 300

N_DFE = 2
W_MIN = 0.10
N_SYM = 600
N_REAL = 5
TASAS = np.linspace(6e9, 32e9, 12)
GB = TASAS / 1e9
CRITERIO = "fixed"
SUBSTR = ("FR-4", "Megtron 6")
# Niveles de jitter: (RJ, DJ) = escala * (0.01, 0.05) UI ; escala 0 = sin jitter
# Rejilla fina cerca del acantilado (DJ 0.30-0.40 UI) para resolver el colapso.
ESCALAS = [0, 1, 2, 3, 4, 5, 6, 6.5, 7, 7.5, 8]
RJ0, DJ0 = 0.010, 0.050
COLORES = {"EQ-RC": "#e53e3e", "γ-opt": "#319795"}


def construir_C(H, sym, Rs, w_z, c_post, a_dc, rng):
    sps = SPS
    ts = 1.0 / (Rs * sps)
    fc = Rs
    ffe = aplicar_ffe(sym, 0.0, 1.0, c_post)
    ych = aplicar_canal_freq(np.repeat(ffe, sps), ts, H)
    nin = rng.normal(0, SIGMA_IN, len(ych))
    return aplicar_rx_lpf(aplicar_ctle_freq(ych + nin, ts, w_z, a_dc), ts, fc)


def cruce_decreciente(gbx, vals, umbral):
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


def bathtubs_clean(H, w_z, c_post, a_dc, bits_pool):
    """Bañeras DFE-aware finas (sin jitter) por realización y tasa. Devuelve (xf, [rep][rate]=bc)."""
    xf = None
    out = []
    for r in range(N_REAL):
        bsl = bits_pool[r * 2 * N_SYM: r * 2 * N_SYM + 2 * N_SYM + 200]
        sym = modular_pam4(bsl)[:N_SYM]
        rng = np.random.default_rng(1000 + r)
        per_rate = []
        for Rs in TASAS:
            yC = construir_C(H, sym, Rs, w_z, c_post, a_dc, rng)
            fase, ber = banera_dfe(yC, sym, SPS, N_DFE)
            xf, bc = refinar(fase, ber)
            per_rate.append(bc)
        out.append(per_rate)
    return xf, out


def viable2d(xf, bths, sigma_rj, dj_pp, modelo="dual_dirac"):
    vs = []
    for per_rate in bths:
        if sigma_rj <= 0 and dj_pp <= 0:
            ews = [eye_width(xf, bc) for bc in per_rate]
        else:
            ews = [eye_width(xf, aplicar_jitter(xf, bc, sigma_rj, dj_pp, modelo)) for bc in per_rate]
        vs.append(cruce_decreciente(GB, ews, W_MIN))
    return float(np.mean(vs)), float(np.std(vs))


def fig_modelo(xf, bths):
    """Compara el acantilado con jitter dual-Dirac vs PJ senoidal (FR-4, EQ-RC)."""
    dj = [DJ0 * k for k in ESCALAS]
    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    for modelo, ls, lab in [("dual_dirac", "-", "Dual-Dirac (hard DJ, worst case)"),
                            ("pj", "--", "Sinusoidal PJ (arcsine)")]:
        med = [viable2d(xf, bths, RJ0 * k, DJ0 * k, modelo)[0] for k in ESCALAS]
        ax.plot(dj, med, ls, marker="o", lw=2, label=lab)
    ax.set_xlabel("Peak-to-peak DJ (UI)  [RJ = DJ/5]")
    ax.set_ylabel(f"2D viable rate (GBaud)  [eye ≥ {W_MIN} UI]")
    ax.grid(True, alpha=0.3); ax.legend(fontsize=9)
    fig.tight_layout()
    ruta = os.path.join(DIR_FIG, "fig_estres_modelo.png")
    fig.savefig(ruta, dpi=DPI, bbox_inches="tight"); plt.close(fig)
    return ruta


def main():
    os.makedirs(DIR_FIG, exist_ok=True)
    print("== Prueba de estrés de jitter (¿cuándo liga el eje horizontal?) ==")
    print(f"Niveles (RJ,DJ) = escala·({RJ0},{DJ0}) UI; escalas {ESCALAS}")

    canales, _ = construir_canales(CRITERIO)
    ruta_g = os.path.join(DIR_FIG, "coopt2d_resultados.json")
    gamma = {}
    if os.path.exists(ruta_g):
        with open(ruta_g, encoding="utf-8") as fh:
            gamma = json.load(fh)

    bits_pool = generar_prbs(N_REAL * 2 * N_SYM + 5000, orden=23)
    resultados = {}   # resultados[nb][eq] = lista de (media,desv) por escala
    filas = []
    fr4_bths = None   # bañeras FR-4/EQ-RC para comparar modelos de jitter
    for nb in SUBSTR:
        H = canales[nb]["H"]
        eqs = {"EQ-RC": (W_Z_RC, C_POST_RC, CTLE_ADC)}
        g = gamma.get(nb)
        if g:
            eqs["γ-opt"] = (2 * np.pi * g["w_z_GHz"] * 1e9, g["c_post"], g["a_dc"])
        resultados[nb] = {}
        for eqname, (wz, cp, adc) in eqs.items():
            print(f"\n[{nb} / {eqname}] precomputando bañeras DFE-aware...")
            xf, bths = bathtubs_clean(H, wz, cp, adc, bits_pool)
            if nb == "FR-4" and eqname == "EQ-RC":
                fr4_bths = (xf, bths)
            serie = []
            for k in ESCALAS:
                m, s = viable2d(xf, bths, RJ0 * k, DJ0 * k)
                serie.append((m, s))
                filas.append({"canal": nb, "eq": eqname, "escala": k,
                              "DJ_pp_UI": round(DJ0 * k, 3), "RJ_UI": round(RJ0 * k, 3),
                              "viable2d_GBaud": round(m, 2), "desv": round(s, 2)})
            resultados[nb][eqname] = serie
            print(f"   {eqname}: " + " | ".join(f"k{ESCALAS[j]}:{serie[j][0]:.1f}" for j in range(len(ESCALAS))))

    # Figura
    dj = [DJ0 * k for k in ESCALAS]
    fig, axes = plt.subplots(1, len(SUBSTR), figsize=(6.2 * len(SUBSTR), 4.6), sharey=True)
    for ax, nb in zip(np.atleast_1d(axes), SUBSTR):
        for eqname, serie in resultados[nb].items():
            med = [v[0] for v in serie]
            sd = [v[1] for v in serie]
            ax.errorbar(dj, med, yerr=sd, marker="o", lw=2, capsize=3,
                        color=COLORES.get(eqname, "k"), label=eqname)
        ref = resultados[nb]["EQ-RC"][0][0]   # k=0 (sin jitter) ~ tasa viable vertical
        ax.axhline(ref, color="0.5", ls="--", lw=1, label="vertical viable (no jitter)")
        ax.set_title(nb, fontsize=11, fontweight="bold")
        ax.set_xlabel("Peak-to-peak DJ (UI)  [RJ = DJ/5]", fontsize=9)
        ax.grid(True, alpha=0.3)
    np.atleast_1d(axes)[0].set_ylabel(f"2D viable rate (GBaud)  [eye ≥ {W_MIN} UI]", fontsize=9)
    np.atleast_1d(axes)[0].legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    ruta = os.path.join(DIR_FIG, "fig_estres_jitter.png")
    fig.savefig(ruta, dpi=DPI, bbox_inches="tight"); plt.close(fig)

    ruta_mod = fig_modelo(fr4_bths[0], fr4_bths[1]) if fr4_bths else None

    cols = ["canal", "eq", "escala", "DJ_pp_UI", "RJ_UI", "viable2d_GBaud", "desv"]
    anchos = {c: max(len(c), *(len(str(f[c])) for f in filas)) for c in cols}
    linea = " | ".join(c.ljust(anchos[c]) for c in cols)
    print("\n" + linea); print("-" * len(linea))
    for f in filas:
        print(" | ".join(str(f[c]).ljust(anchos[c]) for c in cols))

    ruta_csv = os.path.join(DIR_FIG, "tabla_estres_jitter.csv")
    with open(ruta_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(filas)

    print("\nEntregables generados:")
    print(f"  - {ruta}")
    if ruta_mod:
        print(f"  - {ruta_mod}")
    print(f"  - {ruta_csv}")


if __name__ == "__main__":
    main()
