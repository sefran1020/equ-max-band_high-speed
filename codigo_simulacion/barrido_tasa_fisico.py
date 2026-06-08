"""
barrido_tasa_fisico.py — Tasa viable (BER vs GBaud) sobre el canal FÍSICO.

Replica la metodología de robustez del notebook (cell 29): N realizaciones, cada
una con un PATRÓN PRBS distinto + RUIDO distinto, y reporta la tasa de símbolo
viable (BER <= 1e-2, pre-FEC) INTERPOLADA como media ± desviación.

Compara, a igualdad de co-diseño heredado del RC (w_z=2π·2.35 GHz, c_post=-0.301,
DFE-ZF 2 taps), la tasa viable dispersiva (FR-4, Megtron 6) contra la del canal
RC del artículo (referencia checkpoint: 2.79 / 9.56 / 14.11 GBaud para
sin EQ / FFE+CTLE / +DFE).

Genera DOS vistas (decisión del usuario):
  - 'equal_il' : cada línea recortada a la misma pérdida del RC en Nyquist.
  - 'fixed'    : ambas líneas a longitud física fija (10 in).

Entregables en ./figuras/:
  - fig_barrido_tasa_equal_il.png / _fixed.png  : BER vs GBaud (mediana + envolvente)
  - fig_tasa_viable_equal_il.png  / _fixed.png  : tasa viable media ± desv. (barras)
  - tabla_tasa_viable.csv                        : resumen numérico

Uso:
    python barrido_tasa_fisico.py
"""

import csv
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cadena_enlace import (
    generar_prbs, modular_pam4, evaluar_config, viable_interp,
    W_Z_RC, C_POST_RC,
)
from canales import construir_canales

DIR_FIG = os.path.join(os.path.dirname(__file__), "figuras")
DPI = 300

N_REAL = 30                                  # realizaciones (patrón + ruido)
N_SYM = 1000                                 # símbolos por realización
N_DFE = 2                                    # taps de DFE (fiel al artículo)
TASAS = np.linspace(1e9, 32e9, 32)   # hasta 32 GBaud para capturar todos los cruces
GB = TASAS / 1e9
UMBRAL = 1e-2

CONFIGS = [("none", "No EQ", "#e53e3e"),
           ("lin", "FFE+CTLE", "#dd6b20"),
           ("dfe", "FFE+CTLE+DFE", "#319795")]
CHK_RC = {"none": 2.79, "lin": 9.56, "dfe": 14.11}   # referencia checkpoint (RC)


def correr_criterio(criterio):
    canales, il_obj = construir_canales(criterio)
    bits_pool = generar_prbs(N_REAL * 2 * N_SYM + 5000, orden=23)
    curvas = {nb: {"none": [], "lin": [], "dfe": []} for nb in canales}
    viables = {nb: {"none": [], "lin": [], "dfe": []} for nb in canales}

    print(f"\n== Criterio '{criterio}' ({N_REAL} realizaciones) ==")
    for r in range(N_REAL):
        bsl = bits_pool[r * 2 * N_SYM: r * 2 * N_SYM + 2 * N_SYM + 200]
        sym = modular_pam4(bsl)[:N_SYM]
        for nb, info in canales.items():
            rng = np.random.default_rng(1000 + r)   # ruido pareado entre canales
            bn, bl, bd = [], [], []
            for Rs in TASAS:
                res = evaluar_config(info["H"], sym, Rs, W_Z_RC, C_POST_RC, N_DFE, rng)
                bn.append(res["A"]["ber"]); bl.append(res["C"]["ber"]); bd.append(res["D"]["ber"])
            curvas[nb]["none"].append(bn); curvas[nb]["lin"].append(bl); curvas[nb]["dfe"].append(bd)
            viables[nb]["none"].append(viable_interp(bn, GB))
            viables[nb]["lin"].append(viable_interp(bl, GB))
            viables[nb]["dfe"].append(viable_interp(bd, GB))
        print(f"  realización {r + 1:2d}/{N_REAL} completada")
    return canales, il_obj, curvas, viables


def fig_ber_vs_tasa(criterio, canales, curvas):
    nombres = list(canales.keys())
    fig, axes = plt.subplots(1, len(nombres), figsize=(5.2 * len(nombres), 4.6), sharey=True)
    for ax, nb in zip(np.atleast_1d(axes), nombres):
        for key, lab, col in CONFIGS:
            B = np.clip(np.array(curvas[nb][key]), 1e-12, 1.0)   # (n_real, n_rates)
            ax.semilogy(GB, np.median(B, 0), "-", color=col, label=lab)
            ax.fill_between(GB, B.min(0), B.max(0), color=col, alpha=0.18)
        ax.axhline(UMBRAL, color="gray", ls="--", alpha=0.7)
        long_txt = "" if canales[nb]["long_in"] is None else f"\n({canales[nb]['long_in']:.1f} in)"
        ax.set_title(nb + long_txt, fontsize=10, fontweight="bold")
        ax.set_xlabel("Rate (GBaud)", fontsize=9)
        ax.grid(True, which="both", alpha=0.25)
        ax.set_ylim(1e-12, 1)
    np.atleast_1d(axes)[0].set_ylabel("BER (median, envelope)", fontsize=9)
    np.atleast_1d(axes)[0].legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    ruta = os.path.join(DIR_FIG, f"fig_barrido_tasa_{criterio}.png")
    fig.savefig(ruta, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return ruta


def fig_tasa_viable(criterio, canales, viables):
    nombres = list(canales.keys())
    x = np.arange(len(nombres))
    w = 0.26
    fig, ax = plt.subplots(figsize=(2.6 * len(nombres) + 2, 5.2))
    for j, (key, lab, col) in enumerate(CONFIGS):
        med = [np.mean(viables[nb][key]) for nb in nombres]
        sd = [np.std(viables[nb][key]) for nb in nombres]
        ax.bar(x + (j - 1) * w, med, w, yerr=sd, capsize=4, color=col, edgecolor="k",
               alpha=0.88, label=lab)
        for i in range(len(nombres)):
            ax.text(x[i] + (j - 1) * w, med[i] + sd[i] + 0.1, f"{med[i]:.1f}",
                    ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x); ax.set_xticklabels(nombres)
    ax.set_ylabel("Tasa de símbolo viable (GBaud)  [BER ≤ 1e-2]")
    ax.set_title(f"Tasa viable interpolada: media ± desv. sobre {N_REAL} realizaciones "
                 f"— criterio '{criterio}'", fontsize=11, fontweight="bold", color="#1a365d")
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    ruta = os.path.join(DIR_FIG, f"fig_tasa_viable_{criterio}.png")
    fig.savefig(ruta, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return ruta


def main():
    os.makedirs(DIR_FIG, exist_ok=True)
    print("== Recomendación 01 — Barrido de tasa sobre canal físico (30 realizaciones) ==")
    print(f"Co-diseño heredado del RC: w_z=2π·2.35GHz, c_post={C_POST_RC}, DFE {N_DFE} taps")
    print(f"Referencia RC (checkpoint): sin EQ {CHK_RC['none']} | FFE+CTLE {CHK_RC['lin']} | "
          f"+DFE {CHK_RC['dfe']} GBaud")

    filas = []
    rutas = []
    for criterio in ("equal_il", "fixed"):
        canales, il_obj, curvas, viables = correr_criterio(criterio)
        rutas.append(fig_ber_vs_tasa(criterio, canales, curvas))
        rutas.append(fig_tasa_viable(criterio, canales, viables))

        for nb in canales:
            long_txt = "-" if canales[nb]["long_in"] is None else f"{canales[nb]['long_in']:.1f}"
            for key, lab, _ in CONFIGS:
                v = np.array(viables[nb][key])
                filas.append({
                    "criterio": criterio, "canal": nb, "long_pulg": long_txt,
                    "config": lab,
                    "tasa_viable_GBaud_media": round(float(v.mean()), 2),
                    "desv": round(float(v.std()), 2),
                    "ref_RC_checkpoint": CHK_RC[key] if nb.startswith("RC") else "-",
                })

    # Tabla por pantalla
    cols = ["criterio", "canal", "long_pulg", "config",
            "tasa_viable_GBaud_media", "desv", "ref_RC_checkpoint"]
    anchos = {c: max(len(c), *(len(str(f[c])) for f in filas)) for c in cols}
    linea = " | ".join(c.ljust(anchos[c]) for c in cols)
    print("\n" + linea); print("-" * len(linea))
    crit_prev = None
    for f in filas:
        if f["criterio"] != crit_prev:
            print("-" * len(linea)); crit_prev = f["criterio"]
        print(" | ".join(str(f[c]).ljust(anchos[c]) for c in cols))

    ruta_csv = os.path.join(DIR_FIG, "tabla_tasa_viable.csv")
    with open(ruta_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader(); w.writerows(filas)

    print("\nEntregables generados:")
    for r in rutas + [ruta_csv]:
        print(f"  - {r}")


if __name__ == "__main__":
    main()
