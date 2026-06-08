"""
barrido_ojo_horizontal.py — Apertura HORIZONTAL del ojo con jitter (Fase α).

Calcula la curva de bañera (BER vs fase de muestreo) sobre la forma de onda
ecualizada linealmente (entrada del decisor) y la apertura horizontal del ojo
(eye width), sin y con jitter (RJ + DJ, dual-Dirac), para RC N=5, FR-4 y
Megtron 6 (criterio equal_il) en el punto de operación.

Entregables en ./figuras/:
  - fig_banera_canales.png   : bañera (escenario C) por canal, sin/con jitter
  - fig_banera_escenarios.png: bañera FR-4 escenarios A/B/C (apertura horizontal por etapa)
  - tabla_ojo_horizontal.csv : eye width (UI y ps) sin/con jitter

Uso:  python barrido_ojo_horizontal.py
"""

import csv
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cadena_enlace import (
    SPS, SIGMA_IN, W_Z_RC, C_POST_RC,
    generar_prbs, modular_pam4, aplicar_ffe, aplicar_canal_freq,
    aplicar_ctle_freq, aplicar_rx_lpf,
)
from canales import construir_canales
from parametros import JITTER, OPERACION
from jitter import banera_clean, refinar, aplicar_jitter, eye_width, tj_dual_dirac

DIR_FIG = os.path.join(os.path.dirname(__file__), "figuras")
DPI = 300

RS = OPERACION.baud_op       # 4 GBaud
UI_PS = 1e12 / RS            # 250 ps
N_SYM = 3000
N_REAL = 5
N_FINE = 512
SEL = ["RC N=5", "FR-4", "Megtron 6"]
COLORES = {"RC N=5": "k", "FR-4": "#2b6cb0", "Megtron 6": "#dd6b20"}


def construir_ABC(H, sym, rng):
    sps = SPS
    ts = 1.0 / (RS * sps)
    fc = RS
    y_a = aplicar_canal_freq(np.repeat(sym, sps), ts, H)
    yA = aplicar_rx_lpf(y_a + rng.normal(0, SIGMA_IN, len(y_a)), ts, fc)
    ffe = aplicar_ffe(sym, 0.0, 1.0, C_POST_RC)
    y_chan = aplicar_canal_freq(np.repeat(ffe, sps), ts, H)
    yB = aplicar_rx_lpf(y_chan + rng.normal(0, SIGMA_IN, len(y_chan)), ts, fc)
    nin = rng.normal(0, SIGMA_IN, len(y_chan))
    yC = aplicar_rx_lpf(aplicar_ctle_freq(y_chan + nin, ts, W_Z_RC), ts, fc)
    return {"A": yA, "B": yB, "C": yC}


def banera_promedio(H, escenario, bits_pool):
    """Bañera CLEAN promediada (en log10 BER) sobre N_REAL realizaciones."""
    acc = None
    fase = None
    for r in range(N_REAL):
        bsl = bits_pool[r * 2 * N_SYM: r * 2 * N_SYM + 2 * N_SYM + 200]
        sym = modular_pam4(bsl)[:N_SYM]
        rng = np.random.default_rng(2000 + r)
        w = construir_ABC(H, sym, rng)[escenario]
        f, b = banera_clean(w, sym)
        lb = np.log10(b)
        acc = lb if acc is None else acc + lb
        fase = f
    return fase, np.power(10.0, acc / N_REAL)


def procesar(H, escenario, bits_pool):
    """Devuelve (xf, ber_clean, ber_jit, ew_clean_UI, ew_jit_UI)."""
    fase, ber = banera_promedio(H, escenario, bits_pool)
    xf, bc = refinar(fase, ber, N_FINE)
    bj = aplicar_jitter(xf, bc, JITTER.sigma_rj_ui, JITTER.dj_pp_ui)
    return xf, bc, bj, eye_width(xf, bc), eye_width(xf, bj)


def fig_canales(res):
    fig, axes = plt.subplots(1, len(SEL), figsize=(5 * len(SEL), 4.2), sharey=True)
    for ax, nb in zip(np.atleast_1d(axes), SEL):
        xf, bc, bj, ewc, ewj = res[nb]["C"]
        ax.semilogy(xf, bc, color=COLORES[nb], lw=2, label="no jitter")
        ax.semilogy(xf, bj, color=COLORES[nb], lw=2, ls="--", label="with jitter")
        ax.axhline(1e-2, color="gray", ls=":", lw=1)
        ax.set_title(f"{nb}\neye width: {ewc:.3f} → {ewj:.3f} UI", fontsize=10, fontweight="bold")
        ax.set_xlabel("Sampling phase (UI)", fontsize=9)
        ax.set_xlim(-0.5, 0.5); ax.set_ylim(1e-12, 1)
        ax.grid(True, which="both", alpha=0.25)
    np.atleast_1d(axes)[0].set_ylabel("BER", fontsize=9)
    np.atleast_1d(axes)[0].legend(fontsize=8, loc="upper center")
    fig.tight_layout()
    ruta = os.path.join(DIR_FIG, "fig_banera_canales.png")
    fig.savefig(ruta, dpi=DPI, bbox_inches="tight"); plt.close(fig)
    return ruta


def fig_escenarios(res, nb="FR-4"):
    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    for esc, col, lab in [("A", "#e53e3e", "A — no EQ"),
                          ("B", "#dd6b20", "B — FFE"),
                          ("C", "#319795", "C — FFE+CTLE")]:
        xf, bc, _, ewc, _ = res[nb][esc]
        ax.semilogy(xf, bc, color=col, lw=2, label=f"{lab}  (width {ewc:.2f} UI)")
    ax.axhline(1e-2, color="gray", ls=":", lw=1, label="pre-FEC threshold")
    ax.set_xlabel("Sampling phase (UI)"); ax.set_ylabel("BER")
    ax.set_xlim(-0.5, 0.5); ax.set_ylim(1e-12, 1)
    ax.grid(True, which="both", alpha=0.25); ax.legend(fontsize=8)
    fig.tight_layout()
    ruta = os.path.join(DIR_FIG, "fig_banera_escenarios.png")
    fig.savefig(ruta, dpi=DPI, bbox_inches="tight"); plt.close(fig)
    return ruta


def main():
    os.makedirs(DIR_FIG, exist_ok=True)
    print("== Fase α — apertura horizontal del ojo (bañera + jitter) ==")
    print(f"Operación {RS/1e9:.0f} GBaud (UI={UI_PS:.0f} ps) | "
          f"jitter dual-Dirac: RJ={JITTER.sigma_rj_ui:.3f} UI rms, DJ={JITTER.dj_pp_ui:.2f} UI pp")

    canales, _ = construir_canales("equal_il")
    bits_pool = generar_prbs(N_REAL * 2 * N_SYM + 5000, orden=23)

    res = {}
    for nb in SEL:
        H = canales[nb]["H"]
        escenarios = ["A", "B", "C"] if nb == "FR-4" else ["C"]
        res[nb] = {e: procesar(H, e, bits_pool) for e in escenarios}

    r1 = fig_canales(res)
    r2 = fig_escenarios(res, "FR-4")

    # Tabla
    filas = []
    for nb in SEL:
        for esc in res[nb]:
            _, _, _, ewc, ewj = res[nb][esc]
            filas.append({
                "canal": nb, "escenario": esc,
                "eye_width_UI_sinJ": round(ewc, 3), "eye_width_ps_sinJ": round(ewc * UI_PS, 1),
                "eye_width_UI_conJ": round(ewj, 3), "eye_width_ps_conJ": round(ewj * UI_PS, 1),
            })
    cols = ["canal", "escenario", "eye_width_UI_sinJ", "eye_width_ps_sinJ",
            "eye_width_UI_conJ", "eye_width_ps_conJ"]
    anchos = {c: max(len(c), *(len(str(f[c])) for f in filas)) for c in cols}
    linea = " | ".join(c.ljust(anchos[c]) for c in cols)
    print("\n" + linea); print("-" * len(linea))
    for f in filas:
        print(" | ".join(str(f[c]).ljust(anchos[c]) for c in cols))

    ruta_csv = os.path.join(DIR_FIG, "tabla_ojo_horizontal.csv")
    with open(ruta_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(filas)

    # Presupuesto de jitter (dual-Dirac) en el punto de operación
    tj2 = tj_dual_dirac(1e-2, JITTER.sigma_rj_ui, JITTER.dj_pp_ui)
    tj12 = tj_dual_dirac(1e-12, JITTER.sigma_rj_ui, JITTER.dj_pp_ui)
    print(f"\nJitter total (dual-Dirac): TJ@1e-2 = {tj2:.3f} UI | TJ@1e-12 = {tj12:.3f} UI")
    print("\nEntregables generados:")
    for r in (r1, r2, ruta_csv):
        print(f"  - {r}")


if __name__ == "__main__":
    main()
