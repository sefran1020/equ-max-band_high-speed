"""
convergencia_nsym.py — Estudio de convergencia en el número de símbolos por
realización (responde la observación #2 del revisor: "Nsym=1000 es muy pequeño").

Comprueba si las medias/desviaciones por nivel del ojo, el BER semi-analítico y
la TASA VIABLE convergen al aumentar Nsym. Si los valores con Nsym=1000 ya están
dentro de la dispersión de los de Nsym grande, queda JUSTIFICADO; en caso
contrario, indica el Nsym mínimo necesario.

Casos representativos (los más exigentes: ISI larga, canal dispersivo):
  - FR-4 fixed (10 in), escenarios C (FFE+CTLE) y D (+DFE)
  - Megtron 6 fixed (10 in), escenario C

Uso:  python convergencia_nsym.py
"""

import csv
import os
import numpy as np

from cadena_enlace import (
    generar_prbs, modular_pam4, evaluar_config, viable_interp,
    aplicar_ffe, aplicar_canal_freq, aplicar_ctle_freq, aplicar_rx_lpf,
    alinear_full, estadisticas_niveles,
    SPS, SIGMA_IN, CTLE_ADC, W_Z_RC, C_POST_RC,
)
from canales import construir_canales

DIR_FIG = os.path.join(os.path.dirname(__file__), "figuras")

NSYM_LIST = [500, 1000, 2000, 4000, 8000, 16000]
N_REAL = 12
N_DFE = 2
TASAS = np.linspace(1e9, 40e9, 20)
GB = TASAS / 1e9
RS_STAT = 16e9          # tasa para el chequeo de estadísticas por nivel (cerca del umbral)


def stats_por_nivel(H, sym, Rs):
    """Devuelve (mu, sigma) por nivel del escenario C a la tasa Rs."""
    ts = 1.0 / (Rs * SPS); fc = Rs
    ffe = aplicar_ffe(sym, 0.0, 1.0, C_POST_RC)
    y_chan = aplicar_canal_freq(np.repeat(ffe, SPS), ts, H)
    rng = np.random.default_rng(123)
    nin = rng.normal(0, SIGMA_IN, len(y_chan))
    yC = aplicar_rx_lpf(aplicar_ctle_freq(y_chan + nin, ts, W_Z_RC, CTLE_ADC), ts, fc)
    m, tx = alinear_full(yC, sym, SPS)
    m, tx = m[100:], tx[100:]
    mu, s, cnt = estadisticas_niveles(m, tx)
    o = np.argsort(mu)
    return mu[o], s[o]


def viable_para_nsym(H, nsym, key):
    """Tasa viable (media±desv, N_REAL) del escenario `key` (C o D) con `nsym`."""
    bits = generar_prbs(N_REAL * 2 * nsym + 5000, orden=23)
    vs = []
    for r in range(N_REAL):
        bsl = bits[r * 2 * nsym: r * 2 * nsym + 2 * nsym + 200]
        sym = modular_pam4(bsl)[:nsym]
        rng = np.random.default_rng(1000 + r)
        bers = [evaluar_config(H, sym, Rs, W_Z_RC, C_POST_RC, N_DFE, rng)[key]["ber"]
                for Rs in TASAS]
        vs.append(viable_interp(bers, GB))
    return float(np.mean(vs)), float(np.std(vs))


def main():
    os.makedirs(DIR_FIG, exist_ok=True)
    canales, _ = construir_canales("fixed")
    casos = [("FR-4", "C (FFE+CTLE)", "C"),
             ("FR-4", "D (+DFE)", "D"),
             ("Megtron 6", "C (FFE+CTLE)", "C")]

    print("== #2 — Convergencia en Nsym (tasa viable, fixed, 12 realizaciones) ==\n")
    filas = []
    series = {}
    for cname, etiq, key in casos:
        H = canales[cname]["H"]
        print(f"[{cname} — {etiq}]")
        med_list, sd_list = [], []
        for ns in NSYM_LIST:
            med, sd = viable_para_nsym(H, ns, key)
            med_list.append(med); sd_list.append(sd)
            filas.append({"canal": cname, "escenario": etiq, "Nsym": ns,
                          "tasa_viable_GBaud": round(med, 2), "desv": round(sd, 2)})
            print(f"   Nsym={ns:>5d}  viable={med:5.2f} ± {sd:.2f} GBaud")
        series[(cname, etiq)] = (med_list, sd_list)
        ref = med_list[-1]
        d1000 = abs(med_list[1] - ref) / ref * 100
        print(f"   -> |Nsym=1000 - Nsym=16000| = {d1000:.1f}% del valor convergido\n")

    # Convergencia de estadísticas por nivel (FR-4, escenario C, Rs=16 GBaud)
    print("[Estadísticas por nivel] FR-4 escenario C @16 GBaud (mu, sigma):")
    stat_rows = []
    for ns in NSYM_LIST:
        bits = generar_prbs(2 * ns + 600, orden=23)
        sym = modular_pam4(bits)[:ns]
        mu, s = stats_por_nivel(canales["FR-4"]["H"], sym, RS_STAT)
        stat_rows.append((ns, mu, s))
        print(f"   Nsym={ns:>5d}  mu={np.array2string(mu, precision=2)}  "
              f"sigma={np.array2string(s, precision=3)}")
        filas.append({"canal": "FR-4", "escenario": "stats_mu_sigma", "Nsym": ns,
                      "tasa_viable_GBaud": f"mu={np.array2string(mu, precision=3)}",
                      "desv": f"sig={np.array2string(s, precision=3)}"})

    ruta_csv = os.path.join(DIR_FIG, "tabla_convergencia_nsym.csv")
    with open(ruta_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["canal", "escenario", "Nsym",
                                           "tasa_viable_GBaud", "desv"])
        w.writeheader(); w.writerows(filas)

    # Figura
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.4))
        for (cname, etiq), (med, sd) in series.items():
            ax1.errorbar(NSYM_LIST, med, yerr=sd, marker="o", capsize=3,
                         label=f"{cname} {etiq}")
        ax1.axvline(1000, color="red", ls=":", alpha=0.6, label="Nsym=1000 (usado)")
        ax1.set_xscale("log"); ax1.set_xlabel("Nsym (símbolos/realización)")
        ax1.set_ylabel("Tasa viable (GBaud)")
        ax1.set_title("Convergencia de la tasa viable", fontsize=10, fontweight="bold")
        ax1.grid(alpha=0.3, which="both"); ax1.legend(fontsize=7)
        # sigma por nivel vs Nsym (FR-4 C)
        sig_arr = np.array([s for _, _, s in stat_rows])
        for lvl in range(4):
            ax2.plot(NSYM_LIST, sig_arr[:, lvl], marker="s", label=f"nivel {lvl}")
        ax2.axvline(1000, color="red", ls=":", alpha=0.6)
        ax2.set_xscale("log"); ax2.set_xlabel("Nsym")
        ax2.set_ylabel(r"$\sigma$ por nivel (V)")
        ax2.set_title("Convergencia de σ por nivel (FR-4 C @16 GBaud)",
                      fontsize=10, fontweight="bold")
        ax2.grid(alpha=0.3, which="both"); ax2.legend(fontsize=7)
        fig.suptitle("Convergencia del BER semi-analítico con el número de símbolos",
                     fontsize=12, fontweight="bold", color="#1a365d")
        fig.tight_layout()
        ruta_fig = os.path.join(DIR_FIG, "fig_convergencia_nsym.png")
        fig.savefig(ruta_fig, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"\nFigura: {ruta_fig}")
    except Exception as e:
        print(f"(figura omitida: {e})")
    print(f"CSV: {ruta_csv}")


if __name__ == "__main__":
    main()
