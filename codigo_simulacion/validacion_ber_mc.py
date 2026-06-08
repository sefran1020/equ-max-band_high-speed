"""
validacion_ber_mc.py — Validación del BER semi-analítico vs conteo Monte Carlo
(Bloque C1 del plan; responde la observación #5 del revisor).

El estimador semi-analítico (Q-función con ajuste gaussiano por nivel) podría ser
optimista. Aquí se contrasta, en la región donde el conteo directo es factible
(BER ~ 1e-2 a 1e-4), la BER semi-analítica contra la BER por CONTEO DE ERRORES
sobre secuencias largas, usando EXACTAMENTE los mismos umbrales de decisión
(puntos medios entre las medias de nivel) que emplea el estimador.

Si ambas curvas coinciden, el supuesto gaussiano no infla las tasas reportadas.

Canales: FR-4 y Megtron 6, longitud fija 10 in. Escenarios C (FFE+CTLE) y
D (FFE+CTLE+DFE). Co-diseño heredado del RC.

Uso:  python validacion_ber_mc.py
"""

import csv
import os
import numpy as np

from cadena_enlace import (
    generar_prbs, modular_pam4, aplicar_ffe, aplicar_canal_freq,
    aplicar_ctle_freq, aplicar_rx_lpf, alinear_full, evaluar_enlace,
    estadisticas_niveles, _decidir_pam4, NIVELES_PAM4,
    SPS, SIGMA_IN, CTLE_ADC, W_Z_RC, C_POST_RC,
)
from canales import construir_canales

DIR_FIG = os.path.join(os.path.dirname(__file__), "figuras")

N_SYM = 120000          # símbolos por realización (para alcanzar BER ~1e-4)
N_REAL = 4              # realizaciones (patrón + ruido) que se promedian
N_DESCARTE = 200
TASAS_GB = np.array([12, 14, 16, 18, 20, 22, 24])   # GBaud


def _waveform_C(H, sym, Rs, rng, con_dfe):
    """Devuelve la forma de onda en la entrada del decisor (escenario C o D)."""
    ts = 1.0 / (Rs * SPS)
    fc = Rs
    ffe = aplicar_ffe(sym, 0.0, 1.0, C_POST_RC)
    y_chan = aplicar_canal_freq(np.repeat(ffe, SPS), ts, H)
    nin = rng.normal(0, SIGMA_IN, len(y_chan))
    yC = aplicar_rx_lpf(aplicar_ctle_freq(y_chan + nin, ts, W_Z_RC, CTLE_ADC), ts, fc)
    return yC, ts


def ber_montecarlo(yC, sym, con_dfe, n_dfe=2):
    """Cuenta errores de símbolo con los MISMOS umbrales que el semi-analítico."""
    m, tx = alinear_full(yC, sym, SPS)
    m, tx = m[N_DESCARTE:], tx[N_DESCARTE:]
    if len(m) < 100:
        return 1.0, 0, 0

    if con_dfe:
        nt = n_dfe
        N = len(m)
        Xreg = np.column_stack([tx[nt - i: N - i] for i in range(nt + 1)])
        h, *_ = np.linalg.lstsq(Xreg, m[nt:N], rcond=None)
        h0 = h[0] if abs(h[0]) > 1e-12 else 1.0
        post = h[1:]
        hist = [0.0] * nt
        m_eq = np.empty_like(m)
        for kk in range(len(m)):
            fb = sum(post[i] * hist[i] for i in range(nt))
            yk = m[kk] - fb
            m_eq[kk] = yk
            hist = [_decidir_pam4(yk, h0)] + hist[:-1]
        m = m_eq

    mu, s_isi, cnt = estadisticas_niveles(m, tx)
    if np.any(cnt == 0):
        return 1.0, 0, len(m)
    o = np.argsort(mu)
    mu_s = mu[o]
    thr = (mu_s[:-1] + mu_s[1:]) / 2.0          # mismos umbrales que el estimador
    niveles_ord = NIVELES_PAM4[o]               # nivel físico asociado a cada bin

    # Índice de bin de cada muestra y del símbolo verdadero
    dec_bin = np.digitize(m, thr)               # 0..3
    # mapa nivel físico -> índice ordenado
    idx_true = np.searchsorted(np.sort(mu_s), mu_s)  # identidad; usamos tx directo
    # bin verdadero: posición de tx en niveles ordenados
    tx_bin = np.zeros(len(tx), dtype=int)
    for b, lv in enumerate(niveles_ord):
        tx_bin[np.isclose(tx, lv)] = b
    errores = int(np.sum(dec_bin != tx_bin))
    ser = errores / len(m)
    ber = ser / 2.0                              # misma definición que el estimador
    return ber, errores, len(m)


def main():
    os.makedirs(DIR_FIG, exist_ok=True)
    canales, _ = construir_canales("fixed")
    objetivo = {"FR-4": canales["FR-4"], "Megtron 6": canales["Megtron 6"]}

    print("== C1 — Validación BER semi-analítico vs Monte Carlo ==")
    print(f"N_sym={N_SYM}, N_real={N_REAL}, conteo con umbrales del estimador\n")

    filas = []
    resultados = {}     # (canal, esc) -> (rates, semi[], mc[])
    for cname, info in objetivo.items():
        H = info["H"]
        for esc, con_dfe in (("C (FFE+CTLE)", False), ("D (+DFE)", True)):
            semis, mcs, errs, tots = [], [], [], []
            for gb in TASAS_GB:
                Rs = gb * 1e9
                bs, bm, ne, nt = [], [], 0, 0
                for r in range(N_REAL):
                    bits = generar_prbs(2 * N_SYM + 600, orden=23,
                                        semilla=[1] * 23 if r == 0 else None)
                    # ventana distinta por realización para variar el patrón
                    off = r * 777 % 500
                    sym = modular_pam4(bits[2 * off:])[:N_SYM]
                    rng = np.random.default_rng(5000 + r)
                    yC, _ = _waveform_C(H, sym, Rs, rng, con_dfe)
                    rsemi = evaluar_enlace(yC, sym, SPS,
                                           dfe_ntaps=2 if con_dfe else 0,
                                           n_descartar=N_DESCARTE)
                    bmc, e, t = ber_montecarlo(yC, sym, con_dfe)
                    bs.append(rsemi["ber"]); bm.append(bmc); ne += e; nt += t
                semi = float(np.mean(bs)); mc = float(np.mean(bm))
                semis.append(semi); mcs.append(mc); errs.append(ne); tots.append(nt)
                print(f"  {cname:<10} {esc:<13} {gb:>4d} GBaud | "
                      f"semi={semi:.2e}  MC={mc:.2e}  ({ne} errs / {nt} sym)")
                filas.append({"canal": cname, "escenario": esc, "GBaud": int(gb),
                              "BER_semi": f"{semi:.3e}", "BER_MC": f"{mc:.3e}",
                              "errores": ne, "simbolos": nt})
            resultados[(cname, esc)] = (TASAS_GB, semis, mcs)

    # CSV
    ruta_csv = os.path.join(DIR_FIG, "tabla_validacion_ber_mc.csv")
    with open(ruta_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["canal", "escenario", "GBaud",
                                           "BER_semi", "BER_MC", "errores", "simbolos"])
        w.writeheader(); w.writerows(filas)

    # Figura overlay
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
        colores = {"FR-4": "#dd6b20", "Megtron 6": "#319795"}
        for ax, esc in zip(axes, ("C (FFE+CTLE)", "D (+DFE)")):
            for cname in objetivo:
                gb, semi, mc = resultados[(cname, esc)]
                col = colores[cname]
                ax.semilogy(gb, np.clip(semi, 1e-12, 1), "-", color=col,
                            label=f"{cname} semi-analytical")
                ax.semilogy(gb, np.clip(mc, 1e-12, 1), "o", color=col, mfc="white",
                            label=f"{cname} Monte Carlo")
            ax.axhline(1e-2, color="gray", ls="--", alpha=0.7)
            ax.set_title(f"Scenario {esc}")
            ax.set_xlabel("Rate (GBaud)"); ax.grid(True, which="both", alpha=0.25)
        axes[0].set_ylabel("BER"); axes[0].legend(fontsize=7, loc="lower right")
        fig.tight_layout()
        ruta_fig = os.path.join(DIR_FIG, "fig_validacion_ber_mc.png")
        fig.savefig(ruta_fig, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"\nFigura: {ruta_fig}")
    except Exception as e:
        print(f"(figura omitida: {e})")
    print(f"CSV: {ruta_csv}")


if __name__ == "__main__":
    main()
