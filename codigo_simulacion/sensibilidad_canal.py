"""
sensibilidad_canal.py — Análisis de sensibilidad del canal físico (Bloque C2).

Responde la observación #3 del revisor (validación del canal) por la vía de bajo
costo que él mismo sugiere: un ANÁLISIS DE SENSIBILIDAD ante tangente de pérdida
(tan d), longitud y rugosidad superficial. Convierte el modelo dispersivo de un
"modelo sintético" en un ESTUDIO PARAMÉTRICO con barras de variación.

- tan d   : barre la pérdida dieléctrica (de baja pérdida a FR-4 estándar).
- longitud: barre la longitud física de la traza.
- rugosidad: factor de Hammerstad sobre la resistencia de efecto pelicular,
             K_sr(f) = 1 + (2/pi)*arctan(1.4*(Delta/delta_s)^2), con Delta la
             rugosidad RMS y delta_s la profundidad de penetración. Es una
             corrección de primer orden (se declara como tal en el manuscrito).

Métricas por punto: atenuación en Nyquist (dB/in) y tasa de símbolo viable
(BER<=1e-2, escenario FFE+CTLE, co-diseño heredado del RC).

Uso:  python sensibilidad_canal.py
"""

import csv
import os
import numpy as np

from parametros import Substrato, MU0, SIGMA_CU, OPERACION, INCH, NP_TO_DB
from modelo_rlgc import LineaTransmision
from cadena_enlace import (
    generar_prbs, modular_pam4, evaluar_config, viable_interp,
    W_Z_RC, C_POST_RC,
)

DIR_FIG = os.path.join(os.path.dirname(__file__), "figuras")

N_REAL = 8
N_SYM = 800
N_DFE = 2
TASAS = np.linspace(1e9, 40e9, 20)
GB = TASAS / 1e9
F_NYQ = OPERACION.f_nyquist
FR4 = Substrato("FR-4", 4.3, 0.020)
MEG = Substrato("Megtron 6", 3.6, 0.004)


class LineaRugosa(LineaTransmision):
    """Línea con corrección de rugosidad de Hammerstad sobre R_skin(f)."""
    def __init__(self, substrato, rugosidad_rms_m=0.0):
        super().__init__(substrato)
        self.delta_rms = rugosidad_rms_m

    def _k_sr(self, f):
        f = np.asarray(f, dtype=float)
        if self.delta_rms <= 0:
            return np.ones_like(f)
        delta_s = np.sqrt(1.0 / (np.pi * np.maximum(f, 1.0) * MU0 * SIGMA_CU))
        return 1.0 + (2.0 / np.pi) * np.arctan(1.4 * (self.delta_rms / delta_s) ** 2)

    def R(self, f):
        f = np.asarray(f, dtype=float)
        r_skin = self.k_skin * np.sqrt(np.abs(f)) * self._k_sr(f)
        return np.sqrt(self.R_dc ** 2 + r_skin ** 2)


def tasa_viable_C(H):
    """Tasa viable (GBaud) del escenario C (FFE+CTLE) con EQ heredada del RC."""
    bits = generar_prbs(N_REAL * 2 * N_SYM + 4000, orden=23)
    vs = []
    for r in range(N_REAL):
        bsl = bits[r * 2 * N_SYM: r * 2 * N_SYM + 2 * N_SYM + 200]
        sym = modular_pam4(bsl)[:N_SYM]
        rng = np.random.default_rng(1000 + r)
        bl = [evaluar_config(H, sym, Rs, W_Z_RC, C_POST_RC, N_DFE, rng)["C"]["ber"]
              for Rs in TASAS]
        vs.append(viable_interp(bl, GB))
    return float(np.mean(vs)), float(np.std(vs))


def aten_nyq(linea, L_m):
    alpha = float(np.real(linea.gamma(np.array([F_NYQ])))[0])
    return alpha * NP_TO_DB * INCH       # dB/in


def main():
    os.makedirs(DIR_FIG, exist_ok=True)
    print("== C2 — Sensibilidad del canal (tan d, longitud, rugosidad) ==\n")
    L10 = 10.0 * INCH
    filas = []

    # ---- Barrido 1: tangente de pérdida (longitud fija 10 in) ----------------
    print("[1] tan d  (10 in)")
    tand_vals = np.array([0.002, 0.004, 0.008, 0.012, 0.016, 0.020, 0.025])
    s1 = []
    for td in tand_vals:
        ln = LineaTransmision(Substrato("custom", 4.0, float(td)))
        a = aten_nyq(ln, L10)
        v, sd = tasa_viable_C(lambda f, ln=ln: ln.H(f, L10))
        s1.append((td, a, v, sd))
        filas.append({"barrido": "tan_d", "valor": td, "aten_dbin_nyq": round(a, 3),
                      "tasa_viable_GBaud": round(v, 2), "desv": round(sd, 2)})
        print(f"  tan d={td:.3f}  aten={a:.3f} dB/in  viable={v:.1f}±{sd:.1f} GBaud")

    # ---- Barrido 2: longitud (FR-4) ------------------------------------------
    print("\n[2] longitud (FR-4)")
    long_vals = np.array([5, 8, 10, 14, 20, 26, 32])
    s2 = []
    for Lin in long_vals:
        Lm = float(Lin) * INCH
        ln = LineaTransmision(FR4)
        a = aten_nyq(ln, Lm)
        v, sd = tasa_viable_C(lambda f, ln=ln, Lm=Lm: ln.H(f, Lm))
        s2.append((Lin, a, v, sd))
        filas.append({"barrido": "longitud_in", "valor": int(Lin), "aten_dbin_nyq": round(a, 3),
                      "tasa_viable_GBaud": round(v, 2), "desv": round(sd, 2)})
        print(f"  L={Lin:>2d} in  aten={a:.3f} dB/in  viable={v:.1f}±{sd:.1f} GBaud")

    # ---- Barrido 3: rugosidad (FR-4, 10 in) ----------------------------------
    print("\n[3] rugosidad RMS (FR-4, 10 in)")
    rug_vals_um = np.array([0.0, 0.1, 0.3, 0.5, 0.7, 1.0])
    s3 = []
    for ru in rug_vals_um:
        ln = LineaRugosa(FR4, rugosidad_rms_m=float(ru) * 1e-6)
        a = aten_nyq(ln, L10)
        v, sd = tasa_viable_C(lambda f, ln=ln: ln.H(f, L10))
        s3.append((ru, a, v, sd))
        filas.append({"barrido": "rugosidad_um", "valor": ru, "aten_dbin_nyq": round(a, 3),
                      "tasa_viable_GBaud": round(v, 2), "desv": round(sd, 2)})
        print(f"  Delta={ru:.1f} um  aten={a:.3f} dB/in  viable={v:.1f}±{sd:.1f} GBaud")

    # ---- CSV -----------------------------------------------------------------
    ruta_csv = os.path.join(DIR_FIG, "tabla_sensibilidad_canal.csv")
    with open(ruta_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["barrido", "valor", "aten_dbin_nyq",
                                           "tasa_viable_GBaud", "desv"])
        w.writeheader(); w.writerows(filas)

    # ---- Figura --------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(1, 3, figsize=(13.5, 4.2))

        def panel(ax, datos, xlabel, titulo, color):
            x = [d[0] for d in datos]; v = [d[2] for d in datos]; sd = [d[3] for d in datos]
            ax.errorbar(x, v, yerr=sd, marker="o", color=color, capsize=3)
            ax.set_xlabel(xlabel); ax.set_ylabel("Viable rate (GBaud)")
            ax.grid(alpha=0.3)

        panel(axs[0], s1, r"$\tan\delta$", "", "#dd6b20")
        panel(axs[1], s2, "Length (in)", "", "#3182ce")
        panel(axs[2], s3, r"RMS roughness ($\mu$m)", "", "#805ad5")
        fig.tight_layout()
        ruta_fig = os.path.join(DIR_FIG, "fig_sensibilidad_canal.png")
        fig.savefig(ruta_fig, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"\nFigura: {ruta_fig}")
    except Exception as e:
        print(f"(figura omitida: {e})")
    print(f"CSV: {ruta_csv}")


if __name__ == "__main__":
    main()
