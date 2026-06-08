"""
exportar_touchstone.py — Cross-validation numérico (Bloque C5 del plan).

Exporta los parámetros S del modelo de línea dispersiva a archivos Touchstone
(.s2p) referidos a 50 ohm, ejecuta chequeos de PASIVIDAD y de coherencia con la
función de transferencia usada en la cadena de enlace, y deja los .s2p listos
para una futura validación en SPICE/ADS (sustituto, en Fase 1, de la
"validación cruzada" pedida por la Rec. 01.3 del equipo senior, sin necesidad de
una herramienta comercial).

Parámetros S exactos de una línea uniforme de longitud L, impedancia
característica Zc(f) y constante de propagación gamma(f), referidos a Z0=50 ohm:

    Gamma = (Zc - Z0) / (Zc + Z0)
    e2    = exp(-2*gamma*L)
    den   = 1 - Gamma^2 * e2
    S11 = S22 = Gamma * (1 - e2) / den
    S21 = S12 = (1 - Gamma^2) * exp(-gamma*L) / den

Para una red pasiva con pérdidas se cumple |S11|^2 + |S21|^2 <= 1.

Uso:  python exportar_touchstone.py
"""

import os
import numpy as np

from parametros import SUBSTRATOS, INCH, OPERACION
from modelo_rlgc import LineaTransmision
from canales import FIXED_LEN_INCH

DIR_FIG = os.path.join(os.path.dirname(__file__), "figuras")
Z0 = 50.0
N_F = 1001
F_MAX = 50e9            # 50 GHz


def s_parametros(linea: LineaTransmision, L_m, f):
    """Parámetros S (S11,S21,S12,S22) de la línea de longitud L_m, ref. Z0."""
    f = np.asarray(f, dtype=float)
    gamma = linea.gamma(f)
    zc = linea.zc(f)
    Gamma = (zc - Z0) / (zc + Z0)
    egl = np.exp(-gamma * L_m)
    e2 = egl * egl
    den = 1.0 - Gamma**2 * e2
    s11 = Gamma * (1.0 - e2) / den
    s21 = (1.0 - Gamma**2) * egl / den
    return s11, s21, s21.copy(), s11.copy()


def escribir_s2p(ruta, f, s11, s21, s12, s22):
    """Touchstone v1, formato Real-Imag, referencia 50 ohm."""
    with open(ruta, "w", encoding="utf-8") as fh:
        fh.write("! Touchstone v1.0 generado por exportar_touchstone.py\n")
        fh.write("! Linea de transmision dispersiva (skin + dielectrico)\n")
        fh.write("# Hz S RI R 50\n")
        for i in range(len(f)):
            fh.write(f"{f[i]:.6e} "
                     f"{s11[i].real:.6e} {s11[i].imag:.6e} "
                     f"{s21[i].real:.6e} {s21[i].imag:.6e} "
                     f"{s12[i].real:.6e} {s12[i].imag:.6e} "
                     f"{s22[i].real:.6e} {s22[i].imag:.6e}\n")


def main():
    os.makedirs(DIR_FIG, exist_ok=True)
    f = np.linspace(1e7, F_MAX, N_F)        # 10 MHz .. 50 GHz
    L_m = FIXED_LEN_INCH * INCH             # 10 in
    f_nyq = OPERACION.f_nyquist

    print("== C5 — Export Touchstone + cross-check (linea de 10 in, ref. 50 ohm) ==")
    filas = []
    for sub in SUBSTRATOS:
        ln = LineaTransmision(sub)
        s11, s21, s12, s22 = s_parametros(ln, L_m, f)

        # --- Chequeo de pasividad: |S11|^2 + |S21|^2 <= 1 (+ tolerancia num.) ---
        pot = np.abs(s11) ** 2 + np.abs(s21) ** 2
        pasiva = bool(np.all(pot <= 1.0 + 1e-9))
        pot_max = float(np.max(pot))

        # --- Coherencia con H de la cadena (matched, H=exp(-gamma L)) ---------
        H_chain = ln.H(f, L_m)                       # exp(-gamma L)
        # IL en Nyquist por ambas vías
        il_s21 = float(-20 * np.log10(np.abs(np.interp(f_nyq, f, np.abs(s21)))))
        il_H = float(-20 * np.log10(np.abs(ln.H(np.array([f_nyq]), L_m))[0]))
        # discrepancia media en banda util (<= 2*Nyquist) por reflexiones
        banda = f <= 4e9
        dif_db = 20 * np.log10(np.abs(s21[banda]) + 1e-300) - \
                 20 * np.log10(np.abs(H_chain[banda]) + 1e-300)
        dif_med = float(np.mean(np.abs(dif_db)))

        nombre_f = sub.nombre.replace(" ", "").replace("-", "")
        ruta = os.path.join(DIR_FIG, f"{nombre_f}_10in.s2p")
        escribir_s2p(ruta, f, s11, s21, s12, s22)

        print(f"\n[{sub.nombre}]  ->  {os.path.basename(ruta)}")
        print(f"  Pasividad |S11|^2+|S21|^2 <= 1 : {'OK' if pasiva else 'FALLA'} "
              f"(max={pot_max:.4f})")
        print(f"  IL@Nyquist  S21={il_s21:.2f} dB  vs  H_cadena={il_H:.2f} dB")
        print(f"  |dif S21 vs H| media (<=4 GHz) : {dif_med:.3f} dB "
              f"(reflexiones por Zc!=50)")
        filas.append((sub.nombre, pasiva, pot_max, il_s21, il_H, dif_med))

    # --- Figura de verificación ------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
        for sub in SUBSTRATOS:
            ln = LineaTransmision(sub)
            _, s21, _, _ = s_parametros(ln, L_m, f)
            H_chain = ln.H(f, L_m)
            ax1.plot(f / 1e9, 20 * np.log10(np.abs(s21) + 1e-300), label=f"{sub.nombre} S21")
            ax1.plot(f / 1e9, 20 * np.log10(np.abs(H_chain) + 1e-300), "--",
                     label=f"{sub.nombre} H cadena")
            s11, s21b, _, _ = s_parametros(ln, L_m, f)
            ax2.plot(f / 1e9, np.abs(s11) ** 2 + np.abs(s21b) ** 2,
                     label=f"{sub.nombre}")
        ax1.axvline(f_nyq / 1e9, color="gray", ls=":", alpha=0.7)
        ax1.set_xlabel("Frecuencia (GHz)"); ax1.set_ylabel("Magnitud (dB)")
        ax1.set_title("S21 (Touchstone) vs H de la cadena"); ax1.legend(fontsize=7)
        ax1.grid(alpha=0.3); ax1.set_xlim(0, F_MAX / 1e9)
        ax2.axhline(1.0, color="red", ls="--", alpha=0.7, label="límite pasividad")
        ax2.set_xlabel("Frecuencia (GHz)")
        ax2.set_ylabel(r"$|S_{11}|^2+|S_{21}|^2$")
        ax2.set_title("Chequeo de pasividad"); ax2.legend(fontsize=8)
        ax2.grid(alpha=0.3); ax2.set_xlim(0, F_MAX / 1e9); ax2.set_ylim(0, 1.05)
        fig.tight_layout()
        ruta_fig = os.path.join(DIR_FIG, "fig_touchstone_check.png")
        fig.savefig(ruta_fig, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"\nFigura: {ruta_fig}")
    except Exception as e:
        print(f"(figura omitida: {e})")

    print("\nResumen:")
    print(f"  {'Sustrato':<12} {'pasiva':<7} {'pot_max':<9} {'IL_S21':<8} {'IL_H':<8} {'dif_dB':<7}")
    for n, p, pm, a, b, d in filas:
        print(f"  {n:<12} {str(p):<7} {pm:<9.4f} {a:<8.2f} {b:<8.2f} {d:<7.3f}")


if __name__ == "__main__":
    main()
