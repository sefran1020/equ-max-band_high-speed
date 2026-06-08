"""
generar_comparativa.py — Programa principal de la Recomendación 01.

Genera los dos ENTREGABLES pedidos por el equipo senior:

  (1) Atenuación en dB/pulgada vs frecuencia  -> fig_atenuacion_dbporpulgada.png
  (2) Comparativa de respuestas al pulso       -> fig_respuesta_pulso.png
      (canal RC difusivo del artículo  vs  línea de transmisión dispersiva
       FR-4 y Megtron 6, igualadas a la MISMA pérdida en Nyquist)

Además imprime y guarda una tabla resumen (tabla_comparativa.csv) lista para
citar en una futura versión del manuscrito (cierra la limitación declarada en
discusion.tex: "modelo RC (sin inductancia/skin-effect)").

Uso:
    python generar_comparativa.py

Las figuras y la tabla se guardan en ./figuras/.
"""

import csv
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")  # backend sin ventana, para ejecución local headless
import matplotlib.pyplot as plt

from parametros import SUBSTRATOS, OPERACION, REJILLA, INCH
from modelo_rlgc import LineaTransmision
from modelo_rc import CanalRC
from respuesta_pulso import respuesta_al_pulso, metricas_pulso

DIR_FIG = os.path.join(os.path.dirname(__file__), "figuras")
DPI = 300


def asegurar_dir():
    os.makedirs(DIR_FIG, exist_ok=True)


def il_db(mag_complejo):
    """Pérdida de inserción positiva en dB a partir de |H|."""
    return -20.0 * np.log10(np.abs(mag_complejo))


# --------------------------------------------------------------------------- #
# Entregable 1: atenuación en dB/pulgada
# --------------------------------------------------------------------------- #
def figura_atenuacion(lineas):
    f = np.logspace(np.log10(REJILLA.f_min_bode),
                    np.log10(REJILLA.f_max_bode), REJILLA.n_bode)

    fig, ax = plt.subplots(figsize=(7.0, 4.3))
    for ln in lineas:
        ax.plot(f / 1e9, ln.atenuacion_db_por_pulgada(f),
                lw=2, label=ln.substrato.nombre)

    ax.axvline(OPERACION.f_nyquist / 1e9, color="0.4", ls="--", lw=1)
    ax.text(OPERACION.f_nyquist / 1e9 * 1.05, ax.get_ylim()[1] * 0.05,
            "Nyquist 2 GHz", color="0.3", fontsize=8)

    ax.set_xscale("log")
    ax.set_xlabel("Frequency [GHz]")
    ax.set_ylabel("Attenuation  [dB/inch]")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    ruta = os.path.join(DIR_FIG, "fig_atenuacion_dbporpulgada.png")
    fig.savefig(ruta, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return ruta


# --------------------------------------------------------------------------- #
# Entregable 2: comparativa de respuestas al pulso
# --------------------------------------------------------------------------- #
def figura_pulsos(canal_rc, lineas):
    # Pérdida de inserción del RC en Nyquist -> objetivo para igualar longitudes
    il_obj = float(il_db(canal_rc.H(np.array([OPERACION.f_nyquist])))[0])

    # Respuesta del canal RC
    t, x, y_rc = respuesta_al_pulso(canal_rc)

    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    # ventana de tiempo recortada alrededor del pulso para legibilidad
    t0_idx = int(0.10 * len(t))
    t1_idx = int(0.45 * len(t))
    sl = slice(t0_idx, t1_idx)
    t_ns = t[sl] * 1e9

    ax.plot(t_ns, x[sl], color="0.6", lw=1.2, ls=":", label="Input (1 UI)")
    ax.plot(t_ns, y_rc[sl], color="k", lw=2,
            label=f"RC diffusive (tau={OPERACION.tau_rc*1e12:.0f} ps)")

    longitudes = {}
    for ln in lineas:
        long_m = ln.longitud_para_il(il_obj, OPERACION.f_nyquist)
        longitudes[ln.substrato.nombre] = long_m
        _, _, y = respuesta_al_pulso(ln, longitud=long_m)
        ax.plot(t_ns, y[sl], lw=2,
                label=f"{ln.substrato.nombre} ({long_m/INCH:.1f} in, IL={il_obj:.1f} dB @Nyq)")

    ax.set_xlabel("Time [ns]")
    ax.set_ylabel("Amplitude [V]")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    ruta = os.path.join(DIR_FIG, "fig_respuesta_pulso.png")
    fig.savefig(ruta, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return ruta, il_obj, longitudes


# --------------------------------------------------------------------------- #
# Tabla resumen
# --------------------------------------------------------------------------- #
def tabla_resumen(canal_rc, lineas, il_obj, longitudes):
    filas = []

    # Fila RC
    t, x, y_rc = respuesta_al_pulso(canal_rc)
    m = metricas_pulso(t, y_rc)
    filas.append({
        "canal": "RC difusivo (cosh)",
        "substrato": "-",
        "long_pulg": "-",
        "f_3dB_GHz": round(canal_rc.f_3db() / 1e9, 3),
        "atten_Nyq_dB_in": "-",
        "IL_Nyq_dB": round(il_obj, 2),
        "pico_pulso_V": round(m["pico"], 4),
        "fwhm_ps": round(m["fwhm_s"] * 1e12, 1),
    })

    f_nyq = np.array([OPERACION.f_nyquist])
    for ln in lineas:
        long_m = longitudes[ln.substrato.nombre]
        _, _, y = respuesta_al_pulso(ln, longitud=long_m)
        m = metricas_pulso(t, y)
        atten_in = float(ln.atenuacion_db_por_pulgada(f_nyq)[0])
        filas.append({
            "canal": "Línea TX dispersiva",
            "substrato": ln.substrato.nombre,
            "long_pulg": round(long_m / INCH, 2),
            "f_3dB_GHz": "-",
            "atten_Nyq_dB_in": round(atten_in, 3),
            "IL_Nyq_dB": round(il_obj, 2),
            "pico_pulso_V": round(m["pico"], 4),
            "fwhm_ps": round(m["fwhm_s"] * 1e12, 1),
        })

    # Imprimir
    cols = ["canal", "substrato", "long_pulg", "f_3dB_GHz",
            "atten_Nyq_dB_in", "IL_Nyq_dB", "pico_pulso_V", "fwhm_ps"]
    anchos = {c: max(len(c), *(len(str(fila[c])) for fila in filas)) for c in cols}
    linea = " | ".join(c.ljust(anchos[c]) for c in cols)
    print("\n" + linea)
    print("-" * len(linea))
    for fila in filas:
        print(" | ".join(str(fila[c]).ljust(anchos[c]) for c in cols))

    # Guardar CSV
    ruta = os.path.join(DIR_FIG, "tabla_comparativa.csv")
    with open(ruta, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(filas)
    return ruta


# --------------------------------------------------------------------------- #
def main():
    asegurar_dir()
    lineas = [LineaTransmision(s) for s in SUBSTRATOS]
    canal_rc = CanalRC()

    print("== Recomendación 01: Evolución del modelo de canal (física realista) ==")
    print(f"Punto de operación: {OPERACION.baud_op/1e9:.0f} GBaud, "
          f"Nyquist {OPERACION.f_nyquist/1e9:.0f} GHz, "
          f"tau_RC = {OPERACION.tau_rc*1e12:.0f} ps")

    r1 = figura_atenuacion(lineas)
    r2, il_obj, longitudes = figura_pulsos(canal_rc, lineas)
    r3 = tabla_resumen(canal_rc, lineas, il_obj, longitudes)

    print("\nEntregables generados:")
    print(f"  - {r1}")
    print(f"  - {r2}")
    print(f"  - {r3}")


if __name__ == "__main__":
    main()
