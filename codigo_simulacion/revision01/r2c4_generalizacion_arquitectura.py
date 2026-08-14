"""
r2c4_generalizacion_arquitectura.py — Revisor 2, comentario 4.

    "The proposed design rule states that the optimal CTLE zero scales with
     channel bandwidth. Please discuss whether this observation can be
     generalized to different equalizer architectures or whether it is limited
     to the specific CTLE and optimization setup adopted in this work."

Pregunta legítima: la regla se dedujo con UN ecualizador concreto (un cero, dos
polos en 15 y 30 GHz, A_dc = 2.5) y UN objetivo de optimización (tasa viable de
FFE+CTLE). Este script comprueba si sobrevive al cambio de ambos.

Se repite el mismo re-DoE del artículo sobre cinco canales de anchos de banda
muy distintos, para seis configuraciones:

  A0  linea base            : 1 cero, polos 15/30 GHz, A_dc 2.5   (la del artículo)
  A1  polos de banda ancha  : 1 cero, polos 25/50 GHz
  A2  polos de banda angosta: 1 cero, polos  8/16 GHz
  A3  dos ceros             : ceros w_z y 2*w_z, polos 15/30/30 GHz
  A4  ganancia DC alta      : 1 cero, polos 15/30 GHz, A_dc 4.0
  A5  objetivo por apertura : arquitectura A0, pero el DoE maximiza la APERTURA
                              DE OJO en lugar de la tasa viable

A5 es la que responde a la segunda mitad de la pregunta: si la regla dependiera
del montaje de optimización, cambiar el objetivo debería romperla.

Para cada par (arquitectura, canal) se ajusta  w_z* = k * BW^p  y se reporta el
exponente p con su R^2. Si la regla generaliza, p debe ser consistentemente
positivo y de magnitud parecida en las seis configuraciones.

Verificaciones incorporadas
---------------------------
1. La arquitectura A0 evaluada con el mismo rango de búsqueda del artículo debe
   reproducir sus w_z re-optimizados (2.35 / 2.65 / 5.35 / 9.85 GHz). El script
   aborta si no lo hace.
2. El evaluador parametrizado se contrasta contra `evaluar_config()` del
   artículo con los parámetros de A0: deben devolver el MISMO BER.

Nota sobre el rango de búsqueda: el artículo exploró w_z en 1-10 GHz y el
óptimo de Megtron 6 'fixed' salió 9.85 GHz, es decir PEGADO al borde superior.
Para el estudio de escalamiento el rango se amplía a 1-20 GHz, porque un óptimo
censurado por el borde aplanaría artificialmente la pendiente.

Entregables en ./figuras/ y ./tablas/:
  - fig_generalizacion_wz.png       : (a) w_z* vs ancho de banda del canal en
                                      log-log con su ajuste, (b) exponente
                                      ajustado por arquitectura.
  - tabla_r2c4_generalizacion.csv   : w_z* por arquitectura y canal, y el ajuste.

Uso:
    python revision01/r2c4_generalizacion_arquitectura.py
"""

import csv
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from cadena_enlace import (                                     # noqa: E402
    SPS, SIGMA_IN, CTLE_P1, CTLE_P2, CTLE_ADC, CPOST_MIN, CPOST_MAX,
    generar_prbs, modular_pam4, aplicar_ffe, aplicar_canal_freq,
    aplicar_rx_lpf, evaluar_enlace, evaluar_config, viable_interp, lhs_center,
)
from canales import construir_canales                           # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
DIR_FIG = os.path.join(AQUI, "figuras")
DIR_TAB = os.path.join(AQUI, "tablas")
DPI = 300

GHZ = 2 * np.pi * 1e9

# --- Arquitecturas a comparar --------------------------------------------- #
# ceros: multiplicadores de w_z ; polos: en rad/s ; objetivo: 'tasa' o 'ojo'
ARQUITECTURAS = [
    ("A0 baseline", dict(ceros=(1.0,), polos=(CTLE_P1, CTLE_P2), a_dc=CTLE_ADC,
                         objetivo="tasa"), "#2b6cb0"),
    ("A1 wideband poles", dict(ceros=(1.0,), polos=(25 * GHZ, 50 * GHZ),
                               a_dc=CTLE_ADC, objetivo="tasa"), "#319795"),
    ("A2 narrowband poles", dict(ceros=(1.0,), polos=(8 * GHZ, 16 * GHZ),
                                 a_dc=CTLE_ADC, objetivo="tasa"), "#dd6b20"),
    ("A3 two zeros", dict(ceros=(1.0, 2.0), polos=(15 * GHZ, 30 * GHZ, 30 * GHZ),
                          a_dc=CTLE_ADC, objetivo="tasa"), "#805ad5"),
    ("A4 high DC gain", dict(ceros=(1.0,), polos=(CTLE_P1, CTLE_P2), a_dc=4.0,
                             objetivo="tasa"), "#c05621"),
    ("A5 eye objective", dict(ceros=(1.0,), polos=(CTLE_P1, CTLE_P2),
                              a_dc=CTLE_ADC, objetivo="ojo"), "#e53e3e"),
]

# --- Canales: los cuatro del artículo más el RC como ancla de banda baja --- #
CANALES = [("RC N=5", "equal_il"),
           ("FR-4", "equal_il"), ("Megtron 6", "equal_il"),
           ("FR-4", "fixed"), ("Megtron 6", "fixed")]

# w_z re-optimizados publicados (Tabla IV), para la verificación.
WZ_PUBLICADO = {("FR-4", "equal_il"): 2.35, ("Megtron 6", "equal_il"): 2.65,
                ("FR-4", "fixed"): 5.35, ("Megtron 6", "fixed"): 9.85}

N_DFE = 2
N_LHS = 30
N_SYM = 800
# Rejilla de tasa del DoE. La del artículo (2-26 GBaud, 13 puntos) forma parte
# del procedimiento publicado y es la que hay que usar para reproducirlo. Para
# el estudio se amplía, porque en los canales más anchos la tasa viable se sale
# de esa rejilla y el objetivo deja de discriminar entre valores de w_z: el
# óptimo se vuelve un empate y el LHS se queda con el mayor que muestreó.
TASAS_ART = np.linspace(2e9, 26e9, 13)
TASAS_EST = np.linspace(2e9, 56e9, 19)
SATURA_TOL = 1e-6
RS_OJO = 8e9                 # tasa a la que se mide la apertura para A5
DB_BW = -10.0                # nivel de pérdida que define el ancho de banda

WZ_MIN_ART, WZ_MAX_ART = 1e9, 10e9      # rango del artículo (verificación)
WZ_MIN_EST, WZ_MAX_EST = 1e9, 20e9      # rango ampliado (estudio)


# --------------------------------------------------------------------------- #
# CTLE parametrizable y evaluación de una configuración
# --------------------------------------------------------------------------- #
def ctle_arq(f, w_z, ceros, polos, a_dc):
    """H_CTLE(jw) con un número arbitrario de ceros (múltiplos de w_z) y polos."""
    jw = 1j * 2.0 * np.pi * np.asarray(f, dtype=float)
    num = np.ones_like(jw)
    for m in ceros:
        num = num * (1.0 + jw / (m * w_z))
    den = np.ones_like(jw)
    for p in polos:
        den = den * (1.0 + jw / p)
    return a_dc * num / den


def evaluar_arq(H_func, sym, Rs, w_z, c_post, arq, rng):
    """Escenario C (FFE+CTLE) y D (+DFE) con la arquitectura dada.

    Compone las MISMAS primitivas del artículo y termina en `evaluar_enlace()`,
    de modo que el detector y el cálculo del BER no se duplican.
    """
    ts = 1.0 / (Rs * SPS)
    ffe = aplicar_ffe(sym, 0.0, 1.0, c_post)
    y_chan = aplicar_canal_freq(np.repeat(ffe, SPS), ts, H_func)

    nin = rng.normal(0, SIGMA_IN, len(y_chan))       # ruido referido a la entrada
    y = y_chan + nin
    n = len(y)
    f = np.fft.rfftfreq(n, ts)
    y = np.fft.irfft(np.fft.rfft(y) * ctle_arq(f, w_z, arq["ceros"],
                                               arq["polos"], arq["a_dc"]), n=n)
    yC = aplicar_rx_lpf(y, ts, Rs)
    return evaluar_enlace(yC, sym, SPS), evaluar_enlace(yC, sym, SPS, dfe_ntaps=N_DFE)


def objetivo_doe(H_func, sym, w_z, c_post, arq, tasas):
    """Valor que el DoE maximiza: tasa viable de FFE+CTLE, o apertura de ojo."""
    if arq["objetivo"] == "ojo":
        rng = np.random.default_rng(2026)
        rC, _ = evaluar_arq(H_func, sym, RS_OJO, w_z, c_post, arq, rng)
        return rC["eye_min"]
    rng = np.random.default_rng(2026)
    bers = [evaluar_arq(H_func, sym, Rs, w_z, c_post, arq, rng)[0]["ber"]
            for Rs in tasas]
    return viable_interp(bers, tasas / 1e9)


def re_doe(H_func, arq, wz_min, wz_max, tasas):
    """Mismo LHS del artículo (semilla 42, 30 muestras) sobre (w_z, c_post).

    Devuelve además `saturado`: cierto cuando el mejor objetivo coincide con el
    tope de la rejilla de tasa, es decir cuando el enlace ya no cruza el umbral
    dentro del barrido. En ese caso el óptimo NO es identificable —muchos w_z
    empatan y el LHS se queda con el mayor que muestreó— y el punto no debe
    entrar en el ajuste de escalamiento.
    """
    dis = lhs_center(2, N_LHS, np.random.default_rng(42))
    wz = 2 * np.pi * (wz_min + dis[:, 0] * (wz_max - wz_min))
    cp = CPOST_MIN + dis[:, 1] * (CPOST_MAX - CPOST_MIN)
    sym = modular_pam4(generar_prbs(2 * N_SYM + 400, orden=23))[:N_SYM]

    mejor = None
    for i in range(N_LHS):
        v = objetivo_doe(H_func, sym, wz[i], cp[i], arq, tasas)
        if mejor is None or v > mejor["valor"]:
            mejor = {"w_z": float(wz[i]), "c_post": float(cp[i]), "valor": float(v)}
    mejor["saturado"] = (arq["objetivo"] == "tasa"
                         and mejor["valor"] >= tasas[-1] / 1e9 - SATURA_TOL)
    return mejor


# --------------------------------------------------------------------------- #
# Ancho de banda del canal
# --------------------------------------------------------------------------- #
def frecuencia_a_perdida(H_func, db, f_max=60e9, n=20000):
    """Frecuencia [Hz] a la que |H| ha caído `db` dB respecto de la banda baja."""
    f = np.linspace(1e6, f_max, n)
    mag = np.abs(H_func(f))
    ref = mag[0]
    perdida = 20 * np.log10(mag / ref)
    bajo = np.where(perdida <= db)[0]
    if len(bajo) == 0:
        return float("nan")
    i = bajo[0]
    if i == 0:
        return float(f[0])
    y0, y1 = perdida[i - 1], perdida[i]
    frac = (db - y0) / (y1 - y0)
    return float(f[i - 1] + frac * (f[i] - f[i - 1]))


def ajustar_ley_potencia(bw_ghz, wz_ghz):
    """w_z = k * BW^p  ->  regresión lineal en log-log. Devuelve (p, k, R^2)."""
    x = np.log10(np.asarray(bw_ghz, float))
    y = np.log10(np.asarray(wz_ghz, float))
    p, logk = np.polyfit(x, y, 1)
    pred = p * x + logk
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(p), float(10 ** logk), r2


# --------------------------------------------------------------------------- #
def verificar_evaluador(canales_por_criterio):
    """El evaluador parametrizado con A0 debe dar el MISMO BER que el del artículo.

    `evaluar_config()` extrae ruido para los escenarios A y B antes de llegar al
    C, así que aquí se consumen esas dos extracciones para que el generador
    quede alineado y la comparación sea bit a bit, no aproximada.
    """
    H = canales_por_criterio["fixed"]["FR-4"]["H"]
    sym = modular_pam4(generar_prbs(2 * 400 + 400, orden=23))[:400]
    arq = ARQUITECTURAS[0][1]
    w_z, c_post, Rs = 2 * np.pi * 5.35e9, -0.15, 12e9

    rng = np.random.default_rng(7)
    n = len(sym) * SPS
    rng.normal(0, SIGMA_IN, n)          # escenario A del evaluador original
    rng.normal(0, SIGMA_IN, n)          # escenario B
    rC, rD = evaluar_arq(H, sym, Rs, w_z, c_post, arq, rng)

    ref = evaluar_config(H, sym, Rs, w_z, c_post, N_DFE, np.random.default_rng(7))
    return (rC["ber"], ref["C"]["ber"]), (rD["ber"], ref["D"]["ber"])


def main():
    os.makedirs(DIR_FIG, exist_ok=True)
    os.makedirs(DIR_TAB, exist_ok=True)

    canales_por_criterio = {c: construir_canales(c)[0] for c in ("equal_il", "fixed")}

    print("== R2.4 — ¿la regla w_z ~ ancho de banda depende de la arquitectura? ==")

    # --- Descriptor de ancho de banda -------------------------------------- #
    # El artículo hablaba de "ancho de banda del canal" sin precisar a qué nivel
    # de pérdida. Importa: el corte clásico de -3 dB NO ordena los canales igual
    # que su w_z óptimo (FR-4 'equal_il' tiene MÁS banda a -3 dB que Megtron
    # 'equal_il' y sin embargo un w_z MENOR). El descriptor que sí los ordena es
    # la frecuencia a -10 dB, es decir donde el canal ya ha perdido lo bastante
    # como para que la ecualización tenga trabajo que hacer.
    bw, bw3 = {}, {}
    print("\nDescriptores de banda por canal (GHz):")
    print("   %-12s %-9s %10s %10s" % ("canal", "criterio", "f(-3dB)", "f(-10dB)"))
    for nombre, criterio in CANALES:
        H = canales_por_criterio[criterio][nombre]["H"]
        bw3[(nombre, criterio)] = frecuencia_a_perdida(H, -3.0) / 1e9
        bw[(nombre, criterio)] = frecuencia_a_perdida(H, DB_BW) / 1e9
        print("   %-12s %-9s %10.2f %10.2f"
              % (nombre, criterio, bw3[(nombre, criterio)], bw[(nombre, criterio)]))

    # --- Verificación 1: el evaluador parametrizado == el del artículo ----- #
    (bc, bc_ref), (bd, bd_ref) = verificar_evaluador(canales_por_criterio)
    print("\nVerificación: evaluador parametrizado (A0) contra evaluar_config()")
    print("   FFE+CTLE     : %.6e vs %.6e" % (bc, bc_ref))
    print("   FFE+CTLE+DFE : %.6e vs %.6e" % (bd, bd_ref))
    if not (np.isclose(bc, bc_ref, rtol=1e-12) and np.isclose(bd, bd_ref, rtol=1e-12)):
        raise SystemExit("ERROR: el evaluador parametrizado no reproduce el del artículo.")

    # --- Verificación 2: A0 con el rango del artículo --------------------- #
    print("\nVerificación: A0 con el rango de búsqueda del artículo (1-10 GHz)")
    ok = True
    for (nombre, criterio), esperado in WZ_PUBLICADO.items():
        H = canales_por_criterio[criterio][nombre]["H"]
        m = re_doe(H, ARQUITECTURAS[0][1], WZ_MIN_ART, WZ_MAX_ART, TASAS_ART)
        obtenido = m["w_z"] / (2 * np.pi * 1e9)
        marca = "ok" if abs(obtenido - esperado) < 0.01 else "DIFIERE"
        if marca != "ok":
            ok = False
        print("   %-12s %-9s  w_z = %5.2f GHz  (publicado %5.2f)  %s"
              % (nombre, criterio, obtenido, esperado, marca))
    if not ok:
        raise SystemExit("ERROR: no se reproduce el re-DoE publicado.")

    # --- Estudio: seis arquitecturas, rango ampliado ---------------------- #
    print("\nRe-DoE por arquitectura (rango ampliado 1-20 GHz), w_z* en GHz:")
    encabezado = "   %-22s" % "arquitectura" + "".join(
        "%12s" % ("%s/%s" % (n.split()[0], c[:3])) for n, c in CANALES)
    print(encabezado)

    resultados, ajustes, saturados = {}, {}, {}
    for etiqueta, arq, _ in ARQUITECTURAS:
        wz_ghz, sat = [], []
        for nombre, criterio in CANALES:
            H = canales_por_criterio[criterio][nombre]["H"]
            m = re_doe(H, arq, WZ_MIN_EST, WZ_MAX_EST, TASAS_EST)
            wz_ghz.append(m["w_z"] / (2 * np.pi * 1e9))
            sat.append(bool(m["saturado"]))
        resultados[etiqueta] = wz_ghz
        saturados[etiqueta] = sat
        # Solo entran al ajuste los canales cuyo óptimo es identificable.
        x = [bw[c] for c, s in zip(CANALES, sat) if not s]
        y = [v for v, s in zip(wz_ghz, sat) if not s]
        p, k, r2 = ajustar_ley_potencia(x, y)
        ajustes[etiqueta] = (p, k, r2, len(x))
        print("   %-22s" % etiqueta
              + "".join("%11.2f%s" % (v, "*" if s else " ")
                        for v, s in zip(wz_ghz, sat))
              + "  p = %.2f  R2 = %.3f  (n=%d)" % (p, r2, len(x)))
    print("   (*) objetivo saturado: el óptimo no es identificable y se excluye "
          "del ajuste.")

    print("\nExponente p de  w_z* = k * BW^p  (p > 0 => la regla se mantiene):")
    ps = [ajustes[e][0] for e, _, _ in ARQUITECTURAS]
    print("   rango de p: %.2f a %.2f   (media %.2f)" % (min(ps), max(ps),
                                                         float(np.mean(ps))))
    print("   R2 mínimo : %.3f" % min(ajustes[e][2] for e, _, _ in ARQUITECTURAS))

    # --- Figura ------------------------------------------------------------ #
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))

    ax = axes[0]
    x = np.array([bw[c] for c in CANALES])
    xs = np.linspace(x.min() * 0.8, x.max() * 1.2, 50)
    for etiqueta, _, color in ARQUITECTURAS:
        y = np.array(resultados[etiqueta])
        sat = np.array(saturados[etiqueta])
        p, k = ajustes[etiqueta][0], ajustes[etiqueta][1]
        ax.plot(xs, k * xs ** p, color=color, lw=1.2, alpha=0.7)
        ax.plot(x[~sat], y[~sat], "o", color=color, ms=4.5, label=etiqueta)
        if sat.any():
            ax.plot(x[sat], y[sat], "o", color=color, ms=4.5, mfc="none")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"(a) Channel $-10$ dB bandwidth (GHz)")
    ax.set_ylabel(r"Optimal CTLE zero $f_z^{*}$ (GHz)")
    ax.legend(fontsize=6, loc="upper left")
    ax.grid(alpha=0.3, which="both")

    ax = axes[1]
    etiquetas = [e for e, _, _ in ARQUITECTURAS]
    colores = [c for _, _, c in ARQUITECTURAS]
    yy = np.arange(len(etiquetas))
    ax.barh(yy, [ajustes[e][0] for e in etiquetas], color=colores, alpha=0.85)
    for i, e in enumerate(etiquetas):
        ax.text(ajustes[e][0] + 0.02, i, "$R^2$=%.2f" % ajustes[e][2],
                va="center", fontsize=6.5)
    ax.set_yticks(yy)
    ax.set_yticklabels([e.split()[0] for e in etiquetas], fontsize=8)
    ax.invert_yaxis()
    ax.axvline(0.0, color="0.4", lw=0.9)
    ax.set_xlabel(r"(b) Fitted exponent $p$ in $f_z^{*}\propto \mathrm{BW}^{p}$")
    ax.grid(alpha=0.3, axis="x")

    fig.tight_layout()
    ruta_fig = os.path.join(DIR_FIG, "fig_generalizacion_wz.png")
    fig.savefig(ruta_fig, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    # --- Tabla ------------------------------------------------------------- #
    ruta_tab = os.path.join(DIR_TAB, "tabla_r2c4_generalizacion.csv")
    with open(ruta_tab, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["arquitectura"]
                   + ["%s (%s) BW=%.2fGHz" % (n, c, bw[(n, c)]) for n, c in CANALES]
                   + ["exponente_p", "k", "R2", "n_puntos_ajuste"])
        for etiqueta, _, _ in ARQUITECTURAS:
            p, k, r2, n_aj = ajustes[etiqueta]
            w.writerow([etiqueta]
                       + ["%.2f%s" % (v, "*" if s_ else "")
                          for v, s_ in zip(resultados[etiqueta], saturados[etiqueta])]
                       + ["%.3f" % p, "%.3f" % k, "%.4f" % r2, n_aj])

    print("\nEntregables:\n  - %s\n  - %s" % (ruta_fig, ruta_tab))


if __name__ == "__main__":
    main()
