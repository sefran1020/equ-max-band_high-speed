"""
r2c2_baseline_comparativa.py — Revisor 2, comentario 2.

    "The manuscript repeatedly claims a 'co-design framework', but no comparison
     is provided with existing equalizer optimization methodologies. Please
     clarify the differences from previous optimization or co-design approaches
     and discuss the actual novelty of the proposed framework."

El reparo es razonable: el artículo describe su procedimiento pero nunca lo mide
contra los criterios de ajuste habituales. Este script hace esa medición.

Un procedimiento de ajuste tiene dos partes separables —el CRITERIO que se
optimiza y la BÚSQUEDA que lo recorre— y conviene compararlas por separado,
porque solo una de las dos es propia de este trabajo:

  Criterios
    tasa   : tasa viable de la cadena completa con ruido referido a la entrada
             (el del artículo)
    mmse   : error cuadrático medio entre las muestras ecualizadas y los
             símbolos transmitidos, con ganancia ajustada por mínimos cuadrados
    zf     : forzado a cero — anula el primer post-cursor de la respuesta al
             pulso del enlace
    ojo    : apertura del peor ojo (margen vertical 3-sigma)

  Búsquedas
    LHS-30 : muestreo por hipercubo latino de 30 puntos (el del artículo)
    rejilla: búsqueda exhaustiva 13x13 sobre el mismo dominio

Todas las combinaciones se miden después con la MISMA vara —la tasa viable
promediada sobre varias realizaciones— y se anota el coste en evaluaciones de
la cadena, para poder juzgar qué aporta cada parte.

Los criterios de proxy (mmse, zf, ojo) se evalúan a una tasa de diseño por
canal, tomada como la tasa viable que el enlace alcanza con la ecualización
heredada (Tabla I): es el punto de operación al que un diseñador ajustaría.

El canal Megtron 6 'fixed' se excluye: su objetivo satura dentro de la rejilla
de tasa y el óptimo no es identificable (ver r2c4 y la nota de la Tabla IV).

Entregables en ./figuras/ y ./tablas/:
  - fig_baseline_criterios.png     : (a) tasa viable obtenida por combinación,
                                     (b) coste en evaluaciones de la cadena.
  - tabla_r2c2_baseline.csv        : (w_z, c_post), tasa viable y coste.

Uso:
    python revision01/r2c2_baseline_comparativa.py
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
    SPS, SIGMA_IN, CPOST_MIN, CPOST_MAX, WZ_MIN, WZ_MAX,
    generar_prbs, modular_pam4, aplicar_ffe, aplicar_canal_freq,
    aplicar_ctle_freq, aplicar_rx_lpf, alinear_full, evaluar_config,
    viable_interp, lhs_center,
)
from canales import construir_canales                           # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
DIR_FIG = os.path.join(AQUI, "figuras")
DIR_TAB = os.path.join(AQUI, "tablas")
DPI = 300

N_DFE = 2
N_SYM = 800
TASAS = np.linspace(2e9, 26e9, 13)          # rejilla de tasa del DoE del artículo
GB = TASAS / 1e9
N_LHS = 30
N_REJILLA = 13                              # 13x13 = 169 puntos

# Medición final común: tasa viable de FFE+CTLE con varias realizaciones.
N_REAL_MED = 10
N_SYM_MED = 1000
TASAS_MED = np.linspace(1e9, 32e9, 32)
GB_MED = TASAS_MED / 1e9

# (canal, criterio de comparación, tasa de diseño para los criterios de proxy)
CANALES = [("FR-4", "equal_il", 9e9),
           ("Megtron 6", "equal_il", 11e9),
           ("FR-4", "fixed", 14e9)]

CRITERIOS = ("tasa", "mmse", "zf", "ojo")
COMBINACIONES = [("LHS-30", "tasa"), ("LHS-30", "mmse"), ("LHS-30", "zf"),
                 ("LHS-30", "ojo"), ("grid 13x13", "tasa")]
COLORES = {"tasa": "#2b6cb0", "mmse": "#dd6b20", "zf": "#805ad5", "ojo": "#319795"}
CRIT_EN = {"tasa": "viable rate", "mmse": "MMSE", "zf": "ZF", "ojo": "eye opening"}


# --------------------------------------------------------------------------- #
# Criterios
# --------------------------------------------------------------------------- #
def _cadena_hasta_rx(H_func, sym, Rs, w_z, c_post, rng):
    """FFE -> canal -> (+ruido a la entrada) -> CTLE -> LPF."""
    ts = 1.0 / (Rs * SPS)
    ffe = aplicar_ffe(sym, 0.0, 1.0, c_post)
    y_chan = aplicar_canal_freq(np.repeat(ffe, SPS), ts, H_func)
    nin = rng.normal(0, SIGMA_IN, len(y_chan))
    return aplicar_rx_lpf(aplicar_ctle_freq(y_chan + nin, ts, w_z), ts, Rs), ts


def criterio_mmse(H_func, sym, Rs, w_z, c_post):
    """-MSE normalizado entre las muestras ecualizadas y los símbolos.

    La ganancia se ajusta por mínimos cuadrados antes de medir el error, porque
    un ecualizador que solo escala no es mejor ni peor por ello.
    """
    y, _ = _cadena_hasta_rx(H_func, sym, Rs, w_z, c_post, np.random.default_rng(2026))
    m, tx = alinear_full(y, sym, SPS)
    m, tx = m[100:], tx[100:]
    if len(m) < 50:
        return -1e9
    g = float(np.dot(m, tx) / np.dot(tx, tx))
    if abs(g) < 1e-12:
        return -1e9
    err = m / g - tx
    return -float(np.mean(err ** 2) / np.mean(tx ** 2))


def criterio_zf(H_func, sym, Rs, w_z, c_post):
    """-|p1/p0| de la respuesta al pulso del enlace (forzado a cero)."""
    ts = 1.0 / (Rs * SPS)
    x = np.zeros(256)
    x[16] = 1.0
    ffe = aplicar_ffe(x, 0.0, 1.0, c_post)
    y = aplicar_canal_freq(np.repeat(ffe, SPS), ts, H_func)
    y = aplicar_rx_lpf(aplicar_ctle_freq(y, ts, w_z), ts, Rs)
    mejor = None
    for ph in range(SPS):
        p = y[ph::SPS]
        ci = int(np.argmax(np.abs(p)))
        if mejor is None or abs(p[ci]) > abs(mejor[2]):
            mejor = (ph, ci, p[ci])
    ph, ci, cursor = mejor
    if abs(cursor) < 1e-12:
        return -1e9
    p = y[ph::SPS]
    p1 = p[ci + 1] if ci + 1 < len(p) else 0.0
    return -abs(float(p1) / float(cursor))


def criterio_ojo(H_func, sym, Rs, w_z, c_post):
    rng = np.random.default_rng(2026)
    return evaluar_config(H_func, sym, Rs, w_z, c_post, N_DFE, rng)["C"]["eye_min"]


def criterio_tasa(H_func, sym, Rs_diseno, w_z, c_post):
    rng = np.random.default_rng(2026)
    bers = [evaluar_config(H_func, sym, Rs, w_z, c_post, N_DFE, rng)["C"]["ber"]
            for Rs in TASAS]
    return viable_interp(bers, GB)


CRIT_FUN = {"tasa": criterio_tasa, "mmse": criterio_mmse,
            "zf": criterio_zf, "ojo": criterio_ojo}
# Evaluaciones de la cadena que consume una muestra del diseño con cada criterio.
CRIT_COSTE = {"tasa": len(TASAS), "mmse": 1, "zf": 1, "ojo": 1}


# --------------------------------------------------------------------------- #
# Búsquedas
# --------------------------------------------------------------------------- #
def muestras_lhs():
    dis = lhs_center(2, N_LHS, np.random.default_rng(42))
    return [(2 * np.pi * (WZ_MIN + dis[i, 0] * (WZ_MAX - WZ_MIN)),
             CPOST_MIN + dis[i, 1] * (CPOST_MAX - CPOST_MIN)) for i in range(N_LHS)]


def muestras_rejilla():
    wz = 2 * np.pi * np.linspace(WZ_MIN, WZ_MAX, N_REJILLA)
    cp = np.linspace(CPOST_MIN, CPOST_MAX, N_REJILLA)
    return [(float(a), float(b)) for a in wz for b in cp]


BUSQUEDAS = {"LHS-30": muestras_lhs, "grid 13x13": muestras_rejilla}


def optimizar(H_func, sym, Rs_diseno, busqueda, criterio):
    fun = CRIT_FUN[criterio]
    muestras = BUSQUEDAS[busqueda]()
    mejor = None
    for w_z, c_post in muestras:
        v = fun(H_func, sym, Rs_diseno, w_z, c_post)
        if mejor is None or v > mejor[0]:
            mejor = (v, w_z, c_post)
    coste = len(muestras) * CRIT_COSTE[criterio]
    return {"w_z": mejor[1], "c_post": mejor[2], "coste": coste}


# --------------------------------------------------------------------------- #
# Vara común: tasa viable de FFE+CTLE
# --------------------------------------------------------------------------- #
def medir_tasa_viable(H_func, w_z, c_post):
    """Tasa viable de FFE+CTLE y de la cadena completa con realimentación.

    Se miden las dos porque el criterio de ajuste puede ser indiferente para la
    etapa lineal y no serlo una vez que el DFE entra en el lazo.
    """
    bits_pool = generar_prbs(N_REAL_MED * 2 * N_SYM_MED + 5000, orden=23)
    vC, vD = [], []
    for r in range(N_REAL_MED):
        bsl = bits_pool[r * 2 * N_SYM_MED: r * 2 * N_SYM_MED + 2 * N_SYM_MED + 200]
        sym = modular_pam4(bsl)[:N_SYM_MED]
        rng = np.random.default_rng(1000 + r)
        bC, bD = [], []
        for Rs in TASAS_MED:
            res = evaluar_config(H_func, sym, Rs, w_z, c_post, N_DFE, rng)
            bC.append(res["C"]["ber"]); bD.append(res["D"]["ber"])
        vC.append(viable_interp(bC, GB_MED)); vD.append(viable_interp(bD, GB_MED))
    return (float(np.mean(vC)), float(np.std(vC)),
            float(np.mean(vD)), float(np.std(vD)))


def main():
    os.makedirs(DIR_FIG, exist_ok=True)
    os.makedirs(DIR_TAB, exist_ok=True)

    canales = {c: construir_canales(c)[0] for c in ("equal_il", "fixed")}
    sym = modular_pam4(generar_prbs(2 * N_SYM + 400, orden=23))[:N_SYM]

    print("== R2.2 — criterio y búsqueda frente a los ajustes habituales ==")
    print("Vara común: tasa viable de FFE+CTLE, %d realizaciones, umbral 1e-2."
          % N_REAL_MED)

    filas, resumen = [], {}
    for nombre, criterio_canal, Rs_dis in CANALES:
        H = canales[criterio_canal][nombre]["H"]
        print("\n-- %s (%s), tasa de diseño %.0f GBaud"
              % (nombre, criterio_canal, Rs_dis / 1e9))
        print("   %-12s %-7s %9s %9s %14s %15s %8s"
              % ("búsqueda", "criterio", "f_z (GHz)", "c_post",
                 "FFE+CTLE", "+DFE", "coste"))
        for busqueda, criterio in COMBINACIONES:
            r = optimizar(H, sym, Rs_dis, busqueda, criterio)
            mc, sc, md, sd = medir_tasa_viable(H, r["w_z"], r["c_post"])
            clave = (busqueda, criterio)
            resumen.setdefault(clave, []).append((mc, md))
            print("   %-12s %-7s %9.2f %9.3f  %6.2f ± %.2f  %6.2f ± %.2f %8d"
                  % (busqueda, criterio, r["w_z"] / (2 * np.pi * 1e9),
                     r["c_post"], mc, sc, md, sd, r["coste"]))
            filas.append([nombre, criterio_canal, busqueda, criterio,
                          "%.2f" % (r["w_z"] / (2 * np.pi * 1e9)),
                          "%.3f" % r["c_post"], "%.2f" % mc, "%.2f" % sc,
                          "%.2f" % md, "%.2f" % sd, r["coste"]])

    # --- Lectura agregada -------------------------------------------------- #
    # Las tasas absolutas difieren mucho entre canales (9 frente a 24 GBaud),
    # así que se normaliza cada canal contra la combinación del artículo antes
    # de promediar.
    ref = resumen[("LHS-30", "tasa")]
    rel = {}
    for clave, vals in resumen.items():
        rel[clave] = (float(np.mean([v[0] / r[0] for v, r in zip(vals, ref)])) * 100,
                      float(np.mean([v[1] / r[1] for v, r in zip(vals, ref)])) * 100)

    print("\nTasa viable relativa a la combinación del artículo (%), "
          "promedio de los tres canales:")
    print("   %-12s %-7s %12s %12s %10s"
          % ("búsqueda", "criterio", "FFE+CTLE", "+DFE", "coste"))
    for busqueda, criterio in COMBINACIONES:
        n = N_LHS if busqueda == "LHS-30" else N_REJILLA ** 2
        print("   %-12s %-7s %11.1f %11.1f %10d"
              % (busqueda, criterio, rel[(busqueda, criterio)][0],
                 rel[(busqueda, criterio)][1], n * CRIT_COSTE[criterio]))

    print("\nLectura:")
    print("   - La búsqueda exhaustiva mejora %+.1f %% la etapa lineal con %.1fx"
          " el coste: el muestreo LHS ya captura casi todo el óptimo."
          % (rel[("grid 13x13", "tasa")][0] - 100.0,
             N_REJILLA ** 2 / float(N_LHS)))
    for criterio in ("mmse", "ojo"):
        print("   - El criterio '%s' queda en %+.1f %% (lineal) y %+.1f %% (con DFE) "
              "con 1/%d del coste." % (criterio, rel[("LHS-30", criterio)][0] - 100.0,
                                       rel[("LHS-30", criterio)][1] - 100.0,
                                       CRIT_COSTE["tasa"]))
    print("   - El criterio 'zf' queda en %+.1f %% (lineal) y %+.1f %% (con DFE): "
          "anula ISI sin considerar el ruido que el refuerzo amplifica."
          % (rel[("LHS-30", "zf")][0] - 100.0, rel[("LHS-30", "zf")][1] - 100.0))

    # --- Figura ------------------------------------------------------------ #
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
    # Las claves internas están en español; las etiquetas de la figura van al
    # artículo y por tanto en inglés.
    etiquetas = ["%s\n%s" % (b, CRIT_EN[c]) for b, c in COMBINACIONES]
    xx = np.arange(len(COMBINACIONES))

    ax = axes[0]
    ax.bar(xx - 0.19, [rel[k][0] for k in COMBINACIONES], 0.36,
           color="#dd6b20", label="FFE+CTLE")
    ax.bar(xx + 0.19, [rel[k][1] for k in COMBINACIONES], 0.36,
           color="#319795", label="FFE+CTLE+DFE")
    ax.axhline(100.0, color="0.35", lw=1.0, ls="--")
    ax.set_xticks(xx)
    ax.set_xticklabels(etiquetas, fontsize=6.5)
    ax.set_ylabel("Viable rate, relative (%)")
    ax.set_xlabel("(a) Search and criterion")
    ax.set_ylim(0, 125)
    ax.legend(fontsize=6.5, loc="lower left")
    ax.grid(alpha=0.3, axis="y")

    ax = axes[1]
    costes = [(N_LHS if b == "LHS-30" else N_REJILLA ** 2) * CRIT_COSTE[c]
              for b, c in COMBINACIONES]
    ax.bar(xx, costes, 0.6, color=[COLORES[c] for _, c in COMBINACIONES])
    ax.set_yscale("log")
    ax.set_xticks(xx)
    ax.set_xticklabels(etiquetas, fontsize=6.5)
    ax.set_ylabel("Chain evaluations per design")
    ax.set_xlabel("(b) Design cost")
    ax.grid(alpha=0.3, axis="y", which="both")

    fig.tight_layout()
    ruta_fig = os.path.join(DIR_FIG, "fig_baseline_criterios.png")
    fig.savefig(ruta_fig, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    ruta_tab = os.path.join(DIR_TAB, "tabla_r2c2_baseline.csv")
    with open(ruta_tab, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["canal", "criterio_canal", "busqueda", "criterio",
                    "fz_GHz", "c_post", "FFE+CTLE_GBaud", "desv",
                    "+DFE_GBaud", "desv", "evaluaciones"])
        for fila in filas:
            w.writerow(fila)

    print("\nEntregables:\n  - %s\n  - %s" % (ruta_fig, ruta_tab))


if __name__ == "__main__":
    main()
