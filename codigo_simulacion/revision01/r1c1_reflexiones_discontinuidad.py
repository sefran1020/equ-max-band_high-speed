"""
r1c1_reflexiones_discontinuidad.py — Revisor 1, comentario 1.

    "It is quite difficult to realize the practical dispersive physical channel
     in terms of RC parameters. Because practical channel introduces fading
     with number of multipaths."

Dos cosas que aclarar, y una que demostrar.

1. El canal dispersivo del artículo NO se realiza con parámetros RC. Ese es
   precisamente su resultado: la línea RC difusiva ( H = 1/cosh(sqrt(s·R·C)) )
   se sustituye por una línea de transmisión con parámetros RLGC dependientes de
   la frecuencia —efecto pelicular R(f) y pérdida dieléctrica G(f), ec. (3) del
   manuscrito—, y la contribución es medir cuánto se pierde al reutilizar en ella
   un ecualizador ajustado sobre el modelo RC.

2. En un enlace GUIADO en cobre no hay desvanecimiento por multitrayecto. No hay
   dispersión de trayectos ni estadística de canal: H(f) es determinista y
   estacionaria, fijada por la geometría y el laminado. El fenómeno análogo al
   multitrayecto son las REFLEXIONES en discontinuidades de impedancia —vías con
   muñón (*stub*), conectores, cambios de sección—, que sí generan réplicas
   retardadas del pulso. La diferencia esencial: esas réplicas son repetibles,
   no aleatorias, y por eso se combaten con diseño y ecualización y no con
   diversidad.

3. El manuscrito declaraba las reflexiones fuera de alcance sin cuantificarlas.
   Este script las cuantifica: modela dos discontinuidades canónicas sobre la
   misma traza FR-4 de 10 pulgadas y mide cuánta tasa viable cuestan.

Modelo
------
Cascada de matrices ABCD, terminada en Zs = Zl = 50 ohm (a diferencia de la
línea adaptada del artículo, aquí las reflexiones SÍ se propagan):

  - Tramo de línea de longitud l:
        [[cosh(gamma·l),        Zc·sinh(gamma·l)],
         [sinh(gamma·l)/Zc,     cosh(gamma·l)   ]]

  - Muñón de vía (*stub*) en circuito abierto, en derivación:
        Z_stub = Zc_v · coth(gamma_v · L_s)   ->   [[1, 0], [1/Z_stub, 1]]
    Resuena en cuarto de onda a  f_res = v_p / (4·L_s), donde la línea principal
    ve un cortocircuito y aparece el hundimiento característico de |S21|.

  - Tramo desadaptado (conector / zona de pads): un tramo corto de línea con
    impedancia característica distinta (35 ohm), que produce el eco doble a dos
    veces su retardo de tránsito.

Y la transferencia total, con reflexiones de fuente y carga incluidas:
        H = Zl / ( A·Zl + B + Zs·(C·Zl + D) )

Entregables en ./figuras/ y ./tablas/:
  - fig_reflexiones.png        : (a) respuesta al pulso con ecos, (b) |S21| con
                                 el hundimiento resonante, (c) tasa viable vs
                                 longitud del muñón.
  - tabla_r1c1_reflexiones.csv : tasa viable por caso y configuración.

Uso:
    python revision01/r1c1_reflexiones_discontinuidad.py
"""

import csv
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from parametros import FR4, GEOMETRIA, INCH, Geometria, Substrato   # noqa: E402
from modelo_rlgc import LineaTransmision                            # noqa: E402
from cadena_enlace import (                                         # noqa: E402
    W_Z_RC, C_POST_RC, generar_prbs, modular_pam4, evaluar_config,
    viable_interp,
)
from canales import construir_canales                               # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
DIR_FIG = os.path.join(AQUI, "figuras")
DIR_TAB = os.path.join(AQUI, "tablas")
DPI = 300

LONGITUD_IN = 10.0                     # traza FR-4 del criterio 'fixed'
POS_DISC_IN = 5.0                      # discontinuidad a mitad de traza
MIL = 0.0254e-3                        # 1 mil en metros

STUBS_MIL = (0.0, 20.0, 40.0, 60.0, 90.0, 120.0)
Z_DESADAPTADO = 35.0                   # tramo de conector / pads [ohm]
LARGO_DESADAPTADO_IN = 0.5

# Barrido de tasa (mismos ajustes que el barrido del artículo, con menos
# realizaciones porque aquí interesa la diferencia entre casos, no el valor
# absoluto con su desviación publicada).
N_REAL = 15
N_SYM = 1000
N_DFE = 2
TASAS = np.linspace(1e9, 32e9, 32)
GB = TASAS / 1e9


# --------------------------------------------------------------------------- #
# Línea del muñón: una vía está ENTERRADA y ve el dieléctrico completo, no el
# eps_eff de microstrip. El modelo calcula eps_eff = (eps_r + 1)/2, así que se
# instancia con eps_r' = 2·eps_r - 1 para que su eps_eff resulte igual al eps_r
# del laminado. Se reutiliza así la misma matemática RLGC ya validada.
# --------------------------------------------------------------------------- #
SUB_VIA = Substrato(nombre="%s (via)" % FR4.nombre,
                    eps_r=2.0 * FR4.eps_r - 1.0, tan_d=FR4.tan_d)


def abcd_linea(linea, f, longitud):
    """Matriz ABCD de un tramo de línea, como cuatro arreglos sobre f."""
    gl = linea.gamma(f) * longitud
    zc = linea.zc(f)
    ch, sh = np.cosh(gl), np.sinh(gl)
    return ch, zc * sh, sh / zc, ch


def abcd_derivacion(y_shunt):
    """Matriz ABCD de una admitancia en derivación."""
    unos = np.ones_like(y_shunt)
    ceros = np.zeros_like(y_shunt)
    return unos, ceros, y_shunt, unos


def cascada(*matrices):
    """Producto ordenado de matrices ABCD dadas como tuplas de arreglos."""
    A, B, C, D = matrices[0]
    for (a, b, c, d) in matrices[1:]:
        A, B, C, D = (A * a + B * c, A * b + B * d,
                      C * a + D * c, C * b + D * d)
    return A, B, C, D


def h_desde_abcd(A, B, C, D, z0=50.0):
    """S21 de la red, referida a z0 en ambos puertos.

    Se usa S21 —y no la transferencia tensión->tensión, que incluye el divisor
    2:1 de fuente y carga— para que el caso sin discontinuidad sea directamente
    comparable con el canal adaptado del artículo ( H = exp(-gamma·l) ): en una
    línea perfectamente adaptada ambas expresiones coinciden. Es además la misma
    magnitud que el artículo exporta a Touchstone.
    """
    return 2.0 / (A + B / z0 + C * z0 + D)


class CanalConDiscontinuidad:
    """Traza FR-4 con una discontinuidad opcional a mitad de camino."""

    def __init__(self, stub_mil=0.0, desadaptado=False):
        self.linea = LineaTransmision(FR4)
        self.stub_mil = stub_mil
        self.desadaptado = desadaptado
        self.l1 = POS_DISC_IN * INCH
        self.l2 = (LONGITUD_IN - POS_DISC_IN) * INCH

        if stub_mil > 0:
            self.via = LineaTransmision(SUB_VIA)
            self.largo_stub = stub_mil * MIL
            self.f_res = self.via.v_p / (4.0 * self.largo_stub)
        else:
            self.via = None
            self.f_res = None

        if desadaptado:
            geo_d = Geometria(z0=Z_DESADAPTADO, ancho_w=GEOMETRIA.ancho_w,
                              espesor_t=GEOMETRIA.espesor_t)
            self.linea_d = LineaTransmision(FR4, geo_d)
            self.largo_d = LARGO_DESADAPTADO_IN * INCH
            # El tramo desadaptado sustituye parte de la traza, no la alarga.
            self.l2 -= self.largo_d
        else:
            self.linea_d = None

    def H(self, f):
        f = np.asarray(f, dtype=float)
        seguro = np.where(f == 0.0, 1.0, f)     # se repara el DC al final

        tramos = [abcd_linea(self.linea, seguro, self.l1)]
        if self.via is not None:
            gl = self.via.gamma(seguro) * self.largo_stub
            z_stub = self.via.zc(seguro) / np.tanh(gl)   # circuito abierto
            tramos.append(abcd_derivacion(1.0 / z_stub))
        if self.linea_d is not None:
            tramos.append(abcd_linea(self.linea_d, seguro, self.largo_d))
        tramos.append(abcd_linea(self.linea, seguro, self.l2))

        H = h_desde_abcd(*cascada(*tramos))

        dc = np.atleast_1d(f) == 0.0
        if np.any(dc):
            H = np.atleast_1d(H).astype(complex)
            # En DC la red es una resistencia serie entre dos cargas de 50 ohm;
            # el muñón en circuito abierto no conduce y no interviene.
            r_serie = self.linea.R_dc * (self.l1 + self.l2)
            H[dc] = 2.0 * 50.0 / (2.0 * 50.0 + r_serie)
        return H


def canal_referencia():
    """Traza de 10 in sin discontinuidad, terminada en 50 ohm."""
    return CanalConDiscontinuidad(stub_mil=0.0, desadaptado=False)


# --------------------------------------------------------------------------- #
# Respuesta al pulso (para ver los ecos) y barrido de tasa
# --------------------------------------------------------------------------- #
def respuesta_al_pulso(H_func, n=1 << 14, dt=2e-12, ancho_ps=40.0):
    """Excitación con un pulso estrecho (40 ps) para separar las réplicas.

    Un pulso de un intervalo de símbolo completo enmascararía el eco del tramo
    desadaptado, que llega a 2·tau_d ~ 138 ps del cursor.
    """
    t = np.arange(n) * dt
    x = np.zeros(n)
    x[: max(1, int(ancho_ps * 1e-12 / dt))] = 1.0
    f = np.fft.rfftfreq(n, dt)
    y = np.fft.irfft(np.fft.rfft(x) * H_func(f), n=n)
    return t, y


def tasa_viable(H_func, etiqueta):
    bits_pool = generar_prbs(N_REAL * 2 * N_SYM + 5000, orden=23)
    v_lin, v_dfe = [], []
    for r in range(N_REAL):
        bsl = bits_pool[r * 2 * N_SYM: r * 2 * N_SYM + 2 * N_SYM + 200]
        sym = modular_pam4(bsl)[:N_SYM]
        rng = np.random.default_rng(1000 + r)      # semillas pareadas del artículo
        bl, bd = [], []
        for Rs in TASAS:
            res = evaluar_config(H_func, sym, Rs, W_Z_RC, C_POST_RC, N_DFE, rng)
            bl.append(res["C"]["ber"]); bd.append(res["D"]["ber"])
        v_lin.append(viable_interp(bl, GB))
        v_dfe.append(viable_interp(bd, GB))
    print("   %-28s FFE+CTLE %5.2f ± %.2f    +DFE %5.2f ± %.2f  GBaud"
          % (etiqueta, np.mean(v_lin), np.std(v_lin),
             np.mean(v_dfe), np.std(v_dfe)))
    return (float(np.mean(v_lin)), float(np.std(v_lin)),
            float(np.mean(v_dfe)), float(np.std(v_dfe)))


# --------------------------------------------------------------------------- #
# Figura
# --------------------------------------------------------------------------- #
def graficar(pulsos, respuestas, barrido, ruta):
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.6))

    # (a) Se dibuja el RESIDUO respecto de la traza uniforme: aísla las réplicas
    #     retardadas, que a escala del cursor serían invisibles.
    ax = axes[0]
    (_, t_ref, y_ref, _) = pulsos[0]
    t0 = t_ref[int(np.argmax(y_ref))]
    ax.plot((t_ref - t0) * 1e12, y_ref / y_ref.max() * 0.02, color="0.6",
            lw=1.2, label="cursor (scaled $\\times 0.02$)")
    for (etiqueta, t, y, color) in pulsos[1:]:
        ax.plot((t - t0) * 1e12, y - y_ref, color=color, lw=1.6, label=etiqueta)
    ax.axhline(0.0, color="0.8", lw=0.8)
    # Llegada prevista del eco doble del tramo desadaptado: 2 x retardo de tránsito.
    tau_d = LARGO_DESADAPTADO_IN * INCH / LineaTransmision(FR4).v_p
    ax.axvline(2e12 * tau_d, color="#805ad5", ls=":", lw=1.2)
    ax.annotate(r"$2\tau_d$", xy=(2e12 * tau_d, 0), xytext=(4, 6),
                textcoords="offset points", color="#805ad5", fontsize=9)
    ax.set_xlim(-100, 350)
    ax.set_xlabel("(a) Time after cursor (ps)")
    ax.set_ylabel("Reflected replica (V)")
    ax.legend(fontsize=7.5)
    ax.grid(alpha=0.3)

    ax = axes[1]
    for (etiqueta, f, s21, color) in respuestas:
        ax.plot(f / 1e9, 20 * np.log10(np.abs(s21)), color=color, lw=1.6,
                label=etiqueta)
    ax.set_xlim(0, 40)
    ax.set_ylim(-60, 2)
    ax.set_xlabel("(b) Frequency (GHz)")
    ax.set_ylabel(r"$|S_{21}|$ (dB)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[2]
    stubs, lin, lin_e, dfe, dfe_e = barrido
    ax.errorbar(stubs, lin, yerr=lin_e, color="#dd6b20", lw=2, marker="o",
                ms=5, capsize=3, label="FFE+CTLE")
    ax.errorbar(stubs, dfe, yerr=dfe_e, color="#319795", lw=2, marker="s",
                ms=5, capsize=3, label="FFE+CTLE+DFE")
    ax.set_xlabel("(c) Via-stub length (mil)")
    ax.set_ylabel("Viable rate (GBaud)   [BER $\\leq 10^{-2}$]")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(ruta, dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def main():
    os.makedirs(DIR_FIG, exist_ok=True)
    os.makedirs(DIR_TAB, exist_ok=True)

    print("== R1.1 — reflexiones por discontinuidad (FR-4, %.0f in) ==" % LONGITUD_IN)

    ref = canal_referencia()
    stub90 = CanalConDiscontinuidad(stub_mil=90.0)
    desad = CanalConDiscontinuidad(desadaptado=True)

    # Verificación: sin discontinuidad, la cascada ABCD debe reducirse al canal
    # adaptado del artículo. Si esto falla, el modelo nuevo no es comparable.
    canales_art, _ = construir_canales("fixed")
    f_chk = np.array([1e9, 2e9, 5e9, 10e9, 20e9])
    d_max = float(np.max(np.abs(
        20 * np.log10(np.abs(ref.H(f_chk)))
        - 20 * np.log10(np.abs(canales_art[FR4.nombre]["H"](f_chk))))))
    print("Reducción al canal del artículo sin discontinuidad: "
          "desviación máxima %.4f dB hasta 20 GHz" % d_max)
    if d_max > 0.05:
        raise SystemExit("ERROR: el caso de referencia no reproduce el canal publicado.")

    print("Resonancia del muñón de 90 mil: %.1f GHz" % (stub90.f_res / 1e9))
    print("Tramo desadaptado: %.1f ohm sobre %.2f in "
          "(eco a 2x el retardo de tránsito)"
          % (Z_DESADAPTADO, LARGO_DESADAPTADO_IN))

    # --- (a) respuesta al pulso: los ecos son el análogo del multitrayecto -- #
    pulsos = []
    for etiqueta, canal, color in (
            ("uniform trace", ref, "#2b6cb0"),
            ("90 mil via stub", stub90, "#c05621"),
            ("35 $\\Omega$ section", desad, "#805ad5")):
        t, y = respuesta_al_pulso(canal.H)
        pulsos.append((etiqueta, t, y, color))

    # --- (b) |S21|: el hundimiento resonante del muñón --------------------- #
    f = np.linspace(1e7, 40e9, 1600)
    respuestas = [(etiqueta, f, canal.H(f), color)
                  for etiqueta, canal, color in (
                      ("uniform trace", ref, "#2b6cb0"),
                      ("90 mil via stub", stub90, "#c05621"),
                      ("35 $\\Omega$ section", desad, "#805ad5"))]

    # --- (c) coste en tasa viable ------------------------------------------ #
    print("\nTasa viable por caso (%d realizaciones):" % N_REAL)
    filas = []
    stubs, lin, lin_e, dfe, dfe_e = [], [], [], [], []
    for s_mil in STUBS_MIL:
        canal = CanalConDiscontinuidad(stub_mil=s_mil)
        etiqueta = "sin muñón" if s_mil == 0 else "muñón %.0f mil" % s_mil
        ml, sl, md, sd = tasa_viable(canal.H, etiqueta)
        stubs.append(s_mil); lin.append(ml); lin_e.append(sl)
        dfe.append(md); dfe_e.append(sd)
        filas.append(["via_stub", "%.0f mil" % s_mil,
                      "" if canal.f_res is None else "%.1f" % (canal.f_res / 1e9),
                      ml, sl, md, sd])

    ml, sl, md, sd = tasa_viable(desad.H, "tramo %.0f ohm" % Z_DESADAPTADO)
    filas.append(["tramo_desadaptado", "%.0f ohm / %.2f in"
                  % (Z_DESADAPTADO, LARGO_DESADAPTADO_IN), "", ml, sl, md, sd])

    # --- Resumen para la respuesta a los revisores -------------------------- #
    base_lin, base_dfe = lin[0], dfe[0]
    print("\nCoste relativo a la traza uniforme:")
    for i, s_mil in enumerate(STUBS_MIL[1:], start=1):
        print("   muñón %3.0f mil : FFE+CTLE %+6.1f %%   +DFE %+6.1f %%"
              % (s_mil, 100 * (lin[i] - base_lin) / base_lin,
                 100 * (dfe[i] - base_dfe) / base_dfe))
    print("   tramo %.0f ohm : FFE+CTLE %+6.1f %%   +DFE %+6.1f %%"
          % (Z_DESADAPTADO, 100 * (ml - base_lin) / base_lin,
             100 * (md - base_dfe) / base_dfe))

    ruta_fig = os.path.join(DIR_FIG, "fig_reflexiones.png")
    ruta_tab = os.path.join(DIR_TAB, "tabla_r1c1_reflexiones.csv")
    graficar(pulsos, respuestas, (stubs, lin, lin_e, dfe, dfe_e), ruta_fig)
    with open(ruta_tab, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["discontinuidad", "parametro", "f_resonancia_GHz",
                    "FFE+CTLE_GBaud", "desv", "+DFE_GBaud", "desv"])
        for fila in filas:
            w.writerow(fila[:3] + ["%.2f" % v for v in fila[3:]])
    print("\nEntregables:\n  - %s\n  - %s" % (ruta_fig, ruta_tab))


if __name__ == "__main__":
    main()
