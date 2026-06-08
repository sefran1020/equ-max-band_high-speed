"""
parametros.py — Configuración central de la Recomendación 01.

Recomendación 01 (Fase 1, Dra. Elena Ruiz, equipo senior):
    "Evolución del Modelo de Canal (Física Realista)": sustituir la línea RC
    teórica por un modelo de línea de transmisión con pérdidas dieléctricas
    (FR-4 / Megtron 6) y efecto pelicular ( R(f) y G(f) ), pasando de un canal
    puramente DIFUSIVO a uno DISPERSIVO real.

Este archivo NO contiene lógica de cálculo: solo presets físicos y el punto
de operación tomados del checkpoint del proyecto, para que un revisor pueda
auditar/ajustar los supuestos en un solo lugar.
"""

from dataclasses import dataclass

# --------------------------------------------------------------------------- #
# Constantes físicas (SI)
# --------------------------------------------------------------------------- #
C0 = 2.99792458e8        # velocidad de la luz en vacío [m/s]
MU0 = 4.0e-7 * 3.141592653589793  # permeabilidad del vacío [H/m]
SIGMA_CU = 5.8e7         # conductividad del cobre [S/m]

# Conversión de unidades
INCH = 0.0254            # 1 pulgada en metros
NP_TO_DB = 8.685889638   # 1 neper = 20*log10(e) dB


# --------------------------------------------------------------------------- #
# Substratos de PCB (entregable: comparar FR-4 vs Megtron 6)
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Substrato:
    """Material dieléctrico de la PCB."""
    nombre: str
    eps_r: float       # permitividad relativa (Dk)
    tan_d: float       # tangente de pérdidas (Df) -> pérdida dieléctrica

    @property
    def eps_eff(self) -> float:
        """Permitividad efectiva aprox. para microstrip (mitad aire/mitad Dk)."""
        return (self.eps_r + 1.0) / 2.0


# Valores representativos de hoja de datos (banda de ~1-10 GHz)
FR4 = Substrato(nombre="FR-4", eps_r=4.3, tan_d=0.020)
MEGTRON6 = Substrato(nombre="Megtron 6", eps_r=3.6, tan_d=0.004)

SUBSTRATOS = (FR4, MEGTRON6)


# --------------------------------------------------------------------------- #
# Geometría del conductor (microstrip 50 Ohm típico)
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Geometria:
    """Geometría de la traza de cobre."""
    z0: float = 50.0         # impedancia característica objetivo [Ohm]
    ancho_w: float = 150e-6  # ancho de la traza [m] (~6 mil)
    espesor_t: float = 35e-6 # espesor del cobre [m] (1 oz)


GEOMETRIA = Geometria()


# --------------------------------------------------------------------------- #
# Punto de operación (alineado con checkpoint.md / RESULTADOS BLOQUEADOS)
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class PuntoOperacion:
    """Punto de operación del enlace según el artículo."""
    baud_op: float = 4e9       # tasa de símbolo de operación [Baud] (PAM-4)
    f_nyquist: float = 2e9     # frecuencia de Nyquist [Hz] = baud/2
    tau_rc: float = 400e-12    # constante R_t*C_t del canal RC del artículo [s]


OPERACION = PuntoOperacion()


# --------------------------------------------------------------------------- #
# Parámetros de la rejilla de frecuencia / tiempo para las respuestas
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Rejilla:
    """Discretización temporal/frecuencial para FFT de respuesta al pulso."""
    n_muestras: int = 1 << 14   # 16384 puntos
    dt: float = 2.0e-12         # paso temporal [s] -> Fs = 500 GHz
    f_min_bode: float = 1e7     # 10 MHz, inicio del barrido de atenuación
    f_max_bode: float = 4e10    # 40 GHz, fin del barrido de atenuación
    n_bode: int = 800


REJILLA = Rejilla()


# --------------------------------------------------------------------------- #
# Jitter para la apertura horizontal del ojo (Fase α del plan de jitter)
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Jitter:
    """Presupuesto de jitter en el instante de muestreo del RX (en UI).

    Modelo dual-Dirac: jitter aleatorio gaussiano (RJ, RMS) + jitter determinista
    acotado (DJ pico-a-pico, agrupa DCD y PJ). La DDJ ya está presente por la ISI
    del canal y NO se inyecta aquí.
    """
    sigma_rj_ui: float = 0.010   # RJ: desviación RMS en UI (p. ej. 1 % de UI)
    dj_pp_ui: float = 0.05       # DJ: pico-a-pico en UI (DCD + PJ, dual-Dirac)


JITTER = Jitter()
