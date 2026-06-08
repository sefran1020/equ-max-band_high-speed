"""
modelo_rlgc.py — Línea de transmisión con pérdidas dependientes de frecuencia.

Implementa el modelo físico realista pedido en la Recomendación 01:

    Parámetros por unidad de longitud (p.u.l.):
      R(f)  efecto pelicular  : R(f) = sqrt(R_dc^2 + R_skin(f)^2),  R_skin ~ sqrt(f)
      L     inductancia       : L = Z0 / v_p
      G(f)  pérdida dieléctrica: G(f) = omega * C * tan(delta)
      C     capacitancia      : C = 1 / (Z0 * v_p)

    Constante de propagación  : gamma(f) = sqrt( (R + jwL) (G + jwC) ) = alpha + j*beta
    Impedancia característica  : Zc(f)    = sqrt( (R + jwL) / (G + jwC) )

A diferencia de la línea RC difusiva ( H = 1/cosh(sqrt(s*tau)) ), aquí el canal
es DISPERSIVO: la atenuación crece con sqrt(f) (skin) y con f (dieléctrico), y la
velocidad de fase depende de la frecuencia.

Referencias de modelado (ya en references.bib del manuscrito):
  - Wang & Wang, "Modeling of Distributed RLC Interconnect and Transmission Line
    via Closed Forms and Recursive Algorithms", IEEE TVLSI, 2010.  [5152908]
  - Celik & Cangellaris, "Simulation of dispersive multiconductor transmission
    lines by Pade approximation...", IEEE T-MTT, 1996.            [554593]
"""

import numpy as np

from parametros import (
    C0, MU0, SIGMA_CU, NP_TO_DB, INCH,
    Substrato, Geometria, GEOMETRIA,
)


class LineaTransmision:
    """Línea de transmisión RLGC con R(f) (skin) y G(f) (dieléctrico)."""

    def __init__(self, substrato: Substrato, geometria: Geometria = GEOMETRIA):
        self.substrato = substrato
        self.geo = geometria

        # Velocidad de fase e inductancia/capacitancia p.u.l. a partir de Z0 y eps_eff
        self.v_p = C0 / np.sqrt(substrato.eps_eff)      # [m/s]
        self.L = geometria.z0 / self.v_p                # [H/m]
        self.C = 1.0 / (geometria.z0 * self.v_p)        # [F/m]

        # Resistencia DC p.u.l.  R_dc = 1 / (sigma * area)
        area = geometria.ancho_w * geometria.espesor_t
        self.R_dc = 1.0 / (SIGMA_CU * area)             # [Ohm/m]

        # Coeficiente de efecto pelicular:  R_skin(f) = k_skin * sqrt(f)
        #   k_skin = sqrt(pi * mu0 / sigma) / w   (resistencia superficial / ancho)
        self.k_skin = np.sqrt(np.pi * MU0 / SIGMA_CU) / geometria.ancho_w

    # ---- Parámetros p.u.l. dependientes de la frecuencia ------------------- #
    def R(self, f):
        """Resistencia serie p.u.l. [Ohm/m] con efecto pelicular."""
        f = np.asarray(f, dtype=float)
        r_skin = self.k_skin * np.sqrt(np.abs(f))
        return np.sqrt(self.R_dc**2 + r_skin**2)

    def G(self, f):
        """Conductancia de fuga p.u.l. [S/m] por pérdida dieléctrica."""
        f = np.asarray(f, dtype=float)
        return 2.0 * np.pi * f * self.C * self.substrato.tan_d

    # ---- Magnitudes derivadas --------------------------------------------- #
    def gamma(self, f):
        """Constante de propagación gamma(f) = alpha + j*beta [1/m]."""
        f = np.asarray(f, dtype=float)
        w = 2.0 * np.pi * f
        serie = self.R(f) + 1j * w * self.L
        deriv = self.G(f) + 1j * w * self.C
        return np.sqrt(serie * deriv)

    def zc(self, f):
        """Impedancia característica compleja Zc(f) [Ohm]."""
        f = np.asarray(f, dtype=float)
        w = 2.0 * np.pi * f
        serie = self.R(f) + 1j * w * self.L
        deriv = self.G(f) + 1j * w * self.C
        # Evitar 0/0 en DC (deriv->0): Zc(DC) tiende a infinito; se sustituye
        # por un valor grande pero finito solo para estética numérica.
        with np.errstate(divide="ignore", invalid="ignore"):
            z = np.sqrt(serie / deriv)
        return z

    def atenuacion_db_por_pulgada(self, f):
        """Atenuación alpha en dB/pulgada (entregable de la Recomendación 01)."""
        alpha = np.real(self.gamma(f))      # [Np/m]
        return alpha * NP_TO_DB * INCH      # [dB/inch]

    # ---- Función de transferencia ----------------------------------------- #
    def H(self, f, longitud, zs=None, zl=None):
        """
        Transferencia tensión->tensión del canal de longitud `longitud` [m].

        - Si zs y zl son None: línea adaptada, H = exp(-gamma*longitud) (aísla la
          pérdida/dispersión propia del canal; es el caso pedido para comparar
          el comportamiento difusivo vs dispersivo).
        - Si se dan zs (fuente) y zl (carga): se usa la matriz ABCD con
          reflexiones reales.
        """
        f = np.asarray(f, dtype=float)
        gl = self.gamma(f) * longitud

        if zs is None and zl is None:
            return np.exp(-gl)

        zs = self.geo.z0 if zs is None else zs
        zl = self.geo.z0 if zl is None else zl

        zc = self.zc(f)
        A = np.cosh(gl)
        B = zc * np.sinh(gl)
        C = np.sinh(gl) / zc
        D = np.cosh(gl)

        H = zl / (A * zl + B + zs * (C * zl + D))

        # Reparar DC (f=0): la línea es una R serie (R_dc*longitud) hacia la carga.
        f_flat = np.atleast_1d(f)
        dc = f_flat == 0.0
        if np.any(dc):
            r_series = self.R_dc * longitud
            H = np.atleast_1d(H).astype(complex)
            H[dc] = zl / (zl + zs + r_series)
        return H

    # ---- Utilidad: longitud para igualar una pérdida objetivo en Nyquist --- #
    def longitud_para_il(self, il_db_objetivo, f_nyquist):
        """
        Devuelve la longitud [m] cuya pérdida de inserción adaptada (|H| en dB)
        a `f_nyquist` iguala `il_db_objetivo` (valor positivo en dB).

        Permite una comparación JUSTA con el canal RC: misma pérdida en Nyquist,
        distinto perfil de dispersión.
        """
        alpha = float(np.real(self.gamma(np.array([f_nyquist])))[0])  # Np/m
        il_por_metro = alpha * NP_TO_DB                                # dB/m
        return il_db_objetivo / il_por_metro
