"""
modelo_rc.py — Canal RC distribuido (referencia del artículo actual).

Reproduce el canal que el manuscrito usa hoy, para servir de LÍNEA BASE frente
al modelo dispersivo de `modelo_rlgc.py`.

Modelo difusivo (línea RC distribuida, terminación abierta):

    H_rc(s) = 1 / cosh( sqrt( s * tau ) ),   tau = R_t * C_t

Es exactamente el modelo `cosh` validado contra la "escalera" en el notebook
`desarrollo_modelo_distribuido.ipynb` (R0). No tiene inductancia ni efecto
pelicular: la atenuación crece de forma difusiva (~exp(-sqrt(f))), sin el término
lineal en f de la pérdida dieléctrica.
"""

import numpy as np
import scipy.signal as signal

from parametros import NP_TO_DB, OPERACION


class CanalRC:
    """Canal RC distribuido difusivo:  H(s) = 1/cosh(sqrt(s*tau))."""

    def __init__(self, tau: float = OPERACION.tau_rc):
        self.tau = tau

    def H(self, f):
        """Transferencia compleja H(f)."""
        f = np.asarray(f, dtype=float)
        s = 1j * 2.0 * np.pi * f
        arg = np.sqrt(s * self.tau)
        # cosh(0)=1 -> H(0)=1 (DC pasa). Sin singularidades.
        return 1.0 / np.cosh(arg)

    def atenuacion_db(self, f):
        """Atenuación |H| en dB (el canal RC no tiene 'longitud' física)."""
        return 20.0 * np.log10(np.abs(self.H(f)))

    def f_3db(self, f_busqueda=None):
        """Frecuencia de -3 dB estimada por interpolación."""
        if f_busqueda is None:
            f_busqueda = np.linspace(1e7, 5e9, 20000)
        mag_db = self.atenuacion_db(f_busqueda)
        idx = np.where(mag_db <= -3.0)[0]
        if len(idx) == 0:
            return np.nan
        i = idx[0]
        if i == 0:
            return f_busqueda[0]
        # interpolación lineal alrededor del cruce de -3 dB
        f0, f1 = f_busqueda[i - 1], f_busqueda[i]
        m0, m1 = mag_db[i - 1], mag_db[i]
        return f0 + (-3.0 - m0) * (f1 - f0) / (m1 - m0)


class CanalRCEscalera:
    """Canal RC distribuido como ESCALERA de N celdas (la del notebook).

    Es exactamente `crear_canal_rc_distribuido(R, C, N)` del notebook
    `simulacion_canal_rcEjecutadoActual.ipynb` (espacio de estados), pero expuesto
    como H(f) para aplicarlo en el mismo pipeline en frecuencia. Sirve para
    VALIDAR el barrido contra la referencia del checkpoint (2.79/9.56/14.11 GBaud),
    que se obtuvo con esta escalera (N=5).

    Se precomputa la función de transferencia racional H(s)=num(s)/den(s) una vez
    (ss2tf), de modo que H(f) es un simple cociente de polinomios (rápido).
    """

    def __init__(self, R_total: float = 80.0, C_total: float = 5e-12, N: int = 5):
        self.R_total, self.C_total, self.N = R_total, C_total, N
        R_seg = R_total / N
        C_seg = C_total / N
        alpha = 1.0 / (R_seg * C_seg)

        A = np.zeros((N, N))
        for i in range(N - 1):
            A[i, i] = -2 * alpha
            A[i, i + 1] = alpha
            A[i + 1, i] = alpha
        A[N - 1, N - 1] = -1 * alpha   # extremo abierto
        B = np.zeros((N, 1)); B[0, 0] = alpha
        C = np.zeros((1, N)); C[0, N - 1] = 1.0
        D = np.zeros((1, 1))

        num, den = signal.ss2tf(A, B, C, D)
        self.num = np.asarray(num).ravel()
        self.den = np.asarray(den).ravel()

    def H(self, f):
        """Transferencia compleja H(f) de la escalera (DC -> 1)."""
        f = np.asarray(f, dtype=float)
        s = 1j * 2.0 * np.pi * f
        return np.polyval(self.num, s) / np.polyval(self.den, s)

    def atenuacion_db(self, f):
        return 20.0 * np.log10(np.abs(self.H(f)))
