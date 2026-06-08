"""
respuesta_pulso.py — Respuesta al pulso de un canal vía FFT.

Dado un objeto canal con método H(f), calcula la respuesta a un pulso
rectangular de un intervalo unitario (UI) usando rfft/irfft. Se usa para el
entregable "comparativa de respuestas al pulso" de la Recomendación 01.
"""

import numpy as np

from parametros import REJILLA, OPERACION


def eje_tiempo_frecuencia(rejilla=REJILLA):
    """Devuelve (t, f) coherentes para rfft/irfft."""
    n, dt = rejilla.n_muestras, rejilla.dt
    t = np.arange(n) * dt
    f = np.fft.rfftfreq(n, dt)
    return t, f


def pulso_entrada(t, ancho_ui=1, baud=OPERACION.baud_op, t0=None):
    """Pulso rectangular de `ancho_ui` intervalos unitarios, centrado en t0."""
    ui = 1.0 / baud
    if t0 is None:
        t0 = t[-1] * 0.15           # arranca al 15% de la ventana
    ancho = ancho_ui * ui
    x = np.where((t >= t0) & (t < t0 + ancho), 1.0, 0.0)
    return x, ui, t0


def respuesta_al_pulso(canal, longitud=None, rejilla=REJILLA, **hkw):
    """
    Calcula (t, x_in, y_out) para un pulso de 1 UI a través de `canal`.

    `canal.H(f)` o `canal.H(f, longitud, ...)` según el tipo de modelo.
    """
    t, f = eje_tiempo_frecuencia(rejilla)
    x, ui, t0 = pulso_entrada(t)

    X = np.fft.rfft(x)
    if longitud is None:
        H = canal.H(f)
    else:
        H = canal.H(f, longitud, **hkw)

    Y = X * H
    y = np.fft.irfft(Y, n=rejilla.n_muestras)
    return t, x, np.real(y)


def metricas_pulso(t, y):
    """Métricas escalares de una respuesta al pulso para tabular."""
    pico = float(np.max(y))
    i_pico = int(np.argmax(y))
    t_pico = float(t[i_pico])

    # Anchura a media altura (FWHM) -> medida de dispersión/ISI
    mitad = pico / 2.0
    cruces = np.where(y >= mitad)[0]
    if len(cruces) >= 2:
        fwhm = float(t[cruces[-1]] - t[cruces[0]])
    else:
        fwhm = float("nan")

    return {"pico": pico, "t_pico_s": t_pico, "fwhm_s": fwhm}
