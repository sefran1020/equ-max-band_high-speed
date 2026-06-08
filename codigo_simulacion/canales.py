"""
canales.py — Construye los canales (RC + líneas dispersivas) según el criterio
de comparación: igual pérdida en Nyquist ('equal_il') o longitud física fija
('fixed'). Devuelve cada canal como un proveedor de H(f).
"""

import numpy as np

from parametros import SUBSTRATOS, OPERACION, INCH
from modelo_rc import CanalRC, CanalRCEscalera
from modelo_rlgc import LineaTransmision

FIXED_LEN_INCH = 10.0   # longitud física común para el criterio 'fixed'


def construir_canales(criterio):
    """
    criterio: 'equal_il'  -> cada línea se recorta a la longitud que iguala la
                             pérdida del RC en Nyquist (comparación a igual IL).
              'fixed'     -> ambas líneas a FIXED_LEN_INCH (comparación a igual
                             geometría; penaliza al material de más pérdidas).
    Devuelve (canales, il_obj) donde canales[nombre] = {'H':func, 'long_in':float|None}.
    """
    rc = CanalRC()
    rc_esc = CanalRCEscalera(R_total=80.0, C_total=5e-12, N=5)
    il_obj = float(-20 * np.log10(np.abs(rc.H(np.array([OPERACION.f_nyquist]))[0])))

    # Dos referencias RC: cosh (exacta) y escalera N=5 (la del notebook/checkpoint).
    canales = {
        "RC (cosh)": {"H": (lambda f: rc.H(f)), "long_in": None},
        "RC N=5": {"H": (lambda f: rc_esc.H(f)), "long_in": None},
    }
    for sub in SUBSTRATOS:
        ln = LineaTransmision(sub)
        if criterio == "equal_il":
            L = ln.longitud_para_il(il_obj, OPERACION.f_nyquist)
        elif criterio == "fixed":
            L = FIXED_LEN_INCH * INCH
        else:
            raise ValueError(f"Criterio desconocido: {criterio}")
        canales[sub.nombre] = {
            "H": (lambda f, ln=ln, L=L: ln.H(f, L)),
            "long_in": L / INCH,
        }
    return canales, il_obj
