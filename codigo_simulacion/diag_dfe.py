"""
diag_dfe.py — Diagnóstico del DFE: ¿por qué diverge en el canal dispersivo?

Compara RC (funciona) vs FR-4 (falla) en el MISMO punto, imprimiendo taps,
ganancia (escala), BER de C (FFE+CTLE) vs D (+DFE), niveles y alineación.
No genera figuras; solo imprime. Ejecutar:  python diag_dfe.py
"""

import numpy as np

from cadena_enlace import (
    SPS, SIGMA_IN, W_Z_RC, C_POST_RC,
    generar_prbs, modular_pam4, aplicar_ffe, aplicar_canal_freq,
    aplicar_ctle_freq, aplicar_rx_lpf, alinear_full, estadisticas_niveles,
    evaluar_enlace,
)
from canales import construir_canales


def diagnosticar(nombre, H, Rs, sym, nt=2):
    sps = SPS
    ts = 1.0 / (Rs * sps)
    fc = Rs
    rng = np.random.default_rng(0)

    # Construir la forma de onda C (FFE+CTLE), igual que evaluar_config
    ffe = aplicar_ffe(sym, 0.0, 1.0, C_POST_RC)
    y_chan = aplicar_canal_freq(np.repeat(ffe, sps), ts, H)
    nin = rng.normal(0, SIGMA_IN, len(y_chan))
    yC = aplicar_rx_lpf(aplicar_ctle_freq(y_chan + nin, ts, W_Z_RC), ts, fc)

    # Taps DATA-AIDED estimados de los propios datos (misma fase de muestreo)
    m, tx = alinear_full(yC, sym, sps)
    m, tx = m[100:], tx[100:]
    N = len(m)
    Xreg = np.column_stack([tx[nt - i: N - i] for i in range(nt + 1)])
    h, *_ = np.linalg.lstsq(Xreg, m[nt:N], rcond=None)
    mu0, s0, c0 = estadisticas_niveles(m, tx)

    rC = evaluar_enlace(yC, sym, sps)
    rD = evaluar_enlace(yC, sym, sps, dfe_ntaps=nt)

    print(f"\n===== {nombre}  @ {Rs/1e9:.0f} GBaud =====")
    print(f"  LS data-aided h: cursor h0={h[0]:+.4f}  post-cursores={np.round(h[1:], 4)}"
          f"  (ratios={np.round(h[1:]/h[0], 4)})")
    print(f"  niveles tx -> mu0: {np.round(mu0, 3)}   cnt: {c0}")
    print(f"  C (FFE+CTLE):  BER={rC['ber']:.2e}  eye={rC['eye_min']:+.3f}")
    print(f"  D (+DFE):      BER={rD['ber']:.2e}  eye={rD['eye_min']:+.3f}")
    print(f"  m[:6] = {np.round(m[:6], 3)}   tx[:6] = {tx[:6]}")


def main():
    print("== Diagnóstico DFE: RC (ok) vs FR-4 (falla) ==")
    bits = generar_prbs(2 * 1000 + 400, orden=23)
    sym = modular_pam4(bits)[:1000]

    can_eq, _ = construir_canales("equal_il")
    can_fx, _ = construir_canales("fixed")

    for Rs in (1e9, 4e9):
        diagnosticar("RC (cosh)", can_eq["RC (cosh)"]["H"], Rs, sym)
        diagnosticar("FR-4 equal_il", can_eq["FR-4"]["H"], Rs, sym)
        diagnosticar("FR-4 fixed 10in", can_fx["FR-4"]["H"], Rs, sym)


if __name__ == "__main__":
    main()
