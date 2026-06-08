"""
tabla_operacion.py — Tabla de ojo y BER en el PUNTO DE OPERACIÓN (4 GBaud).

Genera, con el MISMO pipeline validado (`cadena_enlace`) que el barrido de tasa y
el re-DoE, la apertura del peor ojo (V) y el BER por escenario (A/B/C/D) y canal
(RC N=5, FR-4, Megtron 6; criterio equal_il), promediando 30 realizaciones.

Sirve para llenar las tablas R2 (integridad de señal) y R3 (fiabilidad) de
`resultados.tex` de forma COHERENTE con la tasa viable R1 (evita números de otro
pipeline que la contradigan).

Uso:  python tabla_operacion.py
"""

import numpy as np

from cadena_enlace import (
    generar_prbs, modular_pam4, evaluar_config, W_Z_RC, C_POST_RC,
)
from canales import construir_canales

RS = 4e9            # punto de operación (Nyquist 2 GHz)
N_SYM = 1500
N_REAL = 30
N_DFE = 2
SEL = ["RC N=5", "FR-4", "Megtron 6"]
ESC = {"A": "A", "B": "B", "C": "C", "D": "D"}


def fmt_ber(b):
    if b < 1e-12:
        return "<1e-12"
    return f"{b:.1e}"


def main():
    canales, _ = construir_canales("equal_il")
    bits_pool = generar_prbs(N_REAL * 2 * N_SYM + 5000, orden=23)

    print(f"== Tabla de operación @ {RS/1e9:.0f} GBaud "
          f"(EQ heredada del RC, {N_REAL} realizaciones) ==\n")
    print(f"{'canal':<10} {'esc':<4} {'apertura_ojo_V':>14} {'BER(mediana)':>14}")
    print("-" * 46)

    for nb in SEL:
        H = canales[nb]["H"]
        accB = {k: [] for k in "ABCD"}
        accE = {k: [] for k in "ABCD"}
        for r in range(N_REAL):
            bsl = bits_pool[r * 2 * N_SYM: r * 2 * N_SYM + 2 * N_SYM + 200]
            sym = modular_pam4(bsl)[:N_SYM]
            rng = np.random.default_rng(1000 + r)
            res = evaluar_config(H, sym, RS, W_Z_RC, C_POST_RC, N_DFE, rng)
            for k in "ABCD":
                accB[k].append(res[k]["ber"])
                accE[k].append(res[k]["eye_min"])
        for k in "ABCD":
            ber_med = float(np.median(accB[k]))
            eye_mean = float(np.mean(accE[k]))
            print(f"{nb:<10} {k:<4} {eye_mean:>+14.3f} {fmt_ber(ber_med):>14}")
        print("-" * 46)


if __name__ == "__main__":
    main()
