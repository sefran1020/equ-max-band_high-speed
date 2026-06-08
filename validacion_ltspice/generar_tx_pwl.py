"""
generar_tx_pwl.py - Genera los estimulos para cadena_completa.cir y el GEMELO
de referencia en Python (mismo PRBS/PAM-4/cadena que el articulo).

Produce, en esta carpeta:
  tx.pwl        : forma de onda PAM-4 de TX (escalera) para 'Vtx ... PWL file=tx.pwl'
  ruido.pwl     : ruido AWGN (sigma=0.10) referido a la entrada del RX
  tx_meta.inc   : .param Tui y .param TSTOP que consume cadena_completa.cir
  ref_python.npz: salida V(out) calculada por el pipeline (gemelo causal) +
                  simbolos y rejilla, para el cruce en eye_desde_ltspice.py

La cadena del gemelo replica EXACTAMENTE la de LTspice (mismos bloques LTI y
ruido identico muestra a muestra), con LPF de RX CAUSAL (Butterworth orden 3),
para que la comparacion LTspice<->Python sea limpia.

Uso:  python generar_tx_pwl.py            (canal RC N=5, 4 GBaud)
      python generar_tx_pwl.py --baud 8e9 --nsym 600 --canal fr4
"""

import argparse
import os
import sys

import numpy as np
import scipy.signal as sig

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from cadena_enlace import (generar_prbs, modular_pam4, aplicar_ffe,
                           ctle_response, C_POST_RC, W_Z_RC, CTLE_ADC, SIGMA_IN)
from modelo_rc import CanalRCEscalera, CanalRC
from modelo_rlgc import LineaTransmision
from parametros import SUBSTRATOS, INCH

AQUI = os.path.dirname(__file__)
NSPB = 125          # muestras por UI -> dt = Tui/125 (= 2 ps a 4 GBaud)
TR = 1e-12          # transicion de la escalera PWL


def H_canal(nombre):
    """Devuelve H(f) del canal elegido, alineado con canales.py."""
    if nombre == "rc":
        return CanalRCEscalera(80.0, 5e-12, 5).H
    if nombre == "rc_cosh":
        return CanalRC().H
    subs = {s.nombre.replace(" ", "").replace("-", "").lower(): s for s in SUBSTRATOS}
    clave = {"fr4": "fr4", "megtron6": "megtron6", "meg": "megtron6"}.get(nombre, nombre)
    sub = subs[clave]
    ln = LineaTransmision(sub)
    L = 10.0 * INCH
    return lambda f, ln=ln, L=L: ln.H(f, L)


def shift(x, k):
    """Retardo causal de k muestras (rellena con ceros al inicio)."""
    if k <= 0:
        return x.copy()
    return np.concatenate([np.zeros(k), x[:-k]])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--baud", type=float, default=4e9)
    ap.add_argument("--nsym", type=int, default=400)
    ap.add_argument("--canal", default="rc",
                    help="rc | rc_cosh | fr4 | megtron6")
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()

    Tui = 1.0 / args.baud
    dt = Tui / NSPB
    fc_rx = args.baud                      # fc del LPF de RX = Rs

    # --- Simbolos PAM-4 desde el MISMO PRBS del pipeline -------------------- #
    bits = generar_prbs(2 * args.nsym, orden=23)
    syms = modular_pam4(bits)[: args.nsym]
    nsym = len(syms)
    N = nsym * NSPB

    # --- Escalera PAM-4 en rejilla fina ------------------------------------ #
    tx_fine = np.repeat(syms, NSPB)

    # --- FFE (mismos taps normalizados que aplicar_ffe(s,0,1,c_post)) ------- #
    norm = abs(C_POST_RC) + 1.0
    ffe_fine = (1.0 * shift(tx_fine, NSPB) + C_POST_RC * shift(tx_fine, 2 * NSPB)) / norm

    # --- Canal (multiplicacion en frecuencia, = LTI en reposo) ------------- #
    Hf = H_canal(args.canal)
    f = np.fft.rfftfreq(N, dt)
    y_chan = np.fft.irfft(np.fft.rfft(ffe_fine) * Hf(f), n=N)

    # --- Ruido AWGN identico para ambos motores (sigma referido a RX) ------ #
    rng = np.random.default_rng(args.seed)
    ruido = rng.normal(0.0, SIGMA_IN, N)
    y_noisy = y_chan + ruido

    # --- CTLE (freq) + LPF de RX CAUSAL (Butterworth 3) -------------------- #
    y_ctle = np.fft.irfft(np.fft.rfft(y_noisy) * ctle_response(f, W_Z_RC, A_dc=CTLE_ADC), n=N)
    b, a = sig.butter(3, min(fc_rx / (0.5 / dt), 0.99), btype="low")
    y_out = sig.lfilter(b, a, y_ctle)      # CAUSAL (no filtfilt)

    t = np.arange(N) * dt

    # ====================== escribir tx.pwl (escalera) ===================== #
    with open(os.path.join(AQUI, "tx.pwl"), "w") as fh:
        fh.write("0 %.6f\n" % syms[0])
        for k in range(nsym):
            t0 = k * Tui
            fh.write("%.6e %.6f\n" % (t0 + Tui, syms[k]))           # mantener
            if k + 1 < nsym:
                fh.write("%.6e %.6f\n" % (t0 + Tui + TR, syms[k + 1]))  # saltar

    # ====================== escribir ruido.pwl ============================= #
    # Muestra a muestra en la rejilla fina -> identico al gemelo de Python.
    with open(os.path.join(AQUI, "ruido.pwl"), "w") as fh:
        lineas = ["%.6e %.6e" % (t[i], ruido[i]) for i in range(N)]
        fh.write("\n".join(lineas) + "\n")

    # ====================== tx_meta.inc ==================================== #
    tstop = nsym * Tui + 5 * Tui
    with open(os.path.join(AQUI, "tx_meta.inc"), "w") as fh:
        fh.write("* generado por generar_tx_pwl.py (no editar a mano)\n")
        fh.write(".param Tui=%.6e\n" % Tui)
        fh.write(".param TSTOP=%.6e\n" % tstop)

    # ====================== referencia para el cruce ======================= #
    np.savez(os.path.join(AQUI, "ref_python.npz"),
             t=t, vout=y_out, syms=syms, Tui=Tui, baud=args.baud,
             nspb=NSPB, canal=args.canal)

    print("Generado para canal=%s, %.0f GBaud, %d simbolos:" %
          (args.canal, args.baud / 1e9, nsym))
    print("  tx.pwl, ruido.pwl (%d pts), tx_meta.inc, ref_python.npz" % N)
    print("  Tui=%.1f ps  TSTOP=%.2f ns" % (Tui * 1e12, tstop * 1e9))
    print("Ahora:  LTspice -b cadena_completa.cir   y luego  python eye_desde_ltspice.py")


if __name__ == "__main__":
    main()
