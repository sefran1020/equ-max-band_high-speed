"""
correr_ltspice.py - Orquestador de la validacion LTspice (end-to-end).

Localiza LTspice, corre los netlists en batch (-b) y encadena el post-proceso
en Python (Bode, cadena completa, ojo). Pensado para reproducir todo de una.

Uso:
    python correr_ltspice.py            # todo: bode + cadena RC + ojo
    python correr_ltspice.py --bode     # solo cruce de Bode (rc + ctle)
    python correr_ltspice.py --cadena   # solo TX->...->ojo (canal por defecto)
    python correr_ltspice.py --vectorfit  # regenerar disp_lines.lib
    python correr_ltspice.py --ltspice "C:\\ruta\\LTspice.exe"
"""

import argparse
import glob
import os
import subprocess
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))

CANDIDATOS = [
    r"C:\Users\%s\AppData\Local\Programs\ADI\LTspice\LTspice.exe" % os.environ.get("USERNAME", ""),
    r"C:\Program Files\ADI\LTspice\LTspice.exe",
    r"C:\Program Files\LTC\LTspiceXVII\XVIIx64.exe",
    r"C:\Program Files\LTC\LTspiceIV\scad3.exe",
]


def hallar_ltspice(override=None):
    if override:
        return override
    for c in CANDIDATOS:
        if os.path.exists(c):
            return c
    raise FileNotFoundError("No se encontro LTspice; pasa --ltspice RUTA")


def correr_netlist(exe, cir):
    ruta = os.path.join(AQUI, cir)
    print(">> LTspice -b %s" % cir)
    r = subprocess.run([exe, "-b", ruta], cwd=AQUI,
                       capture_output=True, text=True)
    raw = ruta[:-4] + ".raw"
    ok = os.path.exists(raw)
    print("   %s  (%s)" % ("OK" if ok else "FALLO", os.path.basename(raw)))
    if not ok:
        log = ruta[:-4] + ".log"
        if os.path.exists(log):
            print("   --- log ---")
            print(open(log, encoding="latin-1", errors="ignore").read()[-800:])
    return ok


def py(script, *a):
    print(">> python %s %s" % (script, " ".join(a)))
    subprocess.run([sys.executable, os.path.join(AQUI, script), *a], cwd=AQUI)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ltspice", default=None)
    ap.add_argument("--bode", action="store_true")
    ap.add_argument("--cadena", action="store_true")
    ap.add_argument("--vectorfit", action="store_true")
    ap.add_argument("--baud", default="4e9")
    ap.add_argument("--canal", default="rc")
    args = ap.parse_args()
    exe = hallar_ltspice(args.ltspice)
    print("LTspice:", exe)

    todo = not (args.bode or args.cadena or args.vectorfit)

    if args.vectorfit or todo:
        py("vectorfit_s2p.py")

    if args.bode or todo:
        correr_netlist(exe, "rc_escalera_N5.cir")
        correr_netlist(exe, "rc_pulse.cir")
        correr_netlist(exe, "rc_bode.cir")
        correr_netlist(exe, "ctle_bode.cir")
        correr_netlist(exe, "disp_bode.cir")     # dispersivo: solo frecuencia
        py("comparar_bode.py")

    if args.cadena or todo:
        # Ojo en transitorio: canal RC (modelo central). El ojo dispersivo se
        # queda en Python (skin ~ sqrt(f) no es racional para LTspice .tran).
        py("generar_tx_pwl.py", "--baud", args.baud, "--canal", "rc")
        correr_netlist(exe, "cadena_completa.cir")
        py("eye_desde_ltspice.py")
        py("comparar_ojos.py")

    print("\nListo. Figuras .png en", AQUI)


if __name__ == "__main__":
    main()
