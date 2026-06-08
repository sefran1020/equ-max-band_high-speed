"""
lt_raw.py - Lector minimo de archivos .raw de LTspice (sin dependencias extra).

Soporta los casos que genera este kit:
  - Binario .tran (real): eje X (tiempo) en double de 8 bytes, resto en float
    de 4 bytes  (esquema por defecto de LTspice).
  - Binario .ac (complejo): cada valor es un par (real, imag) de doubles.
  - Variante "todo double" y formato ASCII (Values:) como respaldo.

El encabezado de LTspice en Windows suele estar en UTF-16LE; se detecta solo.

Uso:
    from lt_raw import leer_raw
    d = leer_raw("rc_escalera_N5.raw")
    t = d["x"]              # eje (tiempo o frecuencia)
    v = d["V(out)"]         # traza por nombre (complejo si es .ac)
"""

import re
import numpy as np


def _decodificar_encabezado(raw_bytes):
    """Devuelve (texto_encabezado, offset_datos, es_utf16)."""
    # Detectar UTF-16LE: el segundo byte de 'T' (Title) es 0x00.
    es_utf16 = len(raw_bytes) > 1 and raw_bytes[1] == 0
    for marca in (b"Binary:\n", b"Values:\n"):
        if es_utf16:
            marca_b = marca.decode("ascii").encode("utf-16-le")
        else:
            marca_b = marca
        idx = raw_bytes.find(marca_b)
        if idx != -1:
            fin_marca = idx + len(marca_b)
            cab = raw_bytes[:idx].decode("utf-16-le" if es_utf16 else "latin-1",
                                         errors="ignore")
            tipo = "binary" if b"Binary" in marca else "ascii"
            return cab, fin_marca, es_utf16, tipo
    raise ValueError("No se encontro 'Binary:' ni 'Values:' en el .raw")


def leer_raw(ruta):
    with open(ruta, "rb") as fh:
        rb = fh.read()
    cab, off, _es_utf16, tipo = _decodificar_encabezado(rb)

    nvars = int(re.search(r"No\. Variables:\s*(\d+)", cab).group(1))
    npts = int(re.search(r"No\. Points:\s*(\d+)", cab).group(1))
    flags = re.search(r"Flags:\s*(.*)", cab).group(1).lower()
    complejo = "complex" in flags

    # Nombres de variables (lineas 'idx  nombre  tipo' tras 'Variables:')
    nombres = []
    capturando = False
    for linea in cab.splitlines():
        if linea.strip().startswith("Variables:"):
            capturando = True
            continue
        if capturando:
            partes = linea.split()
            if len(partes) >= 2 and partes[0].isdigit():
                nombres.append(partes[1])
            if len(nombres) == nvars:
                break

    if tipo == "ascii":
        datos = _leer_ascii(rb[off:], nvars, npts, complejo)
    else:
        datos = _leer_binario(rb[off:], nvars, npts, complejo)

    salida = {"nombres": nombres, "npts": npts, "complejo": complejo}
    salida["x"] = np.abs(datos[0]) if complejo else datos[0].real
    for i, nom in enumerate(nombres):
        salida[nom] = datos[i] if complejo else datos[i].real
    return salida


def _leer_binario(buf, nvars, npts, complejo):
    cols = [np.empty(npts, dtype=complex) for _ in range(nvars)]
    if complejo:
        # Todos los valores: par de doubles (real, imag).
        arr = np.frombuffer(buf[: npts * nvars * 16], dtype="<f8")
        arr = arr.reshape(npts, nvars * 2)
        for i in range(nvars):
            cols[i] = arr[:, 2 * i] + 1j * arr[:, 2 * i + 1]
        return cols
    # Real: detectar esquema por tamano de bloque.
    bloque_mixto = 8 + 4 * (nvars - 1)      # x:double, resto:float
    bloque_double = 8 * nvars                # todo double
    disponible = len(buf)
    if disponible >= npts * bloque_mixto and \
       (disponible < npts * bloque_double or bloque_mixto == bloque_double):
        # Esquema mixto (por defecto en LTspice .tran)
        bruto = buf[: npts * bloque_mixto]
        registro = np.dtype([("x", "<f8")] +
                            [(f"v{i}", "<f4") for i in range(nvars - 1)])
        rec = np.frombuffer(bruto, dtype=registro)
        cols[0] = rec["x"].astype(complex)
        for i in range(1, nvars):
            cols[i] = rec[f"v{i-1}"].astype(complex)
    else:
        arr = np.frombuffer(buf[: npts * nvars * 8], dtype="<f8").reshape(npts, nvars)
        for i in range(nvars):
            cols[i] = arr[:, i].astype(complex)
    return cols


def _leer_ascii(buf, nvars, npts, complejo):
    texto = buf.decode("utf-16-le", errors="ignore") if buf[1:2] == b"\x00" \
        else buf.decode("latin-1", errors="ignore")
    nums = re.findall(r"[-+0-9.eE]+(?:[-+]\d+)?", texto)
    cols = [np.empty(npts, dtype=complex) for _ in range(nvars)]
    # En ASCII cada punto: idx, luego nvars valores (complejo: 'real,imag')
    it = iter(texto.split("\n"))
    p = 0
    for linea in it:
        toks = linea.replace("\t", " ").split()
        if not toks:
            continue
        # Primera columna de un punto empieza con el indice entero
        if toks[0].isdigit() and len(toks) >= 2:
            vals = toks[1:]
            base = 0
        else:
            vals = toks
            base = 0
        # acumular hasta tener nvars valores para este punto
        # (LTspice ascii pone un valor por linea tras el indice)
        # Implementacion simple: leer secuencialmente
        for v in vals:
            if complejo and "," in v:
                r, im = v.split(",")
                cols[base][p] = float(r) + 1j * float(im)
            else:
                cols[base][p] = float(v)
            base += 1
            if base == nvars:
                p += 1
                base = 0
                break
        if p >= npts:
            break
    return cols


if __name__ == "__main__":
    import sys
    d = leer_raw(sys.argv[1])
    print("variables:", d["nombres"])
    print("npts:", d["npts"], "complejo:", d["complejo"])
    print("x:", d["x"][:3], "...", d["x"][-1])
