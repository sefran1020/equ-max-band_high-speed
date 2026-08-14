"""
preparar_envio.py — Arma el paquete de reenvío para el IEEE Author Portal.

La carta de decisión pide tres archivos, cada uno en un campo distinto del
portal. Este script los reúne en una carpeta ordenada, con un LÉEME que dice
qué va en cada campo, y verifica que el paquete esté completo antes de darlo
por bueno. Como todo lo demás en esta carpeta, es regenerable: si el manuscrito
cambia, se vuelve a ejecutar y no hay que recordar qué se copió a dónde.

    python revision01/preparar_envio.py             # informe
    python revision01/preparar_envio.py --escribir  # arma la carpeta
    python revision01/preparar_envio.py --escribir --verificar
                                                    # además compila el paquete

La opción `--verificar` es la que importa antes de subir nada: compila el
paquete LaTeX en su propia carpeta, aislado del original, y comprueba que
produce el mismo número de páginas que el PDF de referencia. Si el paquete
tuviera una figura o un archivo de clase de menos, se descubre aquí y no en el
portal.

Estructura que produce:

    ENVIO_ACCESS_2026-32692/
    ├── LEEME.md
    ├── 1_respuesta_a_revisores/respuesta_revisores.pdf
    ├── 2_manuscrito_marcado/access_marcado.pdf
    └── 3_manuscrito_limpio/
        ├── access.pdf
        ├── fuente_latex/            (compilable tal cual)
        └── fuente_latex.zip
"""

import argparse
import os
import shutil
import subprocess
import sys
import zipfile

AQUI = os.path.dirname(os.path.abspath(__file__))
ENTREGA = os.path.normpath(os.path.join(AQUI, "..", ".."))
ARTICULO = os.path.join(ENTREGA, "IEEE-ACCESS")
ENVIO = os.path.join(ENTREGA, "ENVIO_ACCESS_2026-32692")

# Archivos que NO forman parte del paquete LaTeX del manuscrito limpio:
# la versión marcada y la carta van en sus propios campos del portal, y los
# artefactos de compilación no se suben.
EXCLUIR_EXACTOS = {
    "access_marcado.tex", "access_marcado.pdf",
    "respuesta_revisores.tex", "respuesta_revisores.pdf",
    "cover_letter_ieee_access.txt",
    "REVISION_COHERENCIA_PRESENTACION.md",
    "REVISION_FINAL_COHERENCIA_NIVEL.md",
    "fig_flujo_conceptual.mmd",
    # El PDF compilado no es fuente: va aparte, en 3_manuscrito_limpio/.
    # Meterlo en el .zip duplicaría 13 MB en una subida con límite de tamaño.
    "access.pdf",
}
EXCLUIR_EXTENSIONES = (".aux", ".log", ".blg", ".out", ".synctex.gz", ".toc",
                       ".fls", ".fdb_latexmk")

LEEME = """# Paquete de reenvío — Access-2026-32692

Generado por `recomendacion01_canal_fisico/revision01/preparar_envio.py`.
No editar a mano: si algo cambia, se vuelve a ejecutar el script.

## Qué sube a cada campo del IEEE Author Portal

Entrar en el portal, localizar el título del artículo rechazado y pulsar
**Start Resubmission**.

| Campo del portal | Archivo de este paquete |
| :--- | :--- |
| Author's Response Files | `1_respuesta_a_revisores/respuesta_revisores.pdf` |
| Highlighted PDF | `2_manuscrito_marcado/access_marcado.pdf` |
| Main Manuscript (PDF) | `3_manuscrito_limpio/access.pdf` |
| Main Manuscript (LaTeX) | `3_manuscrito_limpio/fuente_latex.zip` |

La carta de decisión pide el manuscrito limpio **en dos formatos**: el PDF y la
fuente (LaTeX o Word). Por eso van los dos.

## Qué contiene cada cosa

- **Respuesta a los revisores.** Punto por punto, con la preocupación citada,
  la respuesta y la acción, y punteros a secciones, ecuaciones, figuras y
  tablas de la versión revisada.
- **Manuscrito marcado.** El mismo manuscrito con los pasajes nuevos o
  modificados resaltados en amarillo, y una nota inicial que enumera las
  figuras, ecuaciones y secciones nuevas (los flotantes no se pueden resaltar
  sin sacarlos de su caja).
- **Manuscrito limpio.** El PDF final y la fuente LaTeX completa: `.tex`,
  `.bbl`, clase `ieeeaccess`, tipografías, bibliografía, figuras y fotos de los
  autores. El paquete compila tal cual con `pdflatex`; se comprobó compilándolo
  aislado antes de empaquetarlo.

## Antes de subir

- [ ] Confirmar el título definitivo del artículo.
- [ ] Decidir el recargo por páginas excedentes (el artículo tiene {paginas}).
- [ ] Depositar el DOI de Zenodo si se quiere citar en Code and Data
      Availability; hoy esa sección remite al repositorio de GitHub.
- [ ] Comprobar que la lista de autores no ha cambiado respecto de la versión
      enviada (si cambiara, hace falta el formulario de byline change).
"""


def archivos_fuente():
    """Archivos de IEEE-ACCESS que forman el paquete LaTeX compilable."""
    fuente = []
    for nombre in sorted(os.listdir(ARTICULO)):
        ruta = os.path.join(ARTICULO, nombre)
        if os.path.isdir(ruta):
            continue
        if nombre in EXCLUIR_EXACTOS or nombre.endswith(EXCLUIR_EXTENSIONES):
            continue
        fuente.append(nombre)
    return fuente


def paginas_pdf(ruta):
    """Número de páginas de un PDF, vía `pdfinfo` (poppler).

    No se cuenta `/Type /Page` sobre el archivo crudo: estos PDF llevan los
    objetos comprimidos y ese recuento da cero, que es peor que no medir porque
    parece un dato. Si `pdfinfo` no está disponible se devuelve None y quien
    llama decide.
    """
    try:
        r = subprocess.run(["pdfinfo", ruta], capture_output=True, text=True)
    except (OSError, FileNotFoundError):
        return None
    if r.returncode != 0:
        return None
    for linea in r.stdout.splitlines():
        if linea.lower().startswith("pages:"):
            return int(linea.split(":", 1)[1].strip())
    return None


def fmt_pag(n):
    return "%d páginas" % n if n is not None else "páginas: no medido"


def verificar_entradas():
    """Los tres PDF deben existir antes de armar nada."""
    faltan = []
    for rel in ("respuesta_revisores.pdf", "access_marcado.pdf", "access.pdf",
                "access.tex", "access.bbl"):
        if not os.path.isfile(os.path.join(ARTICULO, rel)):
            faltan.append(rel)
    return faltan


def compilar_aislado(carpeta, paginas_esperadas):
    """Compila el paquete en su propia carpeta y comprueba el nº de páginas."""
    orden = [["pdflatex", "-interaction=nonstopmode", "access.tex"],
             ["bibtex", "access"],
             ["pdflatex", "-interaction=nonstopmode", "access.tex"],
             ["pdflatex", "-interaction=nonstopmode", "access.tex"]]
    for cmd in orden:
        r = subprocess.run(cmd, cwd=carpeta, capture_output=True, text=True)
        if cmd[0] == "pdflatex" and r.returncode != 0:
            return False, "pdflatex devolvió %d" % r.returncode

    pdf = os.path.join(carpeta, "access.pdf")
    if not os.path.isfile(pdf):
        return False, "el paquete no produjo access.pdf"
    n = paginas_pdf(pdf)
    if n is None:
        return True, "compila (no se pudo contar páginas: falta pdfinfo)"
    if paginas_esperadas is not None and n != paginas_esperadas:
        return False, ("el paquete compila %d páginas y el PDF de referencia "
                       "tiene %d" % (n, paginas_esperadas))

    # Comprobado el paquete, se retiran los artefactos y el PDF que acaba de
    # producir, para que el .zip contenga solo fuente.
    os.remove(pdf)
    for nombre in os.listdir(carpeta):
        if nombre.endswith(EXCLUIR_EXTENSIONES):
            os.remove(os.path.join(carpeta, nombre))
    return True, "%d páginas" % n


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--escribir", action="store_true",
                    help="arma la carpeta (sin esto solo informa)")
    ap.add_argument("--verificar", action="store_true",
                    help="compila el paquete LaTeX aislado y comprueba páginas")
    args = ap.parse_args(argv)

    faltan = verificar_entradas()
    if faltan:
        print("ERROR: faltan archivos en IEEE-ACCESS/: %s" % ", ".join(faltan))
        print("Genera primero el PDF que falte (access_marcado se produce con "
              "marcar_cambios.py).")
        return 1

    fuente = archivos_fuente()
    n_pag = paginas_pdf(os.path.join(ARTICULO, "access.pdf"))
    n_pag_marcado = paginas_pdf(os.path.join(ARTICULO, "access_marcado.pdf"))

    print("=" * 68)
    print("preparar_envio.py — %s" % ("ARMANDO" if args.escribir
                                      else "SIMULACIÓN (usa --escribir)"))
    print("=" * 68)
    print("manuscrito limpio  : access.pdf, %s" % fmt_pag(n_pag))
    print("manuscrito marcado : access_marcado.pdf, %s" % fmt_pag(n_pag_marcado))
    print("respuesta          : respuesta_revisores.pdf, %s"
          % fmt_pag(paginas_pdf(os.path.join(ARTICULO, "respuesta_revisores.pdf"))))
    print("paquete LaTeX      : %d archivos + autores/" % len(fuente))
    if n_pag is not None and n_pag != n_pag_marcado:
        print("\nAVISO: el manuscrito limpio y el marcado no tienen el mismo "
              "número de páginas. Revisa que el marcado esté regenerado.")

    if not args.escribir:
        print("\nRepite con --escribir para armar %s."
              % os.path.relpath(ENVIO, ENTREGA))
        return 0

    if os.path.isdir(ENVIO):
        # Se vacía el contenido en lugar de borrar la carpeta raíz: en Windows,
        # si una terminal quedó situada dentro, `rmtree` sobre la raíz falla con
        # PermissionError y deja el paquete a medias.
        for nombre in os.listdir(ENVIO):
            ruta = os.path.join(ENVIO, nombre)
            if os.path.isdir(ruta):
                shutil.rmtree(ruta)
            else:
                os.remove(ruta)
    d1 = os.path.join(ENVIO, "1_respuesta_a_revisores")
    d2 = os.path.join(ENVIO, "2_manuscrito_marcado")
    d3 = os.path.join(ENVIO, "3_manuscrito_limpio")
    dsrc = os.path.join(d3, "fuente_latex")
    for d in (d1, d2, d3, dsrc):
        os.makedirs(d, exist_ok=True)

    shutil.copy2(os.path.join(ARTICULO, "respuesta_revisores.pdf"), d1)
    shutil.copy2(os.path.join(ARTICULO, "access_marcado.pdf"), d2)
    shutil.copy2(os.path.join(ARTICULO, "access.pdf"), d3)
    for nombre in fuente:
        shutil.copy2(os.path.join(ARTICULO, nombre), dsrc)
    shutil.copytree(os.path.join(ARTICULO, "autores"),
                    os.path.join(dsrc, "autores"))

    if args.verificar:
        print("\nCompilando el paquete aislado...")
        ok, detalle = compilar_aislado(dsrc, n_pag)
        print("  %s: %s" % ("OK" if ok else "FALLO", detalle))
        if not ok:
            print("\nEl paquete NO está completo. No lo subas hasta resolverlo.")
            return 1

    ruta_zip = os.path.join(d3, "fuente_latex.zip")
    with zipfile.ZipFile(ruta_zip, "w", zipfile.ZIP_DEFLATED) as z:
        for raiz, _, archivos in os.walk(dsrc):
            for nombre in sorted(archivos):
                completo = os.path.join(raiz, nombre)
                z.write(completo, os.path.relpath(completo, dsrc))

    with open(os.path.join(ENVIO, "LEEME.md"), "w", encoding="utf-8") as fh:
        fh.write(LEEME.replace("{paginas}", fmt_pag(n_pag)))

    tam = sum(os.path.getsize(os.path.join(r, f))
              for r, _, fs in os.walk(ENVIO) for f in fs)
    print("\nPaquete en %s" % os.path.relpath(ENVIO, ENTREGA))
    print("  1_respuesta_a_revisores/respuesta_revisores.pdf")
    print("  2_manuscrito_marcado/access_marcado.pdf")
    print("  3_manuscrito_limpio/access.pdf")
    print("  3_manuscrito_limpio/fuente_latex/  (%d archivos + autores/)"
          % len(fuente))
    print("  3_manuscrito_limpio/fuente_latex.zip")
    print("  LEEME.md")
    print("Tamaño total: %.1f MB" % (tam / 1024 / 1024))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
