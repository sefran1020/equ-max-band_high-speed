"""
marcar_cambios.py — Genera la versión marcada del manuscrito para el reenvío.

IEEE Access exige subir, junto al manuscrito limpio, una copia con todos los
cambios resaltados. Este script la produce a partir de `access.tex`, de modo que
la versión marcada no se mantiene a mano y no puede quedar desfasada del
manuscrito: se regenera cuando haga falta.

    python revision01/marcar_cambios.py            # informe de lo que marcaría
    python revision01/marcar_cambios.py --escribir # escribe access_marcado.tex

Mecanismo: se añaden al preámbulo los paquetes `framed` y `xcolor` y cada
pasaje nuevo o modificado se envuelve en un entorno `shaded` con fondo amarillo.
El resaltado se aplica al TEXTO; las figuras nuevas no se envuelven, porque un
flotante dentro de `shaded` se desplaza fuera de la caja. Para no dejarlas sin
declarar, la versión marcada lleva una nota al inicio que las enumera.

Cada pasaje se localiza por una cadena de inicio y una de fin, ambas únicas en
el archivo. Si alguna deja de encontrarse —porque el texto cambió— el script se
detiene en lugar de producir una versión marcada incompleta, que es el peor
resultado posible para un reenvío con una sola oportunidad.
"""

import argparse
import os
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))
ENTREGA = os.path.normpath(os.path.join(AQUI, "..", ".."))
ARTICULO = os.path.join(ENTREGA, "IEEE-ACCESS")
ORIGEN = os.path.join(ARTICULO, "access.tex")
DESTINO = os.path.join(ARTICULO, "access_marcado.tex")

PREAMBULO = r"""
%% --- Añadido por revision01/marcar_cambios.py: resaltado de cambios ------- %%
\usepackage{framed}
\usepackage{xcolor}
\definecolor{shadecolor}{rgb}{1.0,0.94,0.55}
%% ------------------------------------------------------------------------- %%
"""

NOTA = r"""
\begin{shaded}
\noindent\textbf{Note for the review of this revision.} The highlighted
passages are the text added or modified in response to the reviewers' comments.
Figures~2, 8, 9, 15 and 16, equations~(7)--(10), Section~II-F, Section~III-B,
Section~III-I and Section~III-J are new; the footnote of Table~7 and the
reliability item of Section~II-H were amended. Figures were regenerated from
the released code without post-processing, so their appearance differs from the
previous submission although the reported numbers are unchanged; this is
explained in the response document.
\end{shaded}
"""

# (etiqueta, cadena de inicio, cadena de fin). Ambas deben ser únicas.
PASAJES = [
    ("Contribuciones: apertura reescrita (R2.2)",
     "This work quantifies the gap between ideal distributed-RC equalizer tuning and",
     "against the most related literature. Its contributions are the following."),

    ("Introducción: hallazgos no evidentes (R2.1)",
     "That a tuning derived on an idealized channel should be re-derived on a",
     "which is where the noise model earns its place."),

    ("II-B: canal guiado determinista, sin CSI (R1.1, R1.2)",
     "Two properties of this channel class should be made explicit, because they",
     "Their effect is quantified in\nSection~\\ref{sec:results}."),

    ("II-E: implementación de la cadena y coste de cada etapa (R1.2)",
     "\\emph{How the three stages act on the ISI, and how they are implemented.} The",
     "which is what Sections~\\ref{sec:results} onward measure."),

    ("II-F: subsección nueva de detección de símbolo (R1.4)",
     "\\subsection{Symbol detection at the receiver}",
     "$8.0\\times10^{-6}$ to $6.8\\times10^{-7}$ at the operating point."),

    ("II-H: indicador de fiabilidad, ahora referido a II-F (R1.4)",
     "\\item \\textbf{Reliability}: \\emph{semi-analytical} BER, obtained from the",
     "Carlo counting of low BERs \\cite{6248821,7202856,10415105}. This estimate assumes Gaussian"),

    ("III-A: correlación entre motores de simulación (R1.3)",
     "This cross-check is carried to the quantities the results are actually built",
     "properties of the modeled link and not of the tool that renders them."),

    ("III-B: subsección nueva de discontinuidades (R1.1)",
     "\\subsection{Impedance discontinuities: the reflective analogue of multipath}",
     "not judged from the return loss alone."),

    ("Tabla 7: nota al pie sobre el óptimo no resuelto (R2.4)",
     "\\footnotesize{$^{\\ast}$ Saturates at the maximum of the sweep (32~GBaud);\nlower bound.",
     "Section~\\ref{sec:generalidad}.}"),

    ("Regla de diseño: ancho de banda a -10 dB (R2.4)",
     "\\textit{Design rule.} The combined evidence yields a concrete, transferable",
     "tail cancellation to the DFE."),

    ("III-I: subsección nueva de comparación con criterios estándar (R2.2)",
     "\\subsection{Comparison against standard tuning criteria}",
     "addresses."),

    ("III-J: subsección nueva de generalidad de la regla (R2.4)",
     "\\subsection{Generality of the design rule}",
     "begins to droop."),

    ("Discusión: reflexiones acotadas en lugar de fuera de alcance (R1.1)",
     "and length). Reflections are treated as a bounded case rather than a full",
     "\\cite{9134991,5089870}."),

    ("Discusión: ruta a la validación experimental (R2.3)",
     "The most consequential of these limitations is that the study is entirely",
     "which is a weaker claim than verifying that it matches a physical\nboard."),
]


def marcar(texto):
    """Envuelve cada pasaje en `shaded`, verificando que todos existan."""
    fallos = []
    for etiqueta, inicio, fin in PASAJES:
        if texto.count(inicio) != 1:
            fallos.append("%s: la cadena de INICIO aparece %d veces"
                          % (etiqueta, texto.count(inicio)))
            continue
        i = texto.index(inicio)
        j = texto.find(fin, i)
        if j < 0:
            fallos.append("%s: no se encontró la cadena de FIN tras el inicio"
                          % etiqueta)
            continue
        j += len(fin)
        texto = (texto[:i] + "\\begin{shaded}\n" + texto[i:j]
                 + "\n\\end{shaded}" + texto[j:])
    return texto, fallos


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--escribir", action="store_true",
                    help="escribe access_marcado.tex (sin esto solo informa)")
    args = ap.parse_args(argv)

    with open(ORIGEN, encoding="utf-8") as fh:
        texto = fh.read()

    marcado, fallos = marcar(texto)

    print("=" * 68)
    print("marcar_cambios.py — %d pasajes declarados" % len(PASAJES))
    print("=" * 68)
    for etiqueta, _, _ in PASAJES:
        print("  - %s" % etiqueta)
    if fallos:
        print("\nERRORES:")
        for f in fallos:
            print("  ! %s" % f)
        print("\nNo se escribe nada: una versión marcada incompleta es peor que "
              "ninguna.")
        return 1

    # Preámbulo y nota inicial.
    ancla_preambulo = "\\usepackage{url}"
    if ancla_preambulo not in marcado:
        print("ERROR: no se encontró el ancla del preámbulo.")
        return 1
    marcado = marcado.replace(ancla_preambulo,
                              ancla_preambulo + "\n" + PREAMBULO, 1)

    ancla_nota = "\\maketitle"
    if ancla_nota not in marcado:
        print("ERROR: no se encontró \\maketitle.")
        return 1
    marcado = marcado.replace(ancla_nota, ancla_nota + "\n" + NOTA, 1)

    print("\nTodos los pasajes localizados.")
    if not args.escribir:
        print("Repite con --escribir para generar %s."
              % os.path.relpath(DESTINO, ENTREGA))
        return 0

    with open(DESTINO, "w", encoding="utf-8") as fh:
        fh.write(marcado)
    print("Escrito: %s" % os.path.relpath(DESTINO, ENTREGA))
    print("Compilar con: pdflatex -> bibtex -> pdflatex -> pdflatex")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
