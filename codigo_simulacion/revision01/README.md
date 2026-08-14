# Revisión 01 — Respuesta a los revisores de IEEE Access

Manuscrito **Access-2026-32692**, *"When Equalizer Co-Design Does Not Transfer:
From Idealized RC to Dispersive PCB Interconnects"*. Decisión del 08-ago-2026:
rechazo con **una única oportunidad de reenvío**. Observaciones originales en
`../../observaciones.txt`.

Esta carpeta contiene **todo lo que se genera para responder**, y nada más. La
regla de trazabilidad es que el nombre del archivo identifica la observación:

    rNcM_<tema>.py     N = número de revisor,  M = número de comentario

De modo que `r1c4_detector_pdf_niveles.py` responde al comentario 4 del
Revisor 1. Cualquiera —incluido el editor— puede ir de la observación al código
sin intermediarios.

## Mapa de observaciones

| Obs. | Qué pide el revisor | Tipo de respuesta | Entregable |
| :--- | :--- | :--- | :--- |
| **R1.1** | "Difícil realizar el canal dispersivo con parámetros RC; el canal práctico introduce *fading* con multitrayecto" | Corregir la premisa + evidencia nueva | `r1c1_reflexiones_discontinuidad.py` ✔ |
| **R1.2** | Relación de FFE/CTLE con la implementación; cómo se genera el CSI | Redacción (§ cadena de señal) | II-E + II-B ✔ |
| **R1.3** | Correlación entre implementación HW y SW para constelación y BER | Redacción + figura | `r1c3_constelacion.py` ✔ |
| **R1.4** | El mecanismo de detección de símbolo no está claro | Redacción formal + figura | `r1c4_detector_pdf_niveles.py` ✔ |
| **R2.1** | ¿Qué *insight* nuevo hay más allá de confirmar lo esperado? | Redacción (Intro) | tras la lista de contribuciones ✔ |
| **R2.2** | Comparar con metodologías de co-diseño existentes | Evidencia nueva | `r2c2_baseline_comparativa.py` ✔ |
| **R2.3** | Validación solo simulada | Redacción (§ ruta experimental) | final de la Discusión ✔ |
| **R2.4** | ¿La regla ω_z ∝ ancho de banda generaliza a otras arquitecturas? | Evidencia nueva | `r2c4_generalizacion_arquitectura.py` ✔ |

### Sobre R1.1 y R1.2 (marco de referencia del revisor)

Ambos comentarios usan vocabulario de **canal radio** (*multipath fading*, CSI,
*channel statistics*) sobre un enlace **guiado en cobre**. No son errores del
revisor sino un síntoma: el manuscrito no deja explícito que se trata de un
canal determinista y estacionario. La respuesta debe:

1. Aclarar que el canal dispersivo **no** se realiza con parámetros RC — ese es
   justamente el resultado del artículo. Se modela con RLGC dependiente de la
   frecuencia (efecto pelicular `R(f)` y pérdida dieléctrica `G(f)`, ec. 3).
2. Explicar que en un enlace guiado no hay *fading* por multitrayecto; el
   fenómeno análogo son las **reflexiones por discontinuidades de impedancia**
   (vías, *stubs*, conectores), que el manuscrito declara fuera de alcance.
   `r1c1_reflexiones_discontinuidad.py` lo cuantifica en lugar de solo declararlo.
3. Aclarar que **no hay estimación de CSI**: `H(f)` es conocida por construcción
   y la única aleatoriedad del experimento son el patrón PRBS y la semilla de
   ruido (semilla `1000+r` para la realización `r`, pareada entre canales).

## Los tres archivos del reenvío

IEEE Access pide subir tres cosas. Están en `../../IEEE-ACCESS/`:

| Archivo | Qué es | Cómo se produce |
| :--- | :--- | :--- |
| `respuesta_revisores.tex/.pdf` | Respuesta punto por punto (preocupación / respuesta / acción) | Se edita a mano |
| `access_marcado.tex/.pdf` | Manuscrito con los cambios resaltados en amarillo | **Generado** por `marcar_cambios.py` |
| `access.tex/.pdf` | Manuscrito limpio | Se edita a mano |

Y `preparar_envio.py` arma el paquete completo para el portal:

    python revision01/preparar_envio.py --escribir --verificar

El `--verificar` compila el paquete LaTeX **aislado en su propia carpeta** y
comprueba que da el mismo número de páginas que el PDF de referencia. Si al
paquete le faltara una figura o un archivo de clase, se descubre ahí y no en el
portal. El resultado queda en `../../ENVIO_ACCESS_2026-32692/`, con un LÉEME que
dice qué archivo va en cada campo del portal.

`marcar_cambios.py` localiza cada pasaje nuevo o modificado por una cadena de
inicio y otra de fin, y lo envuelve en un entorno `shaded`. Si alguna deja de
encontrarse —porque el texto cambió— **se detiene sin escribir nada**: una
versión marcada incompleta es peor que ninguna en un reenvío que admite una
sola oportunidad. Al tocar el manuscrito hay que regenerarla:

    python revision01/marcar_cambios.py             # informe
    python revision01/marcar_cambios.py --escribir  # escribe access_marcado.tex

Las figuras nuevas no se envuelven (un flotante dentro de `shaded` se sale de la
caja); en su lugar, la versión marcada lleva una nota inicial que las enumera.

## Convenciones

Idénticas a las de la carpeta canónica, para que el código se lea igual:

- Docstrings y comentarios **en español**; etiquetas de las gráficas **en inglés**
  (van al artículo).
- Figuras a `revision01/figuras/`, tablas CSV a `revision01/tablas/`.
- `fig.savefig(..., dpi=300, bbox_inches="tight")`, sin retoque posterior: la
  figura publicada debe ser byte a byte la salida del script.
- Los parámetros físicos se importan de `parametros.py`; no se duplican valores.

## Cómo se publica

Nada se copia a mano. Desde la carpeta canónica:

    python hacer_release.py              # informe de qué está desincronizado
    python hacer_release.py --aplicar    # copia al repo y a IEEE-ACCESS/

El contenido de `revision01/` viaja al repo como
`codigo_simulacion/revision01/`.

## Estado

Ver `BITACORA_REVISION01.md` para el registro de corridas y resultados.
