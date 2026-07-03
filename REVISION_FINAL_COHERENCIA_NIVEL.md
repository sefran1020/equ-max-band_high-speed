# Revision final de consistencia, coherencia y nivel

Manuscrito revisado: `IEEE-ACCESS/access.tex`  
PDF revisado: `IEEE-ACCESS/access.pdf`

## Veredicto editorial

El manuscrito esta en un nivel tecnicamente defendible para presentacion a IEEE Access. La contribucion principal es clara y con potencial de impacto: cuantifica la perdida de transferencia de un co-diseno de ecualizacion optimizado en un canal RC idealizado hacia canales fisicos dispersivos, y muestra que la contabilidad correcta del ruido puede cambiar incluso el signo del beneficio de una etapa lineal.

La estructura esta bien alineada con un articulo de investigacion aplicado: motivacion, contribuciones, metodologia reproducible, validacion, resultados, discusion, limitaciones, disponibilidad de codigo/datos y biografias de autores.

## Correcciones aplicadas en esta revision

- Se corrigio el encabezado heredado para que coincida con el titulo actual del articulo.
- Se corrigio el comentario interno de biografias para que coincida con el bloque final de autores.
- Se verifico que el PDF ya muestra los seis autores incorporados.
- Se mantuvo el DOI como `to be assigned by IEEE Access`, evitando un DOI inventado.

## Consistencia interna

No se observan contradicciones graves entre resumen, contribuciones, resultados, discusion y conclusion. Las cifras principales se mantienen consistentes:

- mejora por re-optimizacion: `+72%` en FR-4 y `+96%` en Megtron 6 para FFE+CTLE;
- mejora de co-optimizacion conjunta frente a ajuste por etapas: `+12.4%`;
- umbral principal: `BER <= 10^-2`;
- ancho de banda RC: `0.97 GHz`;
- criterio `equal_il` separado del criterio `fixed`.

Las aclaraciones sobre validacion ya estan mejor acotadas: el articulo no afirma validacion experimental completa, sino chequeos complementarios contra respuesta analitica, Touchstone y LTspice para la cadena RC.

## Coherencia narrativa

El hilo argumental queda claro:

1. Dos atajos metodologicos sesgan estudios exploratorios: canal RC idealizado y ruido mal referido.
2. Se construye un banco reproducible con ruido referido a entrada y ancho de banda finito.
3. El co-diseno RC no transfiere directamente a canales fisicos.
4. El CTLE puede perjudicar canales benignos por amplificacion de ruido.
5. CTLE y DFE deben co-optimizarse con el objetivo de la cadena completa.

La discusion no contradice los resultados y la conclusion conserva el mensaje central sin sobreprometer alcance experimental.

## Nivel y riesgo ante revisores

Fortalezas:

- Contribucion metodologica concreta, no solo incremental.
- Reproducibilidad fuerte: repositorio, scripts, CSV, Touchstone y trazabilidad.
- Buen manejo de limitaciones: DFE asistido, BER semi-analitico, ausencia de vias/reflexiones/crosstalk y falta de medicion fisica.
- Comparacion PAM-4/PAM-8 y analisis de jitter funcionan como pruebas de alcance.

Riesgos:

- Un revisor de signal integrity podria pedir validacion con mediciones reales o S-parameters medidos. El manuscrito ya lo reconoce como trabajo futuro.
- Las biografias de algunos autores son amplias y no todas estan centradas en interconexiones de alta velocidad; esto no afecta el articulo, pero conviene que el sistema de envio incluya roles/contribuciones claros si la plataforma los solicita.
- La etiqueta de DOI como `to be assigned by IEEE Access` es adecuada para evitar un numero falso; si el sistema de IEEE exige otro formato, debe ajustarse en la carga final.

## Estado tecnico

La compilacion con `pdflatex -interaction=nonstopmode -halt-on-error access.tex` no reporta errores, citas indefinidas ni referencias indefinidas. Persisten avisos de fuente y maquetacion propios de la plantilla IEEE Access, ademas de algunas tablas ajustadas, pero no bloquean la generacion del PDF.

## Recomendacion

El manuscrito esta listo para preparacion de envio. Antes de cargarlo, confirmar con todos los autores:

- aprobacion del manuscrito final;
- orden y nombres de autores;
- afiliacion institucional;
- ausencia de conflicto de interes;
- que el manuscrito no esta publicado ni bajo revision en otra revista;
- licencia/repositorio y commit de disponibilidad de datos.

