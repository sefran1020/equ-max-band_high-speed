# Bitácora — Revisión 01 (IEEE Access, Access-2026-32692)

Registro de corridas y decisiones de la respuesta a los revisores. Mismo formato
que `../bitacora_resultados.md`: qué se hizo, por qué, qué salió y cómo
reproducirlo.

---

## 09-ago-2026 — Reorganización previa de la carpeta

Antes de generar evidencia nueva se saneó la estructura, porque había cuatro
árboles con el mismo material y las figuras del artículo se habían separado del
código que dice generarlas.

**Hallazgo.** Existían tres copias íntegras del repo de reproducibilidad
(`equ-max-band_high-speed/`, `equ-max-band_high-speed-upload/`, `repoArticulo/`),
verificadas idénticas con `diff -rq`. Y la misma figura existía en tres
generaciones distintas:

| Generación | Ubicación | Fecha | Etiquetas |
| :--- | :--- | :--- | :--- |
| 1 | `recomendacion01_canal_fisico/figuras/` | 05-jun | español |
| 2 | `equ-max-band_high-speed/codigo_simulacion/figuras/` | 07-jun | inglés |
| 3 | `IEEE-ACCESS/` (la publicada) | 09-jun | inglés + recorte de márgenes |

Es decir: **la figura del artículo no era la que producía el script del repo**,
lo que contradice la afirmación del manuscrito de que toda figura es reproducible
desde una sola receta.

**Diagnóstico del retoque del 09-jun.** Comparando las cabeceras PNG, la tercera
generación solo difiere en unas decenas de píxeles de margen a igual DPI (300):
`fig_redoe_fixed` 2550×1590 → 2520×1548, `fig_sensibilidad_canal` 4020×1221 →
3984×1175, `fig_bode` 2180×1572 → 2175×1569. Es un **recorte de espacio en
blanco**, no un cambio de contenido. Como los scripts ya guardan con
`bbox_inches="tight"`, ese recorte era redundante.

**Decisiones.**

1. Fuente única de verdad: `recomendacion01_canal_fisico/`. Todo lo demás es
   derivado y se regenera.
2. Se eliminaron las dos copias redundantes del repo.
3. Se unificó el código en una sola versión, con docstrings en español y
   **etiquetas de gráfica en inglés** (las que van al artículo). Se verificó que
   la diferencia entre los dos árboles era exclusivamente de etiquetas, salvo
   dos mejoras que se adoptaron: `ax.grid(alpha=0.3)` en `sensibilidad_canal.py`
   y `scipy`/`pyDOE` descomentados en `environment.yml` y `requirements.txt`.
   Excepción en sentido inverso: `validacion_ltspice/figura_articulo_ojos.py` era
   más reciente en la carpeta canónica (paneles `(a)/(b)/(c)` en inglés) que en
   el repo, que aún apuntaba a una ruta muerta (`IEEE-conference-template-062824Nvo`).
4. Se incorporaron los notebooks al canónico (`notebooks/`), de donde salen
   `fig_bode`, `fig_pam8_ber` y `fig_throughput`.
5. Se añadió `../hacer_release.py`: sincroniza código, figuras y datos hacia el
   repo y hacia `IEEE-ACCESS/` verificando por hash SHA-256. Copiar a mano queda
   prohibido; la deriva de tres generaciones ya no puede repetirse.
6. **Sin retoque posterior**: la figura publicada es, byte a byte, la salida del
   script. Se descarta el recorte manual del 09-jun.

**Regeneración.** Se re-ejecutaron todos los scripts de figuras con el código
unificado. Entorno: Python 3.14.4, NumPy 2.4.4, SciPy 1.17.1, Matplotlib 3.10.8
(coincide con lo declarado en el artículo).

**Verificación de que la regeneración no altera ningún resultado.** Todos los
scripts terminaron con código 0 y reprodujeron exactamente los valores
publicados; solo cambia el trazado (etiquetas en inglés, sin títulos internos):

- `generar_comparativa.py`: atenuación a Nyquist 0.322 dB/in (FR-4) y 0.200
  dB/in (Megtron 6), pérdida de inserción 7.38 dB, longitudes equivalentes
  22.89 y 36.92 pulgadas.
- `barrido_tasa_fisico.py` → Tabla I del artículo: 2.88 / 9.76 / 13.73 (RC
  N=5), 4.58 / 9.31 / 9.71 (FR-4 *equal_il*), 5.16 / 10.83 / 11.88 (Megtron
  *equal_il*), 12.75 / 14.06 / 27.19 (FR-4 *fixed*), 16.27 (Megtron *fixed*
  FFE+CTLE). Coincidencia total.
- `doe_canal_fisico.py` → Tabla IV: ω_z re-optimizado 2.35 / 2.65 / 5.35 /
  9.85 GHz y 14.1 → 24.2 GBaud (+72 %, FR-4 *fixed*), 16.3 → 32 GBaud (+96 %,
  Megtron *fixed*). Coincidencia total.

Tiempos: `validacion_ber_mc.py` 309 s (es la corrida Monte Carlo),
`doe_canal_fisico.py` 128 s, `barrido_tasa_fisico.py` 114 s; el resto por
debajo de 45 s.

**Excepción registrada.** `fig_bode.png`, `fig_pam8_ber.png` y
`fig_throughput.png` provienen del notebook (`notebooks/_repro_bode.py`,
`notebooks/_repro_pam_throughput.py`, este último con dependencia de `pyDOE`) y
no se re-ejecutaron en esta pasada: se publica la salida archivada del notebook,
que ya está en el canónico. Queda como tarea menor re-ejecutarlos para cerrar el
círculo también ahí.

> Nota operativa (Windows): los scripts imprimen caracteres griegos (π, α, β, γ)
> en el resumen de consola. Al redirigir la salida a un archivo hay que exportar
> `PYTHONIOENCODING=utf-8`, o `cp1252` aborta la corrida con `UnicodeEncodeError`.

**Pendiente de esta fase:** ninguno.

---

## Corridas de la revisión

_(Se registra aquí cada script `rNcM_*` conforme se ejecuta: fecha, comando,
tiempo, resultado numérico y figura/tabla producida.)_

### R1.1 — Reflexiones por discontinuidad · 09-ago-2026

    PYTHONIOENCODING=utf-8 python revision01/r1c1_reflexiones_discontinuidad.py

Traza FR-4 de 10 pulgadas (criterio `fixed`), discontinuidad a mitad de camino,
15 realizaciones, barrido hasta 32 GBaud. Tiempo: ~3 min.

El comentario mezcla dos cosas: una premisa equivocada (que el canal dispersivo
se realiza con parámetros RC) y un concepto trasladado del canal radio
(*multipath fading*). Se responde corrigiendo ambas y, sobre todo, cuantificando
el fenómeno que sí es análogo al multitrayecto en un medio guiado: las
reflexiones en discontinuidades de impedancia.

**Modelo.** Cascada de matrices ABCD terminada en 50 Ω, con dos discontinuidades
canónicas: muñón de vía en circuito abierto (Z_stub = Zc·coth(γL_s), resonancia
en cuarto de onda) y tramo de 0.5 in a 35 Ω. La vía se modela con el dieléctrico
completo —está enterrada, no ve el eps_eff de microstrip— instanciando el modelo
con eps_r' = 2·eps_r − 1, de modo que se reutiliza la matemática RLGC ya
validada en lugar de duplicarla.

**Validación del modelo nuevo.** Al retirar la discontinuidad, la cascada ABCD
reproduce el canal adaptado del artículo (H = e^{−γℓ}) con desviación máxima de
**0.0001 dB hasta 20 GHz**, y su tasa viable da 14.05 ± 0.21 y 27.27 ± 0.42
GBaud frente a los 14.06 y 27.19 publicados en la Tabla I. El script aborta si
esa desviación supera 0.05 dB.

**Resultados.** Coste en tasa viable respecto de la traza uniforme:

| Discontinuidad | f_res | FFE+CTLE | +DFE |
| :--- | ---: | ---: | ---: |
| Muñón 20 mil | 71.2 GHz | −0.8 % | −0.8 % |
| Muñón 40 mil | 35.6 GHz | −1.7 % | −3.8 % |
| Muñón 60 mil | 23.7 GHz | −1.8 % | −7.6 % |
| Muñón 90 mil | 15.8 GHz | −3.4 % | −18.2 % |
| Muñón 120 mil | 11.9 GHz | −10.3 % | −30.7 % |

**Interpretación (esto es lo que hay que llevar al artículo).** La asimetría
entre etapas no es casual: el castigo aparece cuando la frecuencia de Nyquist de
operación alcanza la resonancia. La cadena FFE+CTLE opera cerca de 14 GBaud
(Nyquist 7 GHz), muy por debajo del hundimiento de 15.8 GHz del muñón de 90 mil;
la cadena con DFE opera cerca de 27 GBaud (Nyquist 13.6 GHz) y cae justo sobre
él. De ahí un complemento a la regla de diseño: un cero espectral no lo levanta
el refuerzo lineal ni lo cancela la realimentación, así que hay que dimensionar
el muñón para que f_res = v_p/(4·L_s) quede holgadamente por encima del Nyquist
que el enlace ecualizado llegará a alcanzar — restricción que se vuelve activa
precisamente porque la ecualización extiende ese alcance.

**El caso contraintuitivo, ya explicado.** El tramo de 35 Ω da FFE+CTLE
14.93 ± 0.11 (**+6.3 %** sobre la traza uniforme) y +DFE 25.87 ± 0.37 (−5.1 %).
Que una discontinuidad *mejore* el caso lineal exigía verificación antes de
publicarlo:

1. No es artefacto de interpolación: la curva BER-vs-tasa es monótona y el BER
   es genuinamente menor entre 11 y 15 GBaud (1×10⁻⁴ frente a 9×10⁻⁴ a
   11 GBaud).
2. No es más señal: la separación media entre niveles a 11 GBaud *baja* de
   2.413 a 2.332 V (−3.4 %), coherente con que |S21| es algo peor en toda la
   banda.
3. La causa está en σ: la desviación por nivel cae de ≈0.387 a ≈0.308 V
   (−21 %). Como el ruido se inyecta con la misma σ_in y atraviesa el mismo
   CTLE y el mismo LPF en ambos casos, esa caída **no** puede ser de ruido: es
   ISI residual. El eco reflejado cancela parcialmente la cola del canal a las
   tasas en las que llega con signo opuesto.

El beneficio es específico de la tasa y no sobrevive: con DFE, que opera cerca
de 26 GBaud, el mismo tramo cuesta −5.1 %. Con el mecanismo identificado, el
caso sí entra al manuscrito, con la lectura correcta: una discontinuidad hay que
evaluarla a la tasa a la que el enlace realmente opera, no juzgarla por las
pérdidas de retorno.

**En el manuscrito.** Párrafo nuevo en II-B que fija el marco (canal
determinista y guiado, sin *fading* ni estimación de CSI, y el modelo dispersivo
no se construye con parámetros RC) — cubre también parte de R1.2 — y nueva
subsección de resultados **III-B** *"Impedance discontinuities: the reflective
analogue of multipath"* con la figura `fig:reflexiones` a doble columna, situada tras el
contraste con LTspice y antes de los resultados de tasa viable. Se actualizó la
frase de limitaciones de la Discusión, que antes declaraba las reflexiones fuera
de alcance sin más. El documento compila limpio en **16 páginas**.

> Cuidado al insertar: el primer intento colocó la subsección justo después de
> la figura de sensibilidad, lo que dejó el párrafo de LTspice huérfano dentro de
> la sección nueva. Conviene comprobar en el `.aux` que la numeración de
> subsecciones sea la esperada después de cada inserción.

Entregables: `figuras/fig_reflexiones.png`, `tablas/tabla_r1c1_reflexiones.csv`.

---

### Documentos del reenvío · 09-ago-2026

**Respuesta a los revisores** (`IEEE-ACCESS/respuesta_revisores.tex`, 5 páginas,
documento independiente que no usa `ieeeaccess.cls`). Estructura pedida por la
carta de decisión: para cada comentario, la preocupación citada, la respuesta y
la acción, con punteros exactos a secciones, ecuaciones, figuras y tablas de la
versión revisada.

Dos decisiones de redacción que conviene no deshacer:

1. Donde los comentarios del Revisor 1 usan vocabulario de canal radio, la
   respuesta **no dice que el revisor se equivoque**: atribuye el malentendido a
   que el manuscrito no dejó explícita la clase de canal, que es además lo que
   ocurrió. El mismo criterio en R1.3: el reparo sobre el hardware se responde
   reconociendo que el manuscrito no decía con claridad que no lo hay.
2. Las dos divulgaciones incómodas —que el ω_z de Megtron *fixed* de la Tabla 7
   no es un óptimo resuelto, y que el procedimiento de optimización no aventaja
   a MMSE— **se declaran enteras**, porque ambas salieron de correr justo lo que
   los revisores pidieron y ambas quedan visibles en los CSV del repositorio
   público. Ocultarlas mientras se publica el repositorio sería detectable.

**Corrección de la primera versión de la carta (misma fecha).** El borrador
inicial ponía las dos divulgaciones en la apertura, como bloque destacado. Se
movieron al punto que las motiva —la #1 a R2.4, la #2 a R2.2— porque en la
apertura un editor que lea por encima se lleva "los autores admiten dos
problemas" en vez de "los autores respondieron ocho comentarios". La honestidad
no depende de la colocación: la apertura sigue señalando que hay dos elementos
matizados y dónde están. Se hicieron a la vez otros dos ajustes:

- Se bajó el tono de la #2. Decía *"offers no advantage over standard tuning
  criteria"*, más duro de lo que dicen los números (MMSE queda +0.7 % y +2.1 %,
  o sea empate práctico más barato) y sin recoger que el forzado a cero pierde
  20 % **porque** ignora el ruido, lo que respalda el argumento del artículo.
- Se añadió a la #1 lo que faltaba: que el punto afectado no interviene en
  ninguna conclusión, y que la recuperación del +72 % en FR-4 —la medición
  central— no se ve tocada, porque sus óptimos caen dentro de la rejilla.

También lleva una sección final de cambios que **nadie pidió** —figuras
regeneradas desde el código, procedimiento de publicación, y el crecimiento de
14 a 18 páginas— porque las figuras cambian de aspecto respecto de la versión
enviada y es mejor explicarlo que dejar que se note.

**Versión marcada** (`IEEE-ACCESS/access_marcado.tex/.pdf`, 18 páginas).
Generada por `revision01/marcar_cambios.py`, no editada a mano: declara catorce
pasajes por cadena de inicio y de fin, los envuelve en `shaded` con fondo
amarillo, y **aborta sin escribir** si alguno deja de encontrarse. Verificado
visualmente que el resaltado funciona en texto corrido, en ecuaciones, a través
del salto de columna y dentro del flotante de la Tabla 7. Las figuras nuevas no
se envuelven —un flotante dentro de `shaded` se sale de la caja— y en su lugar
la versión marcada abre con una nota que enumera figuras, ecuaciones y
secciones nuevas.

> Al tocar el manuscrito hay que regenerar la versión marcada. Si el script
> aborta, es que un pasaje declarado cambió de texto: hay que actualizar su
> ancla, no quitar el pasaje.

---

### Paquete de envío y limpieza del repositorio · 10-ago-2026

**Endurecimiento de la respuesta a R1.1 (y mención al editor).** Se rehízo la
respuesta al primer comentario del Revisor 1. La versión anterior atribuía todo
el malentendido a una omisión nuestra, y eso concedía de más: el artículo nunca
ofreció un canal con multitrayecto, ni en el primer envío. La clase de canal
—interconexión guiada en cobre sobre PCB— está en el título, en el resumen y en
la Sección II desde el principio. La respuesta ahora lo dice sin rodeos, con
respeto, y pide que la revisión se evalúe en los términos wireline en que está
planteada; a la vez mantiene las dos concesiones que sí corresponden (no
habíamos declarado explícitamente la clase de canal, y las reflexiones estaban
fuera de alcance sin cuantificar). El mismo argumento se añadió a la carta al
editor, como párrafo propio en la apertura.

También se ajustó R1.3 en la misma línea: no hay implementación en hardware
**y nunca se afirmó que la hubiera**; el artículo se presenta como estudio de
simulación y su sección de disponibilidad de código publica un paquete de
simulación, no datos de medición. Se mantiene la concesión de que debió decirse
en una frase.

**Paquete de envío.** Nuevo `revision01/preparar_envio.py`, que arma
`ENVIO_ACCESS_2026-32692/` con los tres archivos que pide la carta, cada uno en
su carpeta numerada, más un LÉEME que dice qué sube a qué campo del portal.

Lo que aporta de verdad es `--verificar`: **compila el paquete LaTeX aislado en
su propia carpeta** y comprueba que da las mismas 18 páginas que el PDF de
referencia. Si al paquete le faltara una figura, una tipografía o un archivo de
clase, se descubre ahí y no en el portal. Pasó a la primera.

Detalles resueltos por el camino:

- El PDF compilado se excluye del `.zip` de fuente: iba dentro y duplicaba
  13 MB en una subida con límite de tamaño. El paquete quedó en 52 MB, con el
  `.zip` en 14 MB y sin ningún artefacto de compilación dentro.
- `shutil.rmtree` sobre la carpeta raíz del paquete falla en Windows con
  `PermissionError` si una terminal quedó situada dentro. El script vacía el
  contenido en lugar de borrar la raíz.
- Contar páginas buscando `/Type /Page` en el PDF crudo da **cero** con estos
  archivos, porque los objetos van comprimidos. Peor que no medir, porque
  parece un dato. Se usa `pdfinfo`, con `None` y aviso si no está disponible.

**Repositorio.** Sin artefactos de compilación, sin `__pycache__` ni `.raw`;
47 MB. Se corrigió el README, que aún llamaba `repoArticulo/` a la raíz —nombre
de una de las copias que se eliminaron el 09-ago— y no mencionaba `revision01/`.
Ahora declara además que el repositorio corresponde a la versión **revisada** y
que ninguna figura se retoca tras guardarse.

Quedan 39 archivos modificados y 13 sin seguimiento respecto del último commit
(`9296a92`). **No se ha hecho commit**: el mensaje y el momento los decide el
autor.

---

### R1.2, R2.1 y R2.3 — Solo redacción · 09-ago-2026

Las tres se responden sin generar evidencia nueva: la que hacía falta ya la
produjeron los cinco scripts anteriores.

**R1.2 — relación de FFE/CTLE con la implementación, y generación del CSI.**
Se añadió a **II-E** un párrafo de implementación explícita: orden de bloques
(TX-FFE → canal → inyección de ruido → RX-CTLE → filtro de recepción → muestreo
→ DFE), FFE como FIR en el dominio de símbolo, canal y CTLE como multiplicación
por su transferencia sobre la secuencia sobremuestreada a 32 muestras por
símbolo —equivalente a una simulación temporal en estado nulo de los mismos
bloques LTI—, y el mapeo de las dos variables del DoE a los bloques (ω_z al
cero del CTLE, c_post al tap post-cursor del FFE). Se añadió después un párrafo
que explica por qué las tres etapas atacan la misma ISI con costes distintos:
el FFE pre-distorsiona en el transmisor y no amplifica ruido pero gasta
amplitud; el CTLE refuerza en el receptor, después del ruido, y por tanto lo
amplifica; el DFE resta post-cursores ya decididos, sin realce de ruido, pero
no puede con el pre-cursor. La parte del CSI ya estaba cubierta por el párrafo
que se añadió a **II-B** al responder R1.1 (canal determinista, conocido por
construcción, sin estimación); ahora se le puso etiqueta y el párrafo de II-E
la referencia.

**R2.1 — qué aporta más allá de confirmar lo esperado.** Párrafo nuevo tras la
lista de contribuciones, que concede primero lo obvio (que un ajuste hecho
sobre un canal idealizado haya que re-derivarlo es esperable) y enumera después
lo que no se sigue de ello: el beneficio de una etapa puede **cambiar de signo**
—en una traza de bajas pérdidas el refuerzo lineal reduce la tasa viable, o sea
que re-optimizar puede significar ecualizar menos—; el error del ajuste por
etapas tiene una **dirección definida** (−8.1 % frente a +12.4 % del conjunto);
el ancho de banda que gobierna el cero óptimo es el de −10 dB y **no** el corte
de −3 dB, que ordena dos canales al revés; la ecualización **crea su propia
restricción** al subir el Nyquist de operación hasta la resonancia del muñón
(3.4 % de coste en la cadena lineal, 18.2 % en la completa); y el procedimiento
de optimización **no es la palanca** —los criterios establecidos lo igualan con
1/13 del coste—, mientras que el criterio que ignora el ruido pierde el 20 %.
Todo con cifras ya medidas, sin adjetivos.

**R2.3 — validación solo simulada.** Párrafo final de la Discusión. Dice qué
haría falta (cupón con trazas de FR-4 y laminado de bajas pérdidas en las
longitudes modeladas, parámetros S con VNA hasta el doble del Nyquist más alto
examinado, fuente PAM-4 con conteo de errores en hardware, y un transceptor con
FFE/CTLE/DFE programables), qué zanjaría (las tasas absolutas, que dependen de
lo que se dejó fuera: encapsulados, conectores, diafonía y un DFE adaptando a
ciegas) y qué probablemente no cambiaría (las comparaciones relativas, por
calcularse entre configuraciones evaluadas en el mismo banco) — señalando que
esto último es un argumento, no una prueba. Y se dice explícitamente que la
evidencia simulada **no** se ofrece como equivalente a una medición: los tres
contrastes verifican que la cadena está implementada como se describe, que es
una afirmación más débil que verificar que coincide con una placa real.

**Estado del manuscrito:** 18 páginas, compila limpio, sin referencias
indefinidas ni cajas desbordadas.

---

### R2.2 — Comparación con los ajustes habituales · 09-ago-2026

    PYTHONIOENCODING=utf-8 python revision01/r2c2_baseline_comparativa.py

Tres canales dispersivos (FR-4 y Megtron 6 *equal_il*, FR-4 *fixed*), diez
realizaciones para la medición común. ~4 min. Se excluye Megtron *fixed* por la
saturación documentada en R2.4.

El reparo es justo: el artículo describe su procedimiento pero nunca lo midió
contra los criterios de ajuste habituales. Se separan las dos partes de un
procedimiento —el **criterio** que se optimiza y la **búsqueda** que lo
recorre— porque solo cabe atribuirse una de ellas, y ni siquiera esa.

Criterios comparados: el del artículo (tasa viable con ruido referido a la
entrada), MMSE con ganancia ajustada por mínimos cuadrados, forzado a cero del
primer post-cursor, y apertura del peor ojo. Búsquedas: LHS de 30 puntos (la
del artículo) y rejilla exhaustiva 13×13. Todo se mide después con la misma
vara.

**Resultado, relativo a la combinación del artículo (promedio de tres canales):**

| Búsqueda | Criterio | FFE+CTLE | +DFE | Evaluaciones |
| :--- | :--- | ---: | ---: | ---: |
| LHS-30 | tasa (artículo) | 100.0 % | 100.0 % | 390 |
| LHS-30 | MMSE | **100.7 %** | **102.1 %** | 30 |
| LHS-30 | forzado a cero | 79.6 % | 84.1 % | 30 |
| LHS-30 | apertura de ojo | 99.6 % | 99.3 % | 30 |
| rejilla 13×13 | tasa | 101.5 % | 101.6 % | 2197 |

**Lectura, sin adornos.** El procedimiento del artículo **no** aventaja a la
práctica establecida para ajustar la etapa lineal: MMSE alcanza la misma tasa
—algo mejor con DFE— con **1/13 del coste**, y la maximización de apertura de
ojo empata dentro de la dispersión de medida. La búsqueda exhaustiva gana 1.5 %
por 5.6× las evaluaciones, así que el muestreo LHS ya llega al entorno del
óptimo. Conclusión que hay que escribir tal cual: **la optimización no es una
aportación de este trabajo.**

Lo que la comparación sí respalda es el argumento del ruido: de los cuatro
criterios, el único que falla claramente (−20.4 % lineal, −15.9 % con DFE) es
el único que ignora el ruido del receptor. El forzado a cero anula el
post-cursor sin mirar el ruido que el refuerzo amplifica y por eso se va a un
cero altísimo (9.85 y 9.55 GHz en los dos casos de FR-4).

Conviene decir también dónde sí importa el criterio: en el ajuste conjunto de
etapa lineal y realimentación (Tabla V), un espacio de diseño que ninguno de
los tres criterios de referencia aborda.

**En el manuscrito.** Nueva subsección **III-I** *"Comparison against standard
tuning criteria"* con figura a una columna, y reformulación del párrafo de
contribuciones: se dice explícitamente que la aportación está en *qué* se mide y
*bajo qué modelo de receptor*, no en *cómo* se busca, y que "co-diseño" designa
el banco de evaluación y el protocolo de comparación, no un optimizador nuevo.
No se reintroduce el descargo de novedad que hundió la versión de AEU: precisar
el alcance de la aportación no es negarla. Sigue en **17 páginas**.

> Detalle de reproducción: la figura se regeneró dos veces. La primera corrida
> arrancó antes de que se tradujeran las etiquetas al inglés y publicó los
> nombres internos en español; conviene comprobar la figura después de editar
> el script, no después de lanzarlo.

Entregables: `figuras/fig_baseline_criterios.png`,
`tablas/tabla_r2c2_baseline.csv`.

---

### R2.4 — ¿La regla ω_z ∝ ancho de banda generaliza? · 09-ago-2026

    PYTHONIOENCODING=utf-8 python revision01/r2c4_generalizacion_arquitectura.py

Cinco canales (RC N=5, FR-4 y Megtron 6 en ambos criterios) × seis
configuraciones. ~5 min.

Configuraciones: **A0** la del artículo (1 cero, polos 15/30 GHz, A_dc 2.5);
**A1** polos 25/50 GHz; **A2** polos 8/16 GHz; **A3** dos ceros (ω_z y 2ω_z,
polos 15/30/30); **A4** A_dc = 4.0; **A5** arquitectura A0 pero el DoE maximiza
la **apertura de ojo** en lugar de la tasa viable — esa es la que responde a la
segunda mitad de la pregunta, sobre el montaje de optimización.

**Verificaciones.** (1) El evaluador parametrizado con A0 devuelve exactamente
el mismo BER que `evaluar_config()` del artículo (3.579549e-01 y 3.579136e-01,
idénticos). (2) A0 con el rango de búsqueda del artículo reproduce sus cuatro
ω_z publicados: 2.35 / 2.65 / 5.35 / 9.85 GHz. El script aborta si falla
cualquiera de las dos.

**Resultado del ajuste** ω_z* = k·BW^p:

| Config. | p | R² | n |
| :--- | ---: | ---: | ---: |
| A0 línea base | 0.89 | 0.988 | 4 |
| A1 polos anchos | 0.94 | 0.770 | 4 |
| A2 polos angostos | 0.58 | 0.904 | 5 |
| A3 dos ceros | 0.30 | 0.417 | 4 |
| A4 ganancia DC alta | 0.89 | 0.988 | 4 |
| A5 objetivo por ojo | 0.68 | 0.728 | 5 |

**Lectura.** Hay que separar dirección de magnitud. La **dirección** generaliza:
p > 0 en las seis, así que el cero óptimo siempre sube al ensanchar el canal.
La **magnitud** no: p va de 0.30 a 0.94. La regla se transfiere como
*dirección y procedimiento*, no como coeficiente — que es la misma conclusión
del artículo sobre el ajuste, aplicada un nivel más arriba: la regla de diseño
también hay que re-derivarla para otra arquitectura. Detalles: A4 ≡ A0
exactamente (la ganancia DC es irrelevante; el cero lo fija el perfil de
pérdidas del canal), y A3 debilita más la dependencia porque su segundo cero
aporta parte del refuerzo en alta frecuencia.

**Hallazgo 1: hay que decir a qué nivel de pérdida se mide el ancho de banda.**
El corte clásico de −3 dB **no** ordena los canales igual que sus ω_z: FR-4
*equal_il* tiene más banda a −3 dB que Megtron *equal_il* (0.59 vs 0.44 GHz) y
sin embargo un ω_z menor (2.35 vs 2.65). La frecuencia a −10 dB sí los ordena
(3.07 vs 3.52 GHz). Tiene sentido físico: el cero se coloca donde el ecualizador
tiene trabajo, no donde el canal empieza a caer. El artículo decía "ancho de
banda del canal" sin precisar; ahora lo precisa.

**Hallazgo 2: el ω_z publicado de Megtron 6 *fixed* no es un óptimo resuelto.**
En ese canal la tasa viable de FFE+CTLE se sale de la rejilla del DoE (2–26
GBaud), de modo que el objetivo **empata** entre muchos ω_z y el LHS se queda
con el mayor que muestreó. Se comprobó ampliando el rango de búsqueda: con
1–20, 1–40 y 1–60 GHz el óptimo sale 19.68, 39.35 y 59.02 GHz y el objetivo se
queda clavado en 26.00 GBaud — es decir, sigue el borde, no un máximo. Al
ampliar la rejilla de tasa a 2–56 GBaud el punto pasa a 11.45 GHz, pero el
objetivo vuelve a saturar. Se añadió detección de saturación a `re_doe()`: esos
puntos se marcan con `*` y **se excluyen del ajuste**. La Tabla IV del
manuscrito lleva ahora una nota al pie que lo declara; su nota anterior ya
advertía la saturación de la *tasa*, pero no la del ω_z.

> Trampa relacionada: la rejilla de tasa del DoE **forma parte del
> procedimiento publicado**. Al ampliarla para el estudio, la verificación
> contra la Tabla IV empezó a fallar (FR-4 *fixed* daba 4.45 en vez de 5.35).
> Hubo que parametrizar la rejilla: la del artículo para verificar, la ampliada
> para el estudio.

**En el manuscrito.** Se precisó el punto (i) de la regla de diseño ("ancho de
banda medido a −10 dB"), se añadió la subsección **III-I** *"Generality of the
design rule"* con la figura `fig:generalidad` a una columna, y se amplió la nota al pie de la
Tabla IV. El documento pasa a **17 páginas**.

Entregables: `figuras/fig_generalizacion_wz.png`,
`tablas/tabla_r2c4_generalizacion.csv`.

---

### R1.3 — Correlación entre implementaciones · 09-ago-2026

    PYTHONIOENCODING=utf-8 python revision01/r1c3_constelacion.py

Canal RC a 4 GBaud, 400 símbolos PAM-4. Segundos.

El comentario da por supuesto que hay una implementación en hardware. **No la
hay**, y hay que decirlo sin rodeos. Lo que sí existe es una segunda
implementación independiente en LTspice, que el manuscrito ya usaba pero solo
comparaba por apertura de ojo. Se lleva la comparación a las dos magnitudes que
el revisor nombra —constelación y BER— reutilizando los conjuntos
`ref_python.npz` y `cadena_completa.raw` que ya produjo el flujo de validación.

Clave metodológica: ambas formas de onda se procesan con el **mismo detector**
del artículo. Se importa `detectar()` de `r1c4_detector_pdf_niveles.py` en lugar
de reimplementarlo, de modo que la única diferencia entre las dos columnas de
resultados es el motor de simulación.

**Resultados.**

| Magnitud | Cadena en frecuencia | LTspice | Diferencia |
| :--- | ---: | ---: | ---: |
| μ nivel −3 / −1 / +1 / +3 (V) | −4.0119 / −1.3433 / +1.3354 / +4.0073 | −3.9957 / −1.3402 / +1.3295 / +3.9900 | ≤ 17 mV |
| σ por nivel (V) | 0.047 – 0.052 | 0.053 – 0.056 | ≤ 6 mV |
| Umbrales τ (V) | −2.6776 / −0.0040 / +2.6713 | −2.6679 / −0.0054 / +2.6597 | ≤ 12 mV |
| Apertura de ojo (V) | 2.3607 | 2.3309 | 30 mV (1.3 %) |
| BER | < 10⁻¹² | < 10⁻¹² | ambos bajo el piso |

Correlación muestra a muestra: **Pearson r = 0.999992**, RMS de diferencia
**17.0 mV = 0.41 %** del fondo de escala.

El BER de esta configuración cae en 10⁻¹⁴⁶ y 10⁻¹²⁷; se reportan como "bajo el
piso de 10⁻¹²" porque ahí la cola gaussiana deja de ser fiable y citar el número
sería falsa precisión. La comparación con sentido es la de la constelación.

**Trampa encontrada (importante para futuros usos).** `alinear_full()` solo
evalúa retardos que dejen sitio a su ventana de correlación (`n_score`, 400 por
defecto). Con una secuencia de exactamente 400 símbolos, el único retardo
admisible es 0 y **la alineación falla en silencio**: da μ ≈ 0 en los cuatro
niveles, σ ≈ 3 V y BER ≈ 0.38, sin lanzar ningún error. Se detectó comparando
contra una búsqueda por fuerza bruta de fase y retardo (que daba niveles
limpios en −4.19 / −1.33 / +1.42 / +4.17). Se añadió `n_score` como parámetro
de `detectar()` y se documentó en su docstring. En los barridos del artículo no
afecta, porque usan 1000 símbolos.

**En el manuscrito.** Párrafo nuevo en III-A tras la Fig. 5 de LTspice, que
declara explícitamente que no hay hardware y que el simulador de circuitos es
una implementación *de software* independiente, más la figura `fig:constelacion` a una sola
columna (constelación superpuesta + correlación). Se eligió una columna en vez
de tres paneles a doble columna para no gastar media página: el documento sigue
en **16 páginas**.

Entregables: `figuras/fig_constelacion_motores.png`,
`tablas/tabla_r1c3_constelacion.csv`.

---

### R1.4 — Mecanismo de detección de símbolo · 09-ago-2026

    PYTHONIOENCODING=utf-8 python revision01/r1c4_detector_pdf_niveles.py

Escenario: FR-4, criterio `equal_il`, 4 GBaud, realización 0 (semilla de ruido
1000). Tiempo: pocos segundos.

El comentario no señalaba un defecto del método sino de la exposición: el
detector estaba implementado en `cadena_enlace.evaluar_enlace()` pero el
manuscrito solo lo describía en prosa, sin ecuaciones. El script lo explicita en
tres etapas —instante de muestreo, realimentación de decisiones, umbrales y
BER— y produce la figura y la tabla de respaldo.

**Resultados.**

| Magnitud | FFE+CTLE | FFE+CTLE+DFE |
| :--- | ---: | ---: |
| Umbrales τ (V) | −2.083 / +0.058 / +2.117 | −2.078 / +0.062 / +2.116 |
| σ típica por nivel (V) | 0.220 – 0.253 | 0.196 – 0.226 |
| SER | 1.606×10⁻⁵ | 1.358×10⁻⁶ |
| BER = SER/2 | 8.030×10⁻⁶ | 6.791×10⁻⁷ |
| Apertura de ojo (V) | 0.528 | 0.712 |

Instante de muestreo: φ* = 20/32 = 0.625 UI con |ρ| = 0.9951, frente a |ρ| =
0.7231 en la peor fase — la elección del instante no es cosmética y conviene
mostrarla.

**Verificación de identidad.** El BER recalculado con las ecuaciones que se
publicarán coincide con el de `evaluar_enlace()` —la función que genera *todos*
los BER del artículo— hasta 12 cifras significativas: 8.030454×10⁻⁶ y
6.790812×10⁻⁷ en ambos caminos. El script aborta si dejan de coincidir, de modo
que la documentación no puede desincronizarse del código sin que se note.

Dos detalles a incorporar al manuscrito, ausentes hoy y que probablemente
motivaron el comentario:

1. Los umbrales **no** son los niveles nominales sino los puntos medios entre
   las medias recibidas, τ_k = (μ_k + μ_{k+1})/2: son adaptativos.
2. La conversión **BER = SER/2** descansa en el mapeo Gray de PAM-4 (un error
   entre niveles contiguos altera un solo bit de los dos). Nunca se enunció.

Entregables: `figuras/fig_detector_niveles.png`, `tablas/tabla_r1c4_detector.csv`.
(La figura se nombra sin el prefijo `r1c4` porque se publica en el artículo; la
trazabilidad la lleva el nombre del script.)

**En el manuscrito.** Nueva subsección II-F *"Symbol detection at the receiver"*
con las ecuaciones (7)–(10) —instante de muestreo, umbral adaptativo, `Pe` por
nivel y `BER = SER/2`— y la Fig. 2 a doble columna. El párrafo *Reliability* de
la sección de indicadores se acortó para que apunte a la nueva subsección en
lugar de repetir la descripción. El documento compila limpio y pasa de 14 a
**15 páginas**; hay que compensar con recortes antes del reenvío (candidatos:
la sección de PAM-8 y parte de la Discusión).
