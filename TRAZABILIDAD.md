# Traceability — from result to code

A map of every figure/table in the manuscript to its **generating script** and
its **data**. It allows any number in the article to be reproduced from source.

Paths are relative to this repo: code in `codigo_simulacion/`, data in
`datos/`, LTspice validation in `validacion_ltspice/`, manuscript in
`manuscrito/`.

## How this repo is produced

This repository is **generated, not edited**. The single source of truth is the
authoring folder `recomendacion01_canal_fisico/`, where the code is written and
every figure is executed. A release script mirrors it here and into the
manuscript folder, checking every file by SHA-256:

```
recomendacion01_canal_fisico/          <- edit and run here, and only here
    hacer_release.py  --aplicar   ──┬──►  equ-max-band_high-speed/   (this repo)
                                    └──►  IEEE-ACCESS/               (manuscript)
```

Two rules follow from this, and both matter for reproducibility:

- **No figure is retouched after `savefig`.** Every article figure is, byte for
  byte, the output of the script listed below, saved at 300 dpi with
  `bbox_inches="tight"`.
- **No file is copied by hand.** Manual copying is what previously allowed three
  different generations of the same figure to coexist.

## Manuscript figures

| Figure (in the `.tex`) | Generator | Data output |
| :--- | :--- | :--- |
| `fig_bode.png` (ladder $N$ vs $1/\cosh$) | `notebooks/simulacion_canal_rcEjecutadoActual.ipynb` | — |
| `fig_atenuacion_dbporpulgada.png` | `generar_comparativa.py` | `datos/tablas/tabla_comparativa.csv` |
| `fig_respuesta_pulso.png` | `generar_comparativa.py` | `tabla_comparativa.csv` |
| `fig_barrido_tasa_equal_il.png` / `_fixed.png` | `barrido_tasa_fisico.py` (arg. `criterio`) | `tabla_tasa_viable.csv` |
| `fig_sensibilidad_canal.png` | `sensibilidad_canal.py` | `tabla_sensibilidad_canal.csv` |
| `fig_ojo_escenarios.png` | `simulacion_ojo.py` | `tabla_ojo_ber.csv` |
| `fig_redoe_fixed.png` | `doe_canal_fisico.py` (arg. `criterio`) | `tabla_redoe.csv` (+ `doe_resultados.json`) |
| `fig_coopt_fixed.png` | `coopt_ctle_dfe.py` (arg. `criterio`) | `tabla_coopt.csv` |
| `fig_banera_escenarios.png` | `barrido_ojo_horizontal.py` | `tabla_ojo_horizontal.csv` |
| `fig_estres_jitter.png` | `estres_jitter.py` | `tabla_estres_jitter.csv` |
| `fig_pam8_ber.png` | `notebooks/simulacion_canal_rcEjecutadoActual.ipynb` | — |
| `fig_throughput.png` | `notebooks/simulacion_canal_rcEjecutadoActual.ipynb` | — |
| `fig_validacion_ber_mc.png` | `validacion_ber_mc.py` | `tabla_validacion_ber_mc.csv` |
| **`fig_ltspice_ojos.png`** (LTspice↔Python cross-check) | `validacion_ltspice/figura_articulo_ojos.py` | `validacion_ltspice/ref_python.npz` + `cadena_completa.raw` |

## Validations cited in the text (no inserted figure)

| Result in text | Generator | Output |
| :--- | :--- | :--- |
| Touchstone cross-check (passivity, $S_{21}$) | `exportar_touchstone.py` | `datos/touchstone/*.s2p`, `fig_touchstone_check.png` |
| $N_{sym}$ convergence (<0.3 %) | `convergencia_nsym.py` | `tabla_convergencia_nsym.csv` |
| Threshold sensitivity / extended sweep (threshold table) | `barrido_extendido.py` | `tabla_extendido_umbral.csv` |
| 2D co-optimization / jitter | `coopt_2d_jitter.py` | `tabla_coopt2d.csv` |

## Manuscript tables

| Table | Source |
| :--- | :--- |
| Viable rate (rate-reach) | `barrido_tasa_fisico.py` + `barrido_extendido.py` (saturations) |
| Worst-eye opening (4 GBaud) | `simulacion_ojo.py` → `tabla_ojo_ber.csv` |
| BER per scenario (4 GBaud) | `simulacion_ojo.py` |
| Viable rate at different thresholds | `barrido_extendido.py` → `tabla_extendido_umbral.csv` |
| Re-DoE (EQ-RC vs re-opt.) | `doe_canal_fisico.py` → `tabla_redoe.csv` |
| Joint co-optimization contribution | `coopt_ctle_dfe.py` → `tabla_coopt.csv` |
| Master parameter table | `manuscrito/tabla_parametros.tex` (static) |

## LTspice cross-validation (where it is)

- **In the article:** Section *Results*, subsection *"Model validation and
  physical-channel characterization"*, in the paragraph after Fig.
  `fig:senscanal`, followed by Fig. `fig:ltspice` (`fig_ltspice_ojos.png`).
  In `manuscrito/access.tex`, the subsection labelled `\label{sec:results}`
  (the figure carries `\label{fig:ltspice}`).
- **Reproducible kit:** `validacion_ltspice/` — see `README_LTSPICE.md`.
  Results: RC step 0.0019; RC/CTLE Bode 0.000 dB; dispersive `.ac` vs
  `.s2p` 0.55 dB; RC chain LTspice↔Python RMS 0.09 %; eye Δopening 0.011 V.

## Revision 1 — response to the reviewers

Material added to answer the IEEE Access review of manuscript
Access-2026-32692. Naming is `rNcM_`, where `N` is the reviewer number and `M`
the comment number, so that each file states which concern it answers. Code
lives in `codigo_simulacion/revision01/`.

| Concern | What the reviewer asked | Generator | Output |
| :--- | :--- | :--- | :--- |
| R1.1 | Realizing the dispersive channel from RC parameters; multipath fading | `revision01/r1c1_reflexiones_discontinuidad.py` | `fig_reflexiones.png` (`fig:reflexiones`), `revision01/tablas/tabla_r1c1_reflexiones.csv`; manuscript Sections II-B and III-B (`sec:reflexiones`) |
| R1.2 | Relation of FFE/CTLE to the implementation; how CSI is generated | text only | manuscript Section II-E (chain implementation and stage costs) and II-B (`sec:canalfisico`, deterministic channel, no CSI estimation) |
| R1.3 | Correlation of hardware and software implementations | `revision01/r1c3_constelacion.py` | `fig_constelacion_motores.png` (`fig:constelacion`), `revision01/tablas/tabla_r1c3_constelacion.csv`; manuscript Section III-A |
| R1.4 | Symbol-detection mechanism at the receiver | `revision01/r1c4_detector_pdf_niveles.py` | `fig_detector_niveles.png` (`fig:detector`), `revision01/tablas/tabla_r1c4_detector.csv`; manuscript Section II-F (`sec:detector`), eqs. `eq:fase`–`eq:ber` |
| R2.1 | New engineering insight beyond the expected phenomenon | text only | manuscript Introduction, paragraph after the contributions list; draws on the results of R1.1, R2.2 and R2.4 |
| R2.2 | Comparison against existing co-design methodologies | `revision01/r2c2_baseline_comparativa.py` | `fig_baseline_criterios.png` (`fig:baseline`), `revision01/tablas/tabla_r2c2_baseline.csv`; manuscript Section III-I (`sec:baseline`) and Contributions paragraph |
| R2.3 | Simulation-only validation | text only | manuscript Discussion, final paragraph (experimental path and what it would settle) |
| R2.4 | Generality of the $\omega_z \propto$ bandwidth rule | `revision01/r2c4_generalizacion_arquitectura.py` | `fig_generalizacion_wz.png` (`fig:generalidad`), `revision01/tablas/tabla_r2c4_generalizacion.csv`; manuscript Section III-I (`sec:generalidad`), design rule (i), Table IV footnote |

Run log, decisions and numerical outcomes: `revision01/BITACORA_REVISION01.md`.

### Resubmission documents

| File (in the submission folder) | What it is | How it is produced |
| :--- | :--- | :--- |
| `respuesta_revisores.pdf` | Point-by-point response to the reviewers | Written by hand |
| `access_marcado.pdf` | Manuscript with the changes highlighted | **Generated** by `revision01/marcar_cambios.py` |
| `access.pdf` | Clean manuscript | Written by hand |

`marcar_cambios.py` locates each added or modified passage by a start and an end
string and wraps it in a shaded environment. If any anchor is no longer found it
stops without writing, so a partially highlighted copy cannot be produced by
accident. Regenerate it after every change to the manuscript.

## Environment

- **Python 3.14.4** verified (NumPy 2.4.4, SciPy 1.17.1, Matplotlib 3.10.8);
  3.11+ expected to work (`codigo_simulacion/requirements.txt`,
  `environment.yml`).
- **LTspice** ADI 26.0.2.1 (batch mode `-b`) for the cross-validation.
- On Windows, export `PYTHONIOENCODING=utf-8` before redirecting a script's
  output to a file: the console summaries print Greek symbols and `cp1252`
  aborts on them.
