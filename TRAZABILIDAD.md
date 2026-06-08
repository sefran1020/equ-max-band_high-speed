# Traceability — from result to code

A map of every figure/table in the manuscript to its **generating script** and
its **data**. It allows any number in the article to be reproduced from source.

Paths are relative to this repo: code in `codigo_simulacion/`, data in
`datos/`, LTspice validation in `validacion_ltspice/`, manuscript in
`manuscrito/`.

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
  In `manuscrito/elsarticle-IJEC.tex`, the subsection labelled
  `\label{sec:results}` (the figure carries `\label{fig:ltspice}`).
- **Reproducible kit:** `validacion_ltspice/` — see `README_LTSPICE.md`.
  Results: RC step 0.0019; RC/CTLE Bode 0.000 dB; dispersive `.ac` vs
  `.s2p` 0.55 dB; RC chain LTspice↔Python RMS 0.09 %; eye Δopening 0.011 V.

## Environment

- **Python 3.14**, `numpy` / `scipy` / `matplotlib`
  (`codigo_simulacion/requirements.txt`, `environment.yml`).
- **LTspice** ADI 26.0.2.1 (batch mode `-b`) for the cross-validation.
