# Simulation Code — Physical Channel and Equalizer Co-Design

Python implementation of the physical-channel extension used in the article. It compares the distributed RC baseline with dispersive FR-4 and Megtron 6 transmission-line models including skin-effect and dielectric loss, and then evaluates equalizer co-design over those channels.

---

## Environment (isolated — does NOT touch the system Python)

Create a **dedicated** environment for this folder. Two options:

**venv (PowerShell, Windows):**
```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

**conda:**
```powershell
conda env create -f environment.yml
conda activate canal-fisico-rec01
```

Dependencies: `numpy`, `scipy`, `matplotlib`, and `pyDOE` for archived notebook checkpoints (see `requirements.txt`).

---

## Execution

```powershell
python generar_comparativa.py     # base deliverables (attenuation + pulse)
python simulacion_ojo.py          # eye + transmitted signal (FFE+CTLE+DFE chain)
python barrido_tasa_fisico.py     # viable BER rate vs GBaud (30 realizations, 2 criteria)
python doe_canal_fisico.py        # re-DoE: (w_z,c_post) re-optimized per substrate
python barrido_ojo_horizontal.py  # horizontal eye opening: bathtub + jitter
python coopt_ctle_dfe.py          # joint CTLE+DFE co-optimization
python coopt_2d_jitter.py         # 2D margin (vertical x horizontal) under jitter
```

`coopt_2d_jitter.py` integrates the vertical and horizontal margin analyses: it optimizes the EQ for the **2D viable rate** (the rate
at which the horizontal eye under jitter falls below `W_MIN` UI, gain-invariant).
It generates `fig_coopt2d_fixed.png`, `fig_contorno_2d.png`, `tabla_coopt2d.csv`,
`coopt2d_resultados.json`.

`coopt_ctle_dfe.py` co-optimizes CTLE+DFE with the objective = viable rate of the chain
**with DFE** (config D), expanding the space to (ω_z, c_post, **A_DC**); it compares
EQ-RC / C-opt / joint. It generates `fig_coopt_{equal_il,fixed}.png`, `tabla_coopt.csv`,
`coopt_resultados.json`.

`barrido_ojo_horizontal.py` (module `jitter.py`) estimates the **horizontal opening**
of the eye: bathtub curve (BER vs sampling phase) and its closure due to **jitter** (Gaussian
RJ + dual-Dirac DJ, parameters in `parametros.py::Jitter`). It generates
`fig_banera_canales.png`, `fig_banera_escenarios.png` and `tabla_ojo_horizontal.csv`.

`barrido_tasa_fisico.py` reproduces the robustness methodology of the notebook:
**30 realizations** (different PRBS pattern + noise) → viable symbol rate
(BER ≤ 1e-2) interpolated as **mean ± standard deviation**, for *no EQ / FFE+CTLE / +DFE*,
comparing RC vs FR-4 vs Megtron 6 under **two length criteria**
(`equal_il` = equal loss at Nyquist; `fixed` = 10 in physical). It generates:

| File | Deliverable |
| :--- | :--- |
| `fig_barrido_tasa_{equal_il,fixed}.png` | BER vs GBaud (median + min–max envelope) |
| `fig_tasa_viable_{equal_il,fixed}.png` | Viable rate mean ± std. (bars) |
| `tabla_tasa_viable.csv` | Summary vs the checkpoint RC reference (2.79/9.56/14.11) |

The sweep includes **two RC references**: `RC (cosh)` (exact) and `RC N=5`
(staircase from the notebook/checkpoint, to validate against 2.79/9.56/14.11 GBaud), and
sweeps up to **32 GBaud**. The DFE is **data-aided** (taps from least squares of
the data at its sampling phase — avoiding the phase mismatch that appears with the
propagation delay of the physical channel).

Support modules: `cadena_enlace.py` (frequency-domain chain: channel+CTLE via FFT,
alignment by FFT correlation, semi-analytic BER, data-aided DFE, custom DoE-LHS)
and `canales.py` (construction of RC cosh / RC N=5 / FR-4 / Megtron per criterion).
`diag_dfe.py` verifies the DFE. The function `cadena_enlace.doe_optimizar(...)` allows
**re-optimizing** (ω_z, c_post) over the physical channel (next step).

> **`bitacora_resultados.md`** records models, methodology, issues/corrections
> (including the DFE bug), results and reproducibility — material for the presentation.

`generar_comparativa.py` generates in `./figuras/`:

| File | Deliverable |
| :--- | :--- |
| `fig_atenuacion_dbporpulgada.png` | Dispersive attenuation (skin + dielectric) vs frequency, FR-4 vs Megtron 6 |
| `fig_respuesta_pulso.png` | Pulse response: the article's diffusive RC channel vs dispersive lines |
| `tabla_comparativa.csv` | Summary (f_-3dB, dB/in @Nyquist, peak, FWHM) to cite in the manuscript |

`simulacion_ojo.py` reproduces the chain of the notebook `simulacion_canal_rcEjecutadoActual.ipynb`
(PRBS23 → PAM-4 → FFE → **channel** → CTLE → RX LPF → DFE → eye/BER) **changing only the
channel** to the dispersive physical model, with the fixed DoE co-design from the checkpoint:

| File | Deliverable |
| :--- | :--- |
| `fig_ojo_senal_transmitida.png` | Transmitted signal x(t) vs received y(t) (transient), scenarios A and C |
| `fig_ojo_escenarios.png` | Eyes A (no EQ) / B (FFE) / C (FFE+CTLE) + sampled levels C vs D(+DFE) |
| `fig_ojo_comparativa_canal.png` | Eye scenario C: diffusive RC vs FR-4 vs Megtron 6 |
| `tabla_ojo_ber.csv` | BER and eye opening per channel and scenario (A/B/C/D) |

The comparison is **fair**: each transmission line is trimmed to the length
that matches the **same insertion loss at Nyquist (2 GHz)** as the RC channel,
so that the visible difference comes from the **dispersion profile**, not the
total loss.

---

## Package structure

| File | Role |
| :--- | :--- |
| `parametros.py` | Physical constants, substrate presets (FR-4/Megtron 6), geometry, operating point (4 GBaud, τ=400 ps) |
| `modelo_rlgc.py` | RLGC transmission line: R(f) skin, G(f) dielectric, γ(f), Zc(f), H(f), attenuation dB/in |
| `modelo_rc.py` | Distributed diffusive RC channel `H=1/cosh(√(sτ))` (article baseline) |
| `respuesta_pulso.py` | Pulse response via FFT and metrics (peak, FWHM) |
| `generar_comparativa.py` | **Main program**: produces the 2 figures + the table |
| `cadenas_busqueda_extension.md` | Search strings **C7–C9** aligned with C1–C6 of the manuscript |
| `requirements.txt` / `environment.yml` | Isolated local environment |

---

## Implemented physics (summary)

Per-unit-length parameters:

```
R(f) = sqrt(R_dc^2 + (k_skin·√f)^2)      # skin effect  (~√f)
L    = Z0 / v_p                          # inductance
G(f) = 2π·f·C·tan(δ)                     # dielectric loss (~f)
C    = 1 / (Z0·v_p)                      # capacitance
γ(f) = sqrt((R+jωL)(G+jωC)) = α + jβ     # propagation constant
```

- **RC baseline** (`modelo_rc.py`): diffusive attenuation (`exp(-√f)`), without
  inductance or skin.
- **New model** (`modelo_rlgc.py`): dispersive attenuation with a term linear
  in `f` (dielectric) + a `√f` term (conductor) — the real behavior of
  a PCB trace.

Modeling references already present in the manuscript: `5152908` (Wang & Wang,
RLC) and `554593` (Celik & Cangellaris, dispersive lines via Padé).

---

## Optional extensions

1. **Jitter:** inject RJ/DJ at the TX to estimate horizontal eye
   opening (natural extension of `respuesta_pulso.py`).
2. **Cross-validation:** export `Zc(f)`/`γ(f)` to Touchstone and
   correlate with SPICE/ADS — supported by search string **C9**.
