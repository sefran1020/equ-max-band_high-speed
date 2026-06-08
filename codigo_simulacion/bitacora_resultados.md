# Logbook — Recommendation 01: Dispersive physical channel
**Project:** RC Channel / IEEE INTERCON 2026 · **Folder:** `recomendacion01_canal_fisico/`
**Last updated:** run with sweep at 32 GBaud + ladder N=5 (checkpoint reference).

> Record document for the **presentation**. Summarizes what was modeled, how, which
> bugs arose and how they were resolved, the results, and how to reproduce everything.

---

## 1. Objective (Recommendation 01 from the senior team)
Replace the **theoretical (diffusive) RC** line of the manuscript with a model of a
**transmission line with frequency-dependent losses** (skin effect
`R(f)` + dielectric loss `G(f)`), moving from a **diffusive** channel to a real
**dispersive** one (FR-4 and Megtron 6). This closes the limitation declared in
`discusion.tex` ("RC model without inductance/skin-effect") and the future work item C3.

---

## 2. Channel models

| Channel | Model | Notes |
| :--- | :--- | :--- |
| **RC (cosh)** | `H(s)=1/cosh(√(sτ))`, τ=R·C=400 ps | EXACT distributed RC (diffusive) |
| **RC N=5** | 5-cell ladder (state space → `H(f)`) | **The one in the notebook/checkpoint**; validation reference |
| **FR-4** | RLGC: `R(f)=√(R_dc²+(k√f)²)`, `G(f)=ωC·tanδ` | Dk=4.3, tanδ=0.020 (dispersive) |
| **Megtron 6** | RLGC, same | Dk=3.6, tanδ=0.004 (low loss) |

`γ(f)=√((R+jωL)(G+jωC))=α+jβ`; attenuation `α` in dB/inch. The channel is applied
in frequency: `y = IFFT(H(f)·FFT(x))`.

**Two length criteria** for the dispersive lines:
- `equal_il`: length that matches the **RC loss at Nyquist** (7.38 dB @2 GHz) → FR-4 22.9 in, Megtron 36.9 in. *Equal-loss comparison; isolates dispersion.*
- `fixed`: common physical length of **10 in**. *Realistic PCB design comparison.*

---

## 3. Link chain (notebook replica)
`PRBS23 → PAM-4 → FFE(3 taps) → CHANNEL → +noise(σ=0.10 V, RX input) → CTLE → LPF RX → DFE → eye/BER`

- Co-design **inherited from the RC** (fixed DoE): ω_z=2π·2.35 GHz, c_post=−0.301; CTLE 1 zero + 2 poles (15/30 GHz), A_dc=2.5.
- **Noise referred to the RX INPUT** (before the CTLE) + finite-BW LPF → the CTLE does not gain SNR for free.
- **Semi-analytical BER** (Q-function, Gaussian fit of the eye), floor 1e-12.
- **Data-aided DFE** (see §5): taps via least squares of the data.
- Compute optimization: channel and CTLE in **frequency** (equivalent to `lsim`) and alignment via **one FFT** → enables 30 realizations × sweep.

---

## 4. Deliverables (scripts and figures)

| Script | Produces |
| :--- | :--- |
| `generar_comparativa.py` | `fig_atenuacion_dbporpulgada.png`, `fig_respuesta_pulso.png`, `tabla_comparativa.csv` |
| `simulacion_ojo.py` | `fig_ojo_senal_transmitida.png`, `fig_ojo_escenarios.png`, `fig_ojo_comparativa_canal.png`, `tabla_ojo_ber.csv` |
| `barrido_tasa_fisico.py` | `fig_barrido_tasa_{equal_il,fixed}.png`, `fig_tasa_viable_{equal_il,fixed}.png`, `tabla_tasa_viable.csv` |
| `diag_dfe.py` | DFE diagnostics (no figures) |

Modules: `parametros.py`, `modelo_rc.py`, `modelo_rlgc.py`, `cadena_enlace.py`, `canales.py`, `respuesta_pulso.py`.

---

## 5. Incident and fix log (methodological rigor)

> Useful for the presentation: shows critical validation of the simulator.

**B1 — Dispersive attenuation validated.** At equal IL at Nyquist, FR-4 loses ~2.5×
more than Megtron 6 at high frequency (tanδ 0.020 vs 0.004). Dispersive pulse with
**real propagation delay** (≈3–5 ns) that the lumped RC does not have.

**B2 — Sampling alignment.** The first fast aligner (full-resolution FFT
correlation) sampled at the **symbol edge** (eye crossing) → BER
~0.5 everywhere. *Fix:* find the coarse delay via FFT and choose the **optimal phase**
(0..sps−1) by correlation with the actual symbols.

**B3 — DFE divergent on the dispersive channel (the big one).** The DFE designed from the
pulse response took the taps at the **phase of the design pulse peak**,
which differs from the **data sampling phase** (due to the propagation delay,
absent in the RC). Result: spurious taps (~−0.30 where the real ISI was ~0) →
the DFE **injected** ISI and diverged to BER 0.5 at every rate. The RC did not fail because
without delay both phases coincide.
*Diagnosis (`diag_dfe.py`):* FR-4 @1 GBaud gave `C BER=1.7e-9` but `D BER=5.6e-2`.
*Definitive fix:* **data-aided DFE** — estimate the post-cursors via least
squares directly from the aligned samples (`m_k ≈ h0·a_k + Σ h_i·a_{k−i}`),
so they match the real ISI by construction. After the fix: `D BER ≤ C BER` always.
*Note for the paper:* the data-aided DFE is more robust than the ZF-from-pulse one and, moreover,
exposes a phenomenon that **only appears with the physical model** (propagation delay).

---

## 6. Results — rate sweep (30 realizations, BER≤1e-2, mean±std, up to 32 GBaud)

> **PIPELINE VALIDATION:** the **RC N=5** ladder reproduces the notebook
> checkpoint within **~3–5%** → validates channel+CTLE in frequency, data-aided DFE,
> and the 30-realization methodology.
>
> | Config | RC N=5 (this pipeline) | Checkpoint | Δ |
> | :--- | --: | --: | --: |
> | No EQ | 2.88 ± 0.01 | 2.79 | +3% |
> | FFE+CTLE | 9.76 ± 0.18 | 9.56 | +2% |
> | FFE+CTLE+DFE | 13.73 ± 0.29 | 14.11 | −3% |

**`equal_il` criterion (equal IL at Nyquist, 7.38 dB @2 GHz):**
| Config | RC cosh | RC N=5 | FR-4 (22.9 in) | Megtron 6 (36.9 in) | checkpoint ref. |
| :--- | --: | --: | --: | --: | --: |
| No EQ | 3.46 | 2.88 | 4.58 | 5.16 | 2.79 |
| FFE+CTLE | 15.21 | 9.76 | 9.31 | 10.83 | 9.56 |
| FFE+CTLE+DFE | 19.38 | 13.73 | 9.71 | 11.88 | 14.11 |

**`fixed` criterion (10 physical in):**
| Config | RC cosh | RC N=5 | FR-4 10in | Megtron 10in |
| :--- | --: | --: | --: | --: |
| No EQ | 3.46 | 2.88 | 12.75 | 32.0 (cap) |
| FFE+CTLE | 15.21 | 9.76 | 14.06 | 16.27 |
| FFE+CTLE+DFE | 19.38 | 13.73 | 27.19 | 32.0 (cap) |

*(typical std 0.0–0.5 GBaud; "cap" = railed at the sweep maximum, 32 GBaud.)*

**Reading:**
- **cosh vs N=5:** the `cosh` (EXACT distributed RC) is more optimistic (15.2/19.4) than the
  N=5 ladder (9.8/13.7): the notebook's N=5 is a low-order **conservative** approximation
  in the band that the EQ exploits. The paper uses N=5 → its figures are valid.
- **equal_il:** without EQ, the dispersive ones tolerate MORE (FR-4 4.6, Megtron 5.2) than the RC
  (2.9–3.5). With EQ they become comparable to RC N=5 (FR-4 9.3/9.7, Megtron 10.8/11.9);
  FR-4 +DFE (9.7) < RC N=5 +DFE (13.7) due to the more abrupt dispersive roll-off.
- **fixed:** real 10 in traces are VERY capable (FR-4 +DFE 27.2; Megtron ≥32).
  **Over-equalization** in Megtron: FFE+CTLE LOWERS the rate (16.3 < 32) due to CTLE noise.

### 6.1 Re-DoE: EQ re-optimized over the physical channel (`doe_canal_fisico.py`)

(ω_z, c_post) is re-optimized per substrate, maximizing the viable rate of **FFE+CTLE**.
Viable rate (GBaud, 30 realizations) — **EQ inherited from the RC** → **re-optimized EQ**:

| Criterion | Channel | ω_z reopt | c_post reopt | FFE+CTLE | +DFE |
| :--- | :--- | --: | --: | :--- | :--- |
| equal_il | FR-4 (22.9 in) | 2.35 GHz | −0.207 | 9.3 → **9.6** | 9.7 → 9.7 |
| equal_il | Megtron (36.9 in) | 2.65 GHz | −0.114 | 10.8 → **11.9** | 11.9 → 12.2 |
| fixed | FR-4 (10 in) | **5.35 GHz** | −0.149 | 14.1 → **24.2** | 27.2 → 24.5 |
| fixed | Megtron (10 in) | **9.85 GHz** | −0.266 | 16.3 → **32.0** (cap) | 32.0 → 32.0 |

**Messages:**
1. **The re-DoE is decisive on realistic channels (fixed):** FR-4 10 in FFE+CTLE
   **14.1 → 24.2** (+72%), Megtron **16.3 → 32** (+96%). The RC's EQ was leaving a lot of
   rate on the table. → *Equalization MUST be co-designed over the physical channel.*
2. **The optimal ω_z grows with the channel bandwidth:** 2.35 (RC) → 2.65 → 5.35
   → 9.85 GHz. Short low-loss lines are wideband and need the
   CTLE peaking at a much higher frequency than the RC's (~1 GHz BW).
3. **On very lossy channels (equal_il) the improvement is modest** (FR-4 +0.3, Megtron
   +1.1): the EQ is already near its ceiling; the bottleneck is the loss.
4. **Nuance (CTLE+DFE co-optimization):** the re-DoE maximized FFE+CTLE, not +DFE. On
   FR-4 fixed, +DFE DROPS with the re-optimized EQ (27.2 → 24.5): a very aggressive CTLE
   amplifies noise that the DFE would avoid. With the DFE present, a SOFTER CTLE is preferable
   → the C-only optimum ≠ the C+DFE optimum. (Future refinement: co-optimize CTLE+DFE.)

Re-optimized parameters stored in `figuras/doe_resultados.json`.

### 6.2 Horizontal eye opening with jitter (Phase α)

Bathtub curve (BER vs sampling phase) over the linearly equalized waveform;
dual-Dirac jitter RJ=0.010 UI rms + DJ=0.05 UI pp; 4 GBaud (UI=250 ps).
Eye width (UI · ps), **without → with** jitter:

| Channel | Scenario | without jitter | with jitter |
| :--- | :--- | :--- | :--- |
| RC N=5 | C | 0.506 · 127 ps | 0.488 · 122 ps |
| FR-4 | A (no EQ) | 0.205 · 51 ps | 0.193 · 48 ps |
| FR-4 | B (FFE) | **0.289 · 72 ps** | 0.275 · 69 ps |
| FR-4 | C (FFE+CTLE) | 0.182 · 45 ps | 0.158 · 40 ps |
| Megtron 6 | C | 0.154 · 39 ps | 0.129 · 32 ps |

**Consistency:** the floor of each bathtub matches the vertical BER of §R3
(RC <1e-12, FR-4 ~5e-6, Megtron ~4.5e-6) → validates the horizontal metric.

**Messages:**
1. **The horizontal eye reflects the EQ mismatch:** RC (CTLE matched) 0.51 UI
   vs dispersive ~0.15–0.18 UI. The EQ inherited from the RC also degrades the timing
   margin on the physical channel.
2. **On FR-4, the RC's CTLE NARROWS the horizontal eye** (B 0.29 → C 0.18 UI): the
   FFE opens it, but the poorly matched CTLE (peaking + noise) closes it again below
   "no EQ" — consistent with the vertical BER (C worse than B). Reinforces the re-DoE.
3. **Modest and uniform jitter cost** (~0.015–0.025 UI ≈ 6–10 ps), dominated by
   the DJ (dual-Dirac); the narrowest eyes (Megtron) lose proportionally more.
4. All eyes remain **open** horizontally with jitter (minimum Megtron
   0.129 UI = 32 ps). The 2D (vertical × horizontal) and the DFE effect remain as
   a refinement (Phase γ).

> Deliverables: `fig_banera_canales.png`, `fig_banera_escenarios.png`,
> `tabla_ojo_horizontal.csv`. Code: `jitter.py`, `barrido_ojo_horizontal.py`.

### 6.3 Joint CTLE+DFE co-optimization (Phase β)

DoE objective = viable rate of the FULL chain (config D), search space expanded to
(ω_z, c_post, A_DC). Viable rate of D (GBaud, mean±std, 15 realizations):

| Criterion | Channel | EQ-RC | C-opt | joint | ω_z/A_DC joint |
| :--- | :--- | --: | --: | --: | :--- |
| equal_il | FR-4 | 9.7 | 9.6 | **9.9** | 1.56 GHz / 2.62 |
| equal_il | Megtron | 11.8 | 12.1 | 12.1 | 3.81 GHz / 2.79 |
| fixed | FR-4 | 27.3 | **25.1** | **28.2** | 1.56 GHz / 2.62 |
| fixed | Megtron | 32.0\* | 32.0\* | 32.0\* | 6.81 GHz / 2.04 |
*(\*railed at 32 GBaud)*

**Messages (hypothesis confirmed):**
1. **Optimizing FFE+CTLE in isolation is counterproductive for the full chain.** On
   FR-4 `fixed`, C-opt gives **25.1**, even worse than the inherited EQ (27.3): maximizing
   the C eye pushes an aggressive CTLE (ω_z=4.94 GHz) that the DFE then suffers from (noise).
2. **The joint co-optimization uses a SOFTER CTLE and wins.** On FR-4 `fixed` it lowers
   ω_z from 4.94 to **1.56 GHz** (less peaking) and reaches **28.2** (the best), leaving the
   tail ISI to the DFE → CTLE+DFE must be tuned **together**, not in cascade.
3. On severe channels (equal_il) the improvement is modest (near the ceiling); on Megtron
   `fixed` all rail at 32 (no discrimination; raise the cap to resolve).

> Deliverables: `fig_coopt_{equal_il,fixed}.png`, `tabla_coopt.csv`,
> `coopt_resultados.json`. Code: `coopt_ctle_dfe.py` (+ `A_DC` as a variable in
> `cadena_enlace.aplicar_ctle_freq`/`evaluar_config`).

### 6.4 2D margin (vertical × horizontal) under jitter (Phase γ)

Metric: **2D viable rate** = rate at which the **DFE-aware** horizontal eye with jitter
(RJ=0.010, DJ=0.05 UI dual-Dirac) falls below 0.10 UI. `fixed` criterion,
5 realizations. Comparison EQ-RC / β-opt / γ-opt (optimizes the 2D margin):

| Channel | EQ-RC | β-opt | γ-opt (2D) |
| :--- | --: | --: | --: |
| FR-4 (10") | 27.3±0.6 | 27.9±0.7 | 27.3±1.2 |
| Megtron (10") | 32.0\* | 32.0\* | 32.0\* |
*(\*railed at 32 GBaud)*

**Messages (honest result):**
1. **With the DFE in the loop, the 2D margin under jitter is practically the same for
   the three EQs** (~27.3 GBaud on FR-4): the data-aided DFE cleans the eye at each sampling
   phase and **washes out the difference** between linear equalizations.
2. **The 2D viable rate ≈ the vertical viable rate** (FR-4 fixed +DFE ~27.2): at this
   jitter level, **the binding limit is the VERTICAL one (amplitude), not the horizontal
   one (timing)**. Moderate jitter is not the bottleneck.
3. **Important methodological revision:** a first version measured the PRE-DFE
   bathtub (input of the linear slicer) and overstated the role of the linear EQ in the
   horizontal margin (it gave 0.0 for EQs that delegate to the DFE); the DFE-aware
   version corrects this. *Lesson for the defense: measuring the horizontal eye at the
   correct point (post-DFE) is decisive.*

**Caveat:** result dependent on the jitter level; with more severe jitter
(larger RJ/DJ) the horizontal axis could bind and the 2D co-optimization gain
relevance → **jitter stress test** as future work.

> Deliverables: `fig_coopt2d_fixed.png`, `fig_contorno_2d.png` (DFE-aware bathtub),
> `tabla_coopt2d.csv`, `coopt2d_resultados.json`. Code: `coopt_2d_jitter.py`,
> `jitter.banera_dfe`.

### 6.5 Jitter stress test

Sweep of the jitter level (RJ, DJ) = scale·(0.01, 0.05) UI up to DJ=0.40 UI;
2D viable rate (DFE-aware) vs level, `fixed` criterion, EQ-RC vs γ-opt.

| DJ (UI) | FR-4 EQ-RC | FR-4 γ-opt | Megtron EQ-RC | Megtron γ-opt |
| --: | --: | --: | --: | --: |
| 0.05 | 27.0 | 27.4 | 32.0 | 32.0 |
| 0.20 | 26.0 | 26.0 | 32.0 | 32.0 |
| 0.30 | 22.9 | 22.6 | 29.0 | 32.0 |
| 0.325 | 19.9 | 2.6 | 11.3 | 32.0 |
| 0.35 | 9.8 | 0.0 | 3.8 | 11.4 |
| 0.40 | 0.0 | 0.0 | 0.0 | 0.0 |

**Conclusions:**
1. **Wide jitter margin:** no rate penalty up to DJ ≈ 0.20–0.25 UI;
   gradual decline to ~0.30 UI; **hard collapse between 0.33 and 0.38 UI** (the Diracs of the
   dual-Dirac at ±DJ/2 exceed the eye half-opening). 0.2–0.25 UI exceeds typical
   SerDes jitter budgets → timing is NOT the bottleneck.
2. **The EQ choice does not matter for the timing margin:** EQ-RC and γ-opt track
   closely up to ~0.30 UI; in the collapse zone they diverge **erratically and with
   huge variance** (not a real advantage) → the 2D co-optimization (Phase γ) provides no
   practical benefit; the **DFE governs** up to the collapse.
3. **Jitter model (dual-Dirac vs sinusoidal PJ):** the abrupt cutoff is due to the
   **dual-Dirac** (hard Diracs at ±DJ/2, worst-case bound). With **sinusoidal PJ**
   jitter (arcsine, more realistic) the cliff **smooths and shifts**:
   at DJ=0.35 UI it gives ~23 GBaud (vs 9.8 dual-Dirac) and at 0.40 UI ~11 GBaud (vs 0). Under
   a realistic model, the timing margin is even more robust.

> Deliverables: `fig_estres_jitter.png`, `fig_estres_modelo.png`,
> `tabla_estres_jitter.csv`. Code: `estres_jitter.py` (precomputed DFE-aware
> bathtubs; jitter via convolution, `jitter.aplicar_jitter(modelo=...)`).

---

## 7. Findings (messages for the presentation)

1. **At equal loss at Nyquist, the diffusive RC is the MOST pessimistic without EQ.** FR-4
   and Megtron tolerate more rate without equalization (5.2 > 4.6 > 3.5 GBaud) thanks to their
   cleaner pulse; the RC has a long ISI tail (exp(−√f)).

2. **The co-design optimized over the RC does not transfer to the dispersive channel.**
   With the same EQ, the RC rises more (15/19) than FR-4/Megtron (9–12): the CTLE is
   matched to the RC roll-off, not to the dispersive one (skin+dielectric, more abrupt).
   → Motivates **re-optimizing the DoE over the physical channel**.

3. **Physical over-equalization.** In `fixed`, Megtron 10 in (almost lossless) is
   viable at ≥24 GBaud without EQ, and the FFE+CTLE **lowers** it (16.3): the EQ adds noise
   where it is not needed. The DFE recovers it.

4. **The physical model reveals phenomena absent in the RC:** propagation delay
   (affects CDR/timing and broke a poorly designed DFE), and substrate-dependent
   dispersion (FR-4 vs Megtron).

---

## 8. Caveats
- Comparison tied to the matching criterion (equal-IL vs fixed); both views are reported.
- EQ **not re-optimized** for the physical channel (uses the RC's) → B2/B3 of §7 are over inherited EQ.
- **Data-aided (genie) DFE**: performance bound with training sequence; more capable than ZF-from-pulse.
- RC cosh vs N=5 ladder: magnitude differences due to model order (the ladder validates against the checkpoint).

---

## 9. Next steps
1. ✅ **Re-DoE over the physical channel** (`doe_canal_fisico.py`) — done (§6.1).
2. **Co-optimize CTLE+DFE jointly** (objective = viable rate of +DFE, not of
   FFE+CTLE) → resolves the nuance in §6.1.4 (softer CTLE when there is a DFE).
3. (Optional) DFE with more taps for long/dispersive lines.
4. Integrate into the manuscript: dB/in figure + search strings C7–C9 (`cadenas_busqueda_extension.md`).

---

## 10. Reproduce everything
```powershell
# entorno aislado (no toca el Python del sistema)
python -m venv .venv ; .\.venv\Scripts\Activate.ps1 ; pip install -r requirements.txt

python generar_comparativa.py     # atenuación + respuesta al pulso
python simulacion_ojo.py          # ojos + señal transmitida
python barrido_tasa_fisico.py     # tasa viable (32 GBaud, 30 realizaciones, 2 criterios)
python diag_dfe.py                # verificación del DFE (opcional)
```

---

## 11. Milestone 2 — Responses to the major revision (simulation, Phase 1)

New scripts that close the reviewer's observations (`recomendacioneMejora.txt`):

| Script | Observation | Key result |
| :--- | :--- | :--- |
| `validacion_ber_mc.py` (C1) | #5 (optimistic BER) | Semi-analytical ≈ Monte Carlo near the threshold (FR-4 C @16 GBaud: 2.4e-2 vs 2.3e-2); below it the semi-analytical is **conservative** (not optimistic). `fig_validacion_ber_mc.png`, `tabla_validacion_ber_mc.csv`. |
| `sensibilidad_canal.py` (C2) | #3 (channel validation) | Viable rate: tanδ 0.002→0.025 ⇒ 16.8→13.7; length 5″→32″ ⇒ 15.1→6.4; roughness (Hammerstad) ≤1µm ⇒ <0.6 GBaud (2nd order). `fig_sensibilidad_canal.png`, `tabla_sensibilidad_canal.csv`. |
| `barrido_extendido.py` (C3+C4) | #4 (threshold), #6 (saturated) | Sweep to **64 GBaud**: Megtron fixed No EQ **40.1±1.3** (resolves the ≥32); +DFE still saturates at 1e-2 (≥64), finite at 1e-4 (56.4). Viable rate at 1e-2/1e-3/1e-4: the **ordering is preserved**. `fig_sensibilidad_umbral.png`, `tabla_extendido_umbral.csv`. |
| `exportar_touchstone.py` (C5) | #3 (cross-validation) | S-params to `.s2p` (ref. 50Ω); **passivity OK** (|S11|²+|S21|²≤0.97); S21≈H of the chain (<0.01 dB ≤4 GHz). `FR4_10in.s2p`, `Megtron6_10in.s2p`, `fig_touchstone_check.png`. |
| `convergencia_nsym.py` (#2) | Nsym=1000 small | **Converged**: 1000→16000 sym changes the viable rate **<0.3%** (FR-4 C 13.78→13.74; +DFE 27.12→27.11; Megtron C 16.73→16.77); per-level σ stable (<5%). Justifies Nsym=1000 (per-level μ/σ converge; +30 realizations ≈7500 sym/level). `fig_convergencia_nsym.png`, `tabla_convergencia_nsym.csv`. |

**Integration into the manuscript:** Table I (Megtron fixed 40.1 / ≥64), new threshold
sensitivity Table, BER validation subsection (Fig. MC), channel sensitivity
(Fig.) and Touchstone cross-check in §III; discussion updated. The manuscript
recompiles cleanly (11 pp., 0 undefined refs).

> **Pending (non-compute):** Block E (verify 2025 refs + integrate
> `scopus_Cad10.bib`/`ieee_Cad10.bib`). **Future work:** Phases 2-3 (AWG/oscilloscope
> bench, FPGA Kria IBERT) — the `.s2p` is the bridge to SPICE/ADS.
