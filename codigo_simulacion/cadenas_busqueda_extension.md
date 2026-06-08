# Search strings — bibliographic extension for Recommendation 01

> **Scope:** this document extends search strings **C1–C6** from
> `../IEEE-conference-template-062824/manuscrito_borrador.md` with three new
> strings (**C7–C9**) specific to the **realistic physical channel** model
> (frequency-dependent losses, skin effect, validation with S-parameters). It
> keeps the same Scopus `TITLE-ABS-KEY` format and adds the IEEE Xplore variant
> (command search).
>
> These strings **underpin** the material produced by this package
> (`generar_comparativa.py`) and close the limitation stated in
> `discusion.tex` ("RC model without inductance/skin-effect") and the future
> work item **C3** of the manuscript ("RLGC with frequency-dependent losses +
> validation against S-parameters").

---

## Coherence map (which claim each new string supports)

| String | Claim to support | Deliverable it backs |
| :--- | :--- | :--- |
| **C7** | The real channel is **dispersive**: R(f) from skin effect + G(f) from dielectric loss | `fig_atenuacion_dbporpulgada.png` |
| **C8** | Characterization of **substrates** (FR-4 vs Megtron 6): Dk/Df, dB/inch | table `tabla_comparativa.csv` (substrates) |
| **C9** | **Cross-validation** of the numerical model against S-parameters / SPICE/ADS | simulation-gap closure (Phase 1.3) |

---

## C7 — Transmission lines with frequency-dependent losses (skin + dielectric)

**Scopus (TITLE-ABS-KEY):**
```
TITLE-ABS-KEY ( ( "transmission line" OR "interconnect" OR "channel" )
  AND ( "skin effect" OR "frequency-dependent" OR "dispersive" OR "RLGC" )
  AND ( "dielectric loss" OR "loss tangent" OR "conductor loss" )
  AND ( "signal integrity" OR "insertion loss" OR "pulse response" ) )
```

**IEEE Xplore (Command Search):**
```
("transmission line" OR "interconnect") AND ("skin effect" OR "RLGC" OR "frequency-dependent")
  AND ("dielectric loss" OR "loss tangent") AND ("signal integrity" OR "insertion loss")
```

**Candidates already present in `references.bib` (reuse, do not search again):**
- `5152908` — Wang & Wang, *Modeling of Distributed RLC Interconnect and
  Transmission Line via Closed Forms and Recursive Algorithms* (IEEE TVLSI, 2010).
- `554593` — Celik & Cangellaris, *Simulation of dispersive multiconductor
  transmission lines by Padé approximation via the Lanczos process* (T-MTT, 1996).

---

## C8 — PCB substrate characterization (FR-4 / Megtron 6) and attenuation in dB/inch

**Scopus (TITLE-ABS-KEY):**
```
TITLE-ABS-KEY ( ( "FR-4" OR "FR4" OR "Megtron" OR "low-loss laminate" OR "PCB substrate" )
  AND ( "dielectric constant" OR "loss tangent" OR "Dk" OR "Df" )
  AND ( "high-speed" OR "high speed" OR "multi-gigabit" OR "Gb/s" )
  AND ( "insertion loss" OR "attenuation" OR "channel" ) )
```

**IEEE Xplore (Command Search):**
```
("FR-4" OR "Megtron" OR "low-loss laminate") AND ("loss tangent" OR "dielectric constant")
  AND ("insertion loss" OR "attenuation") AND ("high-speed" OR "Gb/s")
```

> Note: no specific substrate is cited yet in the manuscript; this string
> provides the missing reference to justify the Dk/Df values used in
> `parametros.py` (FR-4: eps_r=4.3, tanδ=0.020 · Megtron 6: eps_r=3.6, tanδ=0.004).

---

## C9 — Cross-validation: numerical model vs S-parameters / SPICE / ADS

**Scopus (TITLE-ABS-KEY):**
```
TITLE-ABS-KEY ( ( "S-parameter" OR "scattering parameter" OR "Touchstone" OR "Sparameter" )
  AND ( "transmission line" OR "interconnect" OR "channel model" )
  AND ( "SPICE" OR "ADS" OR "field solver" OR "macromodel" OR "vector fitting" )
  AND ( "validation" OR "correlation" OR "passivity" OR "causality" ) )
```

**IEEE Xplore (Command Search):**
```
("S-parameter" OR "Touchstone" OR "vector fitting") AND ("interconnect" OR "channel model")
  AND ("SPICE" OR "ADS" OR "field solver") AND ("validation" OR "correlation")
```

**Candidate already present in `references.bib`:**
- `11346689` — Raghuram et al., *SPICE Based Optimization of Equalization in
  Channel Simulation* (EPEPS, 2025) — supports the Phase 1.3 cross-validation.

---

## C10 — Jitter (RJ/DJ), dual-Dirac, bathtub, horizontal eye / eye width

**Scopus (TITLE-ABS-KEY):**
```
TITLE-ABS-KEY ( ( "random jitter" OR "deterministic jitter" OR "dual-Dirac" OR "bathtub curve" )
  AND ( "eye width" OR "horizontal eye" OR "timing margin" OR "eye opening" )
  AND ( "serial link" OR "SerDes" OR "high-speed" OR "interconnect" ) )
```

**IEEE Xplore (Command Search):**
```
("random jitter" OR "deterministic jitter" OR "dual-Dirac") AND ("eye width" OR "timing margin")
  AND ("serial link" OR "SerDes" OR "high-speed")
```

**Collected reference (added to `references.bib`, verified):**
- `805809` — Li, Wilstrup, Jessen, Petrich, *A new method for jitter decomposition
  through its distribution tail fitting* (Int. Test Conf. 1999, pp. 788–794,
  doi:10.1109/TEST.1999.805809). **Foundational reference for the dual-Dirac model and
  RJ/DJ decomposition**; cited in Method (timing margin), Results R6 and
  Discussion. The horizontal margin is further anchored in the statistical-eye
  strings **C6** (`7202856`, `10584882`).

---

## Suggested integration into the manuscript

1. **Discussion / Limitations:** replace "RC model (without inductance/skin-effect)"
   with a sentence citing **C7** and pointing to the dB/inch attenuation figure as
   evidence of the diffusive→dispersive shift.
2. **Future work (Conclusion):** anchor the item "RLGC with frequency-dependent
   losses + validation against S-parameters" to **C7** and **C9**.
3. **Methodology (future v2):** a "Physical channel model" subsection citing
   **C7–C8**, reusing the already-curated `5152908` and `554593`.
