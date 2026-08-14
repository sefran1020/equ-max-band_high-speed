# Circuit-Model Fidelity in Equalizer Co-Design — Reproducibility Package

Article: *Circuit-Model Fidelity in Equalizer Co-Design: From Distributed RC
to Dispersive PCB Interconnects*
(distributed RC channel + dispersive FR-4/Megtron 6 lines; FFE+CTLE+DFE).
A reproducible simulation study with analytical, Touchstone and LTspice
cross-checks.

This package is **self-contained**: it holds the manuscript, the code that
generates each figure and table, the intermediate data, and the validation in
an independent circuit simulator, so that any result is **traceable and
reproducible** from its source. The exact result→code correspondence is in
[`TRAZABILIDAD.md`](TRAZABILIDAD.md).

The manuscript is written in English and prepared as a **new submission** to
the *International Journal of Circuit Theory and Applications* (Wiley), using
the official 2026 Wiley LaTeX template (`USG.cls`). The analyses originally
developed during an earlier IEEE Access review are preserved in
`codigo_simulacion/revision01/` as scientific provenance and are incorporated
into the Wiley manuscript.

Every simulation-result figure in the article is, byte for byte, the output of
the source listed in `TRAZABILIDAD.md`; none is retouched after being saved.
The conceptual workflow diagram is the only non-numerical illustration.

## Structure

```
equ-max-band_high-speed/
├── README.md                 (this file)
├── TRAZABILIDAD.md           figure/table → script → data matrix
├── manuscrito/               Wiley LaTeX source + PDF + bibliography + figures
│   ├── manuscript.tex        English manuscript (Wiley/IJCTA)
│   ├── manuscript.pdf        compiled review PDF
│   ├── USG.cls               Wiley 2026 document class
│   ├── COMPILAR.txt          exact XeLaTeX/BibTeX build sequence
│   ├── tabla_parametros.tex  master parameter table (\input)
│   ├── references.bib        bibliography
│   ├── figures/              article figures (English)
│   ├── images/               Wiley template assets
│   └── IJCTA_LaTeX_Source.zip  self-contained submission source package
├── codigo_simulacion/        frequency-domain chain (Python) producing the results
│   ├── *.py                  per-experiment scripts (see TRAZABILIDAD)
│   ├── revision01/           peer-review analyses retained and used by the article
│   │   ├── rNcM_*.py         one script per reviewer comment (N = reviewer,
│   │   │                     M = comment); see revision01/README.md
│   │   ├── figuras/, tablas/ their outputs
│   │   └── BITACORA_REVISION01.md  run log and decisions
│   ├── notebooks/            checkpoint notebooks (fig_bode, PAM-8, throughput)
│   ├── bitacora_resultados.md  technical log (models, bugs, results)
│   ├── requirements.txt / environment.yml  environment
│   └── README.md
├── validacion_ltspice/       cross-validation in a circuit simulator
│   ├── *.cir, bloques.lib, disp_lines.lib   netlists and subcircuits
│   ├── *.py                  stimuli, .raw reader, comparators, article figure
│   ├── README_LTSPICE.md     validation methodology and results
│   └── fig_*.png             evidence (eyes, Bode, vector-fit)
└── datos/
    ├── touchstone/*.s2p      dispersive-channel S-parameters (bridge to SPICE/ADS)
    └── tablas/*.csv          numerical tables backing those in the manuscript
```

## Reproduce

### 1. Manuscript (PDF)
```bash
cd manuscrito
xelatex -no-pdf manuscript.tex
bibtex manuscript
xelatex -no-pdf manuscript.tex
xelatex manuscript.tex              # Wiley 2026 template + BibTeX
```
The resulting file is `manuscript.pdf`. See `manuscrito/COMPILAR.txt` for the
same build sequence and template notes.

### 2. Simulation results (Python 3.11+)
```bash
cd codigo_simulacion
pip install -r requirements.txt        # numpy, scipy, matplotlib, pyDOE
python simulacion_ojo.py               # eyes + BER (eye/BER tables)
python barrido_tasa_fisico.py          # viable rate (rate-reach)
python sensibilidad_canal.py           # channel sensitivity
python validacion_ber_mc.py            # semi-analytical vs Monte Carlo BER validation
python exportar_touchstone.py          # S-parameters -> datos/touchstone/*.s2p
# (see TRAZABILIDAD.md for the rest: doe_canal_fisico, coopt_ctle_dfe, estres_jitter, …)
```
Figures are written to `codigo_simulacion/figuras/` and copied into the manuscript.

### 3. Cross-validation in LTspice
Requires LTspice (ADI). Reproduces Bode, the full chain and the eye diagram, and
cross-checks them against Python:
```bash
cd validacion_ltspice
python correr_ltspice.py               # orchestrates everything (locates LTspice; --ltspice PATH if needed)
python figura_articulo_ojos.py         # regenerates the manuscript's fig_ltspice_ojos.png
```
Details and results in `validacion_ltspice/README_LTSPICE.md`.

## Where the LTspice cross-validation appears in the Wiley manuscript

Section **Results → "Model validation and physical-channel
characterization"**: the paragraph after the channel-sensitivity figure, with
the figure `fig_ltspice_ojos.png` (frequency-domain chain vs LTspice vs
overlay). It extends the Touchstone cross-check already present in that
subsection.

## Configured environment
- Python 3.11 recommended; the code was also verified with Python 3.14.4.
- `requirements.txt` declares `numpy>=1.24,<3`, `scipy>=1.10`, `matplotlib>=3.7`, and `pyDOE>=0.3.8`.
- SciPy is optional and only needed if fitting/interpolation extensions are enabled.
- LTspice ADI 26.0.2.1 (batch mode `-b`).
- TeX Live 2024 or later with XeLaTeX and BibTeX for the Wiley manuscript build.

## License

The reproducibility package is released under the MIT License. See
[`LICENSE`](LICENSE).
