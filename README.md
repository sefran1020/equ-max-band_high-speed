# Equalizer Co-Design Transfer — Reproducibility Package

Article: *When Equalizer Co-Design Does Not Transfer: From Idealized RC to
Dispersive PCB Interconnects*
(distributed RC channel + dispersive FR-4/Megtron 6 lines; FFE+CTLE+DFE).
A reproducible simulation study with analytical, Touchstone and LTspice
cross-checks.

This package is **self-contained**: it holds the manuscript, the code that
generates each figure and table, the intermediate data, and the validation in
an independent circuit simulator, so that any result is **traceable and
reproducible** from its source. The exact result→code correspondence is in
[`TRAZABILIDAD.md`](TRAZABILIDAD.md).

The manuscript is written in English and formatted with the IEEE Access
LaTeX class (`ieeeaccess`). This repository corresponds to the IEEE Access submission version.

## Structure

```
repoArticulo/
├── README.md                 (this file)
├── TRAZABILIDAD.md           figure/table → script → data matrix
├── manuscrito/               LaTeX + PDF + bibliography + figures (.png)
│   ├── access.tex            English manuscript (IEEE Access)
│   ├── access.pdf            compiled PDF
│   ├── tabla_parametros.tex  master parameter table (\input)
│   ├── references.bib        bibliography
│   ├── autores/              author photographs for IEEE biographies
│   └── fig_*.png             article figures (English)
├── codigo_simulacion/        frequency-domain chain (Python) producing the results
│   ├── *.py                  per-experiment scripts (see TRAZABILIDAD)
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
pdflatex access && bibtex access && \
pdflatex access && pdflatex access   # IEEE Access + BibTeX
```

### 2. Simulation results (Python 3.11)
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

## Where the LTspice cross-validation appears in the article

Section **Results → "Model validation and physical-channel
characterization"**: the paragraph after the channel-sensitivity figure, with
the figure `fig_ltspice_ojos.png` (frequency-domain chain vs LTspice vs
overlay). It extends the Touchstone cross-check already present in that
subsection.

## Configured environment
- Python 3.11 recommended; `requirements.txt` declares `numpy>=1.24,<3`, `scipy>=1.10`, `matplotlib>=3.7`, and `pyDOE>=0.3.8`.
- SciPy is optional and only needed if fitting/interpolation extensions are enabled.
- LTspice ADI 26.0.2.1 (batch mode `-b`).
- TeX Live 2025 for the IEEE Access manuscript build.

## License

The reproducibility package is released under the MIT License. See
[`LICENSE`](LICENSE).
