# repoArticulo — Reproducibility package

Article: *Equalization Co-Design for Maximizing Viable Symbol Rate in
Bandwidth-Limited High-Speed Interconnect Links*
(distributed RC channel + dispersive FR-4/Megtron 6 lines; FFE+CTLE+DFE).
A **simulation** study (Phase 1), with cross-validation in LTspice.

This package is **self-contained**: it holds the manuscript, the code that
generates each figure and table, the intermediate data, and the validation in
an independent circuit simulator, so that any result is **traceable and
reproducible** from its source. The exact result→code correspondence is in
[`TRAZABILIDAD.md`](TRAZABILIDAD.md).

The manuscript is written in English and formatted with Elsevier's
`elsarticle` document class (numbered / Vancouver reference style) for
submission to *AEÜ – International Journal of Electronics and Communications*.

## Structure

```
repoArticulo/
├── README.md                 (this file)
├── TRAZABILIDAD.md           figure/table → script → data matrix
├── LICENSE                   MIT for code; CC BY 4.0 for data/figures/manuscript
├── .gitignore
├── manuscrito/               LaTeX + PDF + bibliography + figures (.png)
│   ├── elsarticle-IJEC.tex   English manuscript (elsarticle, numbered/Vancouver)
│   ├── elsarticle-IJEC.pdf   compiled PDF
│   ├── tabla_parametros.tex  master parameter table (\input)
│   ├── references.bib        bibliography
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
pdflatex elsarticle-IJEC && bibtex elsarticle-IJEC && \
pdflatex elsarticle-IJEC && pdflatex elsarticle-IJEC   # elsarticle + bibtex (elsarticle-num)
# or simply: latexmk -pdf elsarticle-IJEC.tex
```

### 2. Simulation results (Python 3.14)
```bash
cd codigo_simulacion
pip install -r requirements.txt        # numpy, scipy, matplotlib
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

## Verified environment
- Python 3.14 · numpy 2.4 · scipy 1.17 · matplotlib 3.10.
- LTspice ADI 26.0.2.1 (batch mode `-b`).
- TeX Live 2025 · elsarticle 3.4c (numbered / `elsarticle-num` style).

## License
- **Code** (`codigo_simulacion/`, `validacion_ltspice/`): MIT.
- **Data, figures, and manuscript** (`datos/`, `*.png`, `manuscrito/`):
  CC BY 4.0.

See [`LICENSE`](LICENSE) for the full terms.

## Citing this package / data availability
The journal's data policy (Option C) asks for the underlying data and code to
be deposited in a repository with a persistent identifier and cited in the
article. To complete this:
1. Deposit this `repoArticulo/` in **Zenodo** or **Mendeley Data**.
2. Obtain the assigned **DOI**.
3. Add the DOI to the manuscript's *Data availability* statement and cite it as
   a `[dataset]` reference.
