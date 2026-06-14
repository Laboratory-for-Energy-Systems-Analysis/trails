# Publication Reproduction Guide

This directory contains the notebooks, scripts, figures, and input data used for
the reproducibility material accompanying:

**A graph-matrix hybrid approach for deep temporalisation in time-explicit LCA**

Romain Sacchi, Tom Terlouw, Arthur Jakobs, Karin Treyer, Alvaro Hahn-Menacho,
and Christian Bauer.

The files below support the figures, tables, validation checks, and
supplementary calculations reported in the main article and supplementary
information. Run commands from the repository root unless stated otherwise.

## Environment

Follow the installation instructions in the main repository `README.md`. A
typical setup for reproducing these materials is:

```bash
conda create -n trails-publication python=3.11
conda activate trails-publication
pip install trails
```

The publication workflows also use notebooks, plotting, Excel import, Brightway,
and sparse-solver packages:

```bash
pip install jupyter matplotlib pandas plotly openpyxl bw2io psutil
```

If you use a Mac with Apple Silicon, you may need to install `scikit-umfpack` from conda-forge
with:
```bash
conda install -c conda-forge scikit-umfpack
```

To generate a prospective Premise-made scenario package, install and configure `premise` 
with access to the required licensed ecoinvent data.

Developers working from a local source checkout can instead install the current
checkout with:

```bash
pip install -e .[testing]
```

## Scenario Package Generation

Open and run:

```text
dev/publication/2.1. generate Trails data package.ipynb
```

Purpose:

- Creates a TRAILS-compatible data package from a configured Brightway/premise
  setup.
- Produces the scenario ZIP used by the case-study and stress-test workflows,
  depending on the notebook configuration.

Required resources:

- licensed ecoinvent data
- a configured Brightway project
- `premise`
- the scenario settings defined in the notebook

This notebook is a data-preparation recipe. It is not expected to run on a clean
machine without a complete local LCA database setup.

## Data Resources

`LCIs/example data package/`

Small Frictionless data package used by the gasoline car example and the
limiting-case validation. It contains `datapackage.json`, classifications, and
year-specific `A` and `B` matrix files for 2005, 2020, 2050, and 2100.

`LCIs/lci-pass_cars.xlsx`

Additional LCIs are available under LCIs/ to reproduce Figure 5.
Foreground inventory for the battery-electric passenger car case.

`LCIs/lci-case-study-ccu_polyol_delayed_release.xlsx`

Foreground inventory for the captured-CO2 polyol case.

`LCIs/lci-case-study-marine_fuel_switch.xlsx`

Foreground inventory for the marine fuel-switch case.

`LCIs/lci-case-study-daccs_storage_risk.xlsx`

Foreground inventory for the direct-air-capture case.

`LCIs/lci-case-study-biomass_growth_vs_gas_heat.xlsx`

Additional foreground inventory for biomass growth timing experiments. It is
not part of the default four-case Figure 5 notebook.

`LCIs/lci-case-study-fertilizer_n2o_timing.xlsx`

Additional foreground inventory for fertilizer N2O timing experiments. It is
not part of the default four-case Figure 5 notebook.

`trails_remind_SSP2-PkBudg1000.zip`

TRAILS-compatible REMIND SSP2-PkBudg1000 data package used by the main
case-study notebook and solver stress test. If this file is absent, recreate it
with `2.1. generate Trails data package.ipynb` or obtain the corresponding
publication data package.

`trails_SSP2-NPi.zip`


## Main Article

### Figure 1

With the conda environment activated, run:

```bash
python dev/publication/generate_algorithm_flow_diagram.py
```

Input resources: none beyond the plotting dependencies.

Outputs:

- `dev/publication/algorithm_flow_diagram.png`
- `dev/publication/algorithm_flow_diagram.svg`

The PNG is the publication figure included in this directory. The SVG is useful
for editing or high-resolution export.

### Figures 2, 3, and 4

With the conda environment activated, open ``jupyter lab`` and run:

```text
dev/publication/1. simple numerical example.ipynb
```

Input resources:

- `dev/publication/LCIs/example data package/datapackage.json`
- bundled LCIA and FaIR resources installed with `trails`

Outputs:

- Figure 2: temporal supply-chain graph for the gasoline car example, exported
  by the notebook as `trails_graph.html`.
- Figure 3: annual and cumulative GWP100 scores for the gasoline car example,
  displayed in the notebook.
- Figure 4: FaIR radiative forcing and temperature response for the gasoline car
  example, displayed in the notebook.
- Additional diagnostic output: `trails_adaptive_sankey.html`.

### Table 1

Table 1 summarises the functional units, system boundaries, and temporal
specifications of the case-study inventories. The corresponding foreground
inventory files are:

- `dev/publication/LCIs/lci-pass_cars.xlsx`
- `dev/publication/LCIs/lci-case-study-ccu_polyol_delayed_release.xlsx`
- `dev/publication/LCIs/lci-case-study-marine_fuel_switch.xlsx`
- `dev/publication/LCIs/lci-case-study-daccs_storage_risk.xlsx`

### Figure 5

With the conda environment activated, open ``jupyter lab`` and run:

```text
dev/publication/3. static vs foreground vs deep routing cases.ipynb
```

Input resources:

- `dev/publication/trails_remind_SSP2-PkBudg1000.zip`
- the four case-study Excel inventories listed for Table 1
- `scikit-umfpack` for direct sparse solves where selected by the notebook

Outputs:

- `dev/publication/notebook_runs/temporal_lci_routing_comparison_clean/routing_comparison_scores.csv`
- PNG and HTML routing-comparison figures under
  `dev/publication/notebook_runs/temporal_lci_routing_comparison_clean/<case>/`

The Figure 5 panels are selected from the generated routing-comparison figures:

- polyol case: EF v3.1 material resources, metals/minerals
- battery-electric vehicle case: EF v3.1 freshwater ecotoxicity
- marine fuel-switch case: EF v3.1 ozone depletion
- direct-air-capture case: EF v3.1 human toxicity, carcinogenic

The notebook writes both stacked and unstacked versions for each case and LCIA
method. Use the stacked PNGs when reproducing the main Figure 5 layout with
annual embodied impacts grouped by first-tier activity.

## Supplementary Information

### Section 5: Validation Against Limiting Cases

Run:

```bash
python dev/publication/validate_limiting_cases.py
```

Input resources:

- `dev/publication/LCIs/example data package/datapackage.json`

Outputs:

- `dev/publication/validation_limiting_cases_results.json`
- `dev/publication/validation_limiting_cases_summary.csv`

This script reproduces the limiting-case checks reported in the supplementary
information: static equivalence, foreground-only routing, identical matrix
slices, and the anchor-year versus target-year amount interpretation.

### Section 6: Depth Sensitivity for the Direct-Air-Capture Case

Run:

```bash
python dev/publication/daccs_pm_routing_diagnostics.py depth-sweep
```

Input resources:

- `dev/publication/LCIs/lci-case-study-daccs_storage_risk.xlsx`
- `dev/publication/trails_remind_SSP2-PkBudg1000.zip`

Outputs:

- `dev/notebook_runs/daccs_pm_depth_sweep/daccs_pm_depth_sweep.csv`

This script regenerates the fixed-depth direct-air-capture particulate-matter
sensitivity table reported in the supplementary information. To reduce runtime
while testing the setup, pass a smaller depth list, for example:

```bash
python dev/publication/daccs_pm_routing_diagnostics.py depth-sweep --depths 1 2
```

### Section 7: Coupling Trails Inventories With FaIR

The simple numerical example notebook demonstrates the FaIR coupling used for
main-text Figure 4:

```text
dev/publication/1. simple numerical example.ipynb
```

The supplementary information documents the model coupling, flow-to-species
mapping, and sign conventions. The default FaIR input files and mapping are
provided by the installed `trails` package.

### Section 9: Solver Stress Test

Run:

```bash
python dev/publication/solver_stress_test.py
```

Input resources:

- `dev/publication/trails_remind_SSP2-PkBudg1000.zip`
- `scikit-umfpack`

Outputs:

- `dev/publication/monolithic_vs_trails_remind_SSP2-PkBudg1000.csv`

`dev/publication/esi_solver_stress_test.md` contains the publication-ready
narrative and table derived from this benchmark. Update that Markdown file after
rerunning the script if matrix dimensions, timings, memory estimates, or solver
failure thresholds change.

### Section 10: Adaptive-Routing Sankey Diagrams

The supplementary information includes four adaptive-routing Sankey diagrams:

- Figure S10.1: polyol case, EF v3.1 material resources, metals/minerals
- Figure S10.2: battery-electric passenger car case, EF v3.1 freshwater
  ecotoxicity
- Figure S10.3: direct-air-capture case, EF v3.1 ozone depletion
- Figure S10.4: marine shipping case, EF v3.1 human toxicity, carcinogenic

Run:

```bash
python dev/publication/generate_depth_sankeys.py
```

Input resources:

- `dev/publication/trails_remind_SSP2-PkBudg1000.zip`
- the four case-study Excel inventories listed for Table 1

Outputs:

- HTML Sankey diagrams under
  `dev/notebook_runs/temporal_lci_depth_sweep_runner/sankey/`
- `dev/notebook_runs/temporal_lci_depth_sweep_runner/sankey/sankey_summary.csv`

For the direct-air-capture adaptive-routing diagnostics, first run the Section 6
depth sweep, then run:

```bash
python dev/publication/daccs_pm_routing_diagnostics.py adaptive --write-sankey
```

This writes adaptive routing CSV summaries and optional activity-year Sankey
HTML files under `dev/notebook_runs/daccs_pm_depth_sweep/`.

## Reproducibility Notes

- The heavy workflows depend on ecoinvent-derived scenario packages. These files
  may not be redistributable with the source repository and may need to be
  obtained separately or regenerated locally.
- Annual interpolation can create a cache in the platform-specific user data
  directory used by `trails`. First-run performance may differ from reruns that
  reuse cached interpolation data.
- Existing notebook outputs may show paths from the machine on which they were
  last executed. Rerunning the notebooks refreshes those displayed paths and
  generated outputs.
