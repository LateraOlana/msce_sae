# Subnational measles vaccination coverage in Ethiopia

Code and derived outputs for

> Tesfaye et al. *Subnational measles vaccination coverage in Ethiopia:
> evaluation of model-based estimates from the 2019 Demographic and Health
> Survey.* (journal and DOI to be added on publication)

The workflow estimates first-dose measles vaccine (MCV1) coverage among
children aged 12 to 23 months at the first (region), second and third
administrative levels of Ethiopia from the 2019 Ethiopia Mini Demographic and
Health Survey (mini-EDHS). It produces the design-based (direct) estimates, the
unit-level beta-binomial BYM2 model estimates fitted with R-INLA, the
leave-one-area-out cross-validation, the comparison with the IHME geostatistical
estimates, the sensitivity analyses and every figure and table in the paper and
its appendices.

## What is in the repository

| path | content |
|---|---|
| `run_all.R` | runs the stages in order (see below) |
| `R/00_config.R` | every path and every analysis switch |
| `R/01_setup.R` | packages, helpers, plotting theme |
| `R/02_data_prep.R` ... `R/12_sensitivity.R` | the analysis stages |
| `R/sae_functions.R` | the INLA model (adjacency graph, fit, posterior draws, post-stratification) |
| `R/plot_functions.R` | the figure panels |
| `data/raw/README.md` | the raw inputs that have to be obtained and where to put them |
| `output/` | the derived outputs of the run reported in the paper (estimates, summary tables, figures) |

The raw inputs are **not** included. The survey microdata and the cluster
locations are distributed by the DHS Program under its own terms, so nothing
below the administrative-unit level is stored here. `data/raw/README.md` lists
each input, its source and the columns the scripts expect.

## Requirements

R 4.5.1 with the packages below (the versions used for the paper in brackets).
INLA is installed from `https://inla.r-inla-download.org/R/stable`.

dplyr (1.1.4), tidyr (1.3.1), readr (2.1.5), purrr (1.0.4), stringr (1.5.2),
tibble (3.3.0), forcats (1.0.0), sf (1.0.21, GEOS 3.13.1, GDAL 3.11.0, PROJ 9.6.0),
terra (1.8.54), exactextractr (0.10.0), spdep (1.3.13), survey (4.4.8),
INLA (25.6.7), ggplot2 (3.5.2), patchwork (1.3.2), ggpubr (0.6.1),
ggnewscale (0.5.2), ggrepel (0.9.6), ggtext (0.1.2), viridis (0.6.5),
scales (1.4.0), msm (1.8.2), hexbin (1.28.5), mgcv (1.9.3), knitr (1.50).

## Running

```
Rscript run_all.R            # every stage (about one hour the first time; stage 07 is about 50 minutes)
Rscript run_all.R fast       # everything except cross-validation (07) and the sensitivity analyses (12)
Rscript run_all.R 02 03      # only the stages whose numbers are given
Rscript run_all.R figures    # only the figure stages (09, 10)
```

Run from the repository root, or open `ethiopia-mcv1-sae.Rproj` in RStudio and
`source("run_all.R")`. Each stage reads what the earlier stages wrote to
`output/`, so once those files exist any stage can be rerun on its own.

Raw inputs are read from `data/raw/` by default. To keep them elsewhere, set the
environment variable `MCV1_DATA_DIR` to a folder with the same layout, or create
`R/config_local.R` (ignored by git) that overwrites entries of `cfg$paths`.

## Stages

| stage | script | what it does | main outputs |
|---|---|---|---|
| 02 | `02_data_prep.R` | children aged 12 to 23 months, cluster locations, point-in-polygon assignment of each cluster to the 2021 administrative units | `output/data/*.rds`, `unit_summary_admin{2,3}.csv`, `tables/data_description.csv`, `tables/flow_analytic_sample.csv`, `tables/weight_spread_by_stratum.csv` |
| 03 | `03_direct_estimates.R` | Horvitz-Thompson estimates with the survey package, Taylor standard errors, bootstrap check, reliability summaries | `estimates/direct_admin{1,2,3}.csv`, `tables/direct_*.csv` |
| 04 | `04_urban_fraction.R` | GHSL settlement model and population grids to a population-weighted, survey-calibrated urban fraction per unit (plus the areal and the uncalibrated alternatives) | `data/urban_fraction.csv`, `tables/urban_calibration.csv` |
| 05 | `05_unit_level_models.R` | four INLA beta-binomial BYM2 models (second and third level, with and without the urban variable), 4,000 joint posterior draws, post-stratification | `estimates/sae_<level>_<model>.csv`, `sae_draws_*.rds`, `tables/sae_hyperparameters.csv` |
| 06 | `06_validation.R` | agreement between direct, model-based and geostatistical estimates, regional benchmark | `estimates/comparison_admin{2,3}.csv`, `tables/validation_correlations.csv`, `tables/agreement_in_sample.csv`, `tables/benchmark_regional*.csv` |
| 07 | `07_loocv.R` | leave-one-area-out cross-validation, one refit per unit with a reliable direct estimate (cached per unit) | `estimates/loocv_*.csv`, `tables/loocv_summary.csv` |
| 08 | `08_stability.R` | probability of being in the bottom coverage quintile, exceedance of 90%, threshold sensitivity | `estimates/stability_admin{2,3}.csv`, `tables/stability_summary.csv`, `tables/stability_thresholds.csv` |
| 12 | `12_sensitivity.R` | refits under rook adjacency, other prior scales, the independent unit term added or removed, a survey-stratum effect, the two alternative urban fractions; posterior predictive check | `tables/sensitivity_summary.csv`, `tables/pit_summary.csv` |
| 09 | `09_figures_main.R` | Fig 1 to 5 | `figures/main/Fig*.png` (and `.tiff`) |
| 10 | `10_figures_supp.R` | Fig B1 to B6 (S2 Appendix), Fig C1 to C2 (S3 Appendix) | `figures/S2/`, `figures/S3/` |
| 11 | `11_summary.R` | every number used in the text, in one Markdown file | `output/results_summary.md` |

Stage 12 runs before the figure stages because Fig B6 uses its output.

## Outputs

* Estimates are stored as **proportions (0 to 1)** in every CSV; the figures convert them to percent.
* Units are identified by the 2021 P-codes of the Central Statistical Agency (`ADM2_PCODE`, `ADM3_PCODE`). `area_id` is the row index used for the INLA graph (units sorted by P-code).
* `output/estimates/comparison_admin{2,3}.csv` is the one-row-per-unit table behind every figure: direct estimate and standard error, both model specifications, geostatistical estimate, urban fraction.
* `output/estimates/sae_draws_*.rds` hold the joint posterior draws of the four models (units by draws), from which every later stage is reproducible.
* `output/tables/` holds the summary numbers that appear in the text, the tables and the appendices.
* The three figure folders hold the PNG versions of every figure and panel.

## Reproducibility

Every stage is deterministic except the posterior sampling in stage 05: the
INLA sampler in version 25.6.7 ignores its seed, so a rerun of stage 05 moves
the model medians by Monte Carlo noise (about half a percentage point) and can
move a few units across the 70% threshold of the bottom-quintile probability.
The draws of the reported run are saved in `output/estimates/sae_draws_*.rds`
and all later stages are reproducible from them. The sensitivity table includes
a refit of the main model that shows the size of this noise.

The cross-validation (stage 07) is cached per held-out unit with a fingerprint
of the model settings, so a rerun only refits units whose settings changed.

## Licence

The code and the derived outputs are released under the MIT licence (see
`LICENSE`). The raw inputs keep the terms of their providers.
