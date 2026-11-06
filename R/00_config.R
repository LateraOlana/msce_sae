# =============================================================================
# 00_config.R - central configuration for the MCV1 small-area-estimation
# workflow (Tesfaye et al., Ethiopia 2019 mini-EDHS).
#
# Everything that points to a file, or that changes what the analysis does, is
# set here. The scripts in R/ never hard-code paths.
#
# The default settings are those of the main analysis in the paper. Every
# alternative that appears in the sensitivity analyses (S2 Appendix, Table B1)
# is a named switch below with a comment saying what it does.
# =============================================================================

cfg <- list()

# --- locations ---------------------------------------------------------------
# This file is sourced with the working directory set to the repository root
# (run_all.R does this). The raw inputs are not part of the repository; see
# data/raw/README.md for what to obtain and where to put it. By default they
# are read from data/raw/. To read them from somewhere else either set the
# environment variable MCV1_DATA_DIR to that folder (same layout) or create
# R/config_local.R (ignored by git) and overwrite any entry of cfg$paths there.
cfg$final_dir <- normalizePath(getwd(), winslash = "/")
cfg$data_dir  <- Sys.getenv("MCV1_DATA_DIR", unset = file.path(cfg$final_dir, "data", "raw"))

ghsl_tiles <- c("R8_C22", "R8_C23", "R9_C22", "R9_C23")   # the four 1-km tiles covering Ethiopia
cfg$paths <- list(
  # survey extract: one row per child, see data/raw/README.md for the columns
  survey_csv = file.path(cfg$data_dir, "survey", "mini_edhs_2019_children.csv"),
  # DHS cluster locations of the 2019 mini-EDHS (ETGE81FL)
  dhs_gps    = file.path(cfg$data_dir, "dhs_gps", "ETGE81FL.shp"),
  # 2021 administrative boundaries of the Central Statistical Agency (woreda level)
  woreda_shp = file.path(cfg$data_dir, "boundaries", "Eth_Wereda_Bofed_CSA_2021.shp"),
  # GHSL settlement model (SMOD) and population (POP) 1-km tiles, Mollweide, same grid.
  # The 2020 epoch is the one used (survey year 2019); the 2030 SMOD tiles are only
  # read when cfg$urban$epoch is set to 2030.
  ghsl_smod_2030 = file.path(cfg$data_dir, "ghsl",
                             paste0("GHS_SMOD_E2030_GLOBE_R2023A_54009_1000_V2_0_", ghsl_tiles, ".tif")),
  ghsl_smod_2020 = file.path(cfg$data_dir, "ghsl",
                             paste0("GHS_SMOD_E2020_GLOBE_R2023A_54009_1000_V2_0_", ghsl_tiles, ".tif")),
  ghsl_pop_2020  = file.path(cfg$data_dir, "ghsl",
                             paste0("GHS_POP_E2020_GLOBE_R2023A_54009_1000_V1_0_", ghsl_tiles, ".tif")),
  # IHME model-based geostatistical (MBG) vaccine coverage estimates aggregated to
  # the second- and third-level units (used only in the comparison of stage 06)
  mbg_admin2 = file.path(cfg$data_dir, "mbg", "agg_summaries_all_admin2.csv"),
  mbg_admin3 = file.path(cfg$data_dir, "mbg", "agg_summaries_all.csv")
)

# All results are written under output/
cfg$out <- list(
  root      = "output",
  data      = "output/data",        # cleaned inputs (rds / csv)
  estimates = "output/estimates",   # direct and model-based estimates
  tables    = "output/tables",      # summary numbers used in the text
  fig_main  = "output/figures/main",
  fig_s2    = "output/figures/S2",
  fig_s3    = "output/figures/S3",
  logs      = "output/logs"
)

# --- data --------------------------------------------------------------------
# Children aged 12-23 months (age_year == 1). NULL keeps every child in the extract.
cfg$age_year <- 1

cfg$seed           <- 2025     # R RNG (bootstrap, Beta draws). INLA's sampler ignores seeds: see README
cfg$n_post_samples <- 4000     # joint posterior draws per model
cfg$near_zero_se   <- 0.005    # direct SE below this (proportion scale, 0.5 pp) = "near zero"
cfg$small_n        <- 5        # "very small sample" threshold

# --- design-based (direct) estimates -----------------------------------------
# Stratification used in svydesign(), per level. Options:
#   "none"  - unstratified
#   "urban" - urban/rural factor (main analysis, every level)
#   "dhs"   - the 21 survey strata (region x urban/rural) from the extract
cfg$direct <- list(
  strata = c(admin1 = "urban", admin2 = "urban", admin3 = "urban"),
  n_boot = 500,           # bootstrap replicates for the Taylor-vs-bootstrap SE check
  boot_type = "bootstrap" # survey::as.svrepdesign type: "bootstrap" (Canty-Davison) or "subbootstrap" (Rao-Wu)
)

# --- urban fraction used for post-stratification ----------------------------------
# epoch:   2020 (main analysis; survey year 2019) or 2030 (areal only)
# weighting:
#   "area"                   share of the unit's area in SMOD classes >= 21 (sensitivity)
#   "population"             share of the unit's GHS-POP population living in
#                            `urban_classes` (sensitivity, uncalibrated)
#   "population_calibrated"  as "population", then shifted on the logit scale within
#                            each 2019 survey region so that the region's population-
#                            weighted urban share equals the survey's design-weighted
#                            urban share of children (main analysis: it reproduces
#                            the design-based regional estimates when aggregated)
# urban_classes: SMOD classes counted as urban for the population fractions
#   (23 = dense urban cluster, 30 = urban centre; these match the survey's regional
#   urban shares best before calibration)
cfg$urban <- list(epoch = 2020, weighting = "population_calibrated", urban_classes = c(23, 30))

# --- unit-level beta-binomial BYM2 model (INLA) -------------------------------
cfg$model <- list(
  # Neighbourhood definition for the BYM2 graph: "queen" (shared border or vertex,
  # main analysis) or "rook" (shared border only; sensitivity).
  adjacency = "queen",
  # Scale U of the PC prior on the BYM2 precision, P(sd > U) = 0.01. 1 = INLA default
  # (main analysis). Other values are passed explicitly (sensitivity).
  bym2_prec_U = 1,
  # Extra f(area_iid, model = "iid") term, per fitted model. The main analysis
  # carries it in every model except the second-level model with the urban variable
  # (S1 Appendix explains why it is redundant with the unstructured BYM2 component).
  iid_term = c(admin2_urban = FALSE, admin2_nourban = TRUE,
               admin3_urban = TRUE,  admin3_nourban = TRUE),
  # FALSE: area predictions use the BYM2 effect only. TRUE adds the iid effect (its
  # precision is so large that the difference is negligible).
  iid_in_prediction = FALSE,
  # The model without urban/rural is fitted with NO intercept (y ~ -1 + ...),
  # so the BYM2 effect carries the overall level. TRUE adds an intercept.
  nourban_intercept = FALSE,
  # "inla_default" = INLA defaults (main analysis): bym2 precision PC(1, 0.01),
  #   mixing PC(0.5, 0.5); iid log-gamma(1, 5e-5); beta-binomial over-dispersion
  #   N(0, 0.4) on the logit scale.
  # "pc_as_described" = PC(1, 0.01) on the iid precision too.
  priors = "inla_default"
)

# --- leave-one-area-out cross-validation --------------------------------------
cfg$loocv <- list(
  run = TRUE,
  n_post_samples = 1000,
  # Only areas with a usable direct estimate (se > near_zero_se) are held out.
  levels = c("admin2", "admin3"),
  models = "urban"           # add "nourban" to cross-validate the model without the covariate too
)

# --- sensitivity analyses (stage 12) ------------------------------------------------
cfg$sensitivity <- list(
  run = TRUE,
  n_post_samples = 2000,            # draws per sensitivity fit (fewer than the main run)
  prec_U = c(0.5, 2),               # alternative PC-prior scales for the BYM2 precision
  highlight_prob = 0.70
)

# --- bottom-quintile stability ---------------------------------------------------
cfg$stability <- list(
  # "posterior"   = use the joint posterior draws saved in stage 05 (main analysis,
  #   as the Methods describe).
  # "beta_approx" = moment-match a Beta distribution to each unit's posterior median
  #   and SD and draw independently (B draws).
  method         = "posterior",
  n_draws        = 2000,
  highlight_prob = 0.70,     # units at or above this are "consistently under-performing"
  threshold_90   = 0.90      # exceedance probability P(coverage > 90%)
)

# --- figures -------------------------------------------------------------------
cfg$figures <- list(
  dpi      = 400,            # PNG resolution
  tiff     = TRUE,           # also write TIFFs of the main figures (journal format)
  tiff_dpi = 600,
  base_size = 16,
  shrinkage_max_n = NULL     # e.g. 23 to restrict Fig 2 to small samples (NULL = all)
)

# --- local overrides (not in git) ---------------------------------------------------
if (file.exists("R/config_local.R")) source("R/config_local.R")
