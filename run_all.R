#!/usr/bin/env Rscript
# =============================================================================
# run_all.R - run the MCV1 small-area-estimation workflow.
#
#   Rscript run_all.R                 # every stage, in order
#   Rscript run_all.R 02 03           # only the stages whose number is given
#   Rscript run_all.R figures         # the two figure stages (09, 10)
#   Rscript run_all.R fast            # everything except 07 (cross-validation) and 12 (sensitivity)
#
# In RStudio: open ethiopia-mcv1-sae.Rproj, then source("run_all.R") or run the
# stages in R/ one after another (each stage reads what the previous ones saved
# under output/, so they can be run individually once their inputs exist).
#
# Stages
#   02_data_prep          survey extract + GPS + woreda polygons -> analysis data
#   03_direct_estimates   design-based (HT) estimates, SE reliability checks
#   04_urban_fraction     GHSL urban fraction per unit (post-stratification)
#   05_unit_level_models  beta-binomial BYM2 models in INLA (4 models)
#   06_validation         SE comparisons, correlations, MBG comparison
#   07_loocv              leave-one-area-out cross-validation (cached; ~50 min only when
#                         results are missing or the model settings changed)
#   08_stability          bottom-quintile stability, exceedance probabilities
#   12_sensitivity        adjacency / prior / stratum / urban-fraction sensitivity, PIT (~15 min)
#   09_figures_main       Fig 1-5
#   10_figures_supp       S2 Appendix (B1-B4) and S3 Appendix (C1-C2)
#   11_summary            results_summary.md: every number used in the text
# =============================================================================

# Work from the directory that holds this file (the repository root).
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
if (length(this_file)) setwd(dirname(normalizePath(this_file)))
if (!file.exists("R/00_config.R")) stop("run_all.R must be run from the repository root")

# 12_sensitivity runs before the figure stages (Fig B6 uses its PIT values)
stages <- c("02_data_prep", "03_direct_estimates", "04_urban_fraction",
            "05_unit_level_models", "06_validation", "07_loocv", "08_stability", "12_sensitivity",
            "09_figures_main", "10_figures_supp", "11_summary")

args <- commandArgs(trailingOnly = TRUE)
selected <- if (!length(args)) stages else {
  if (identical(args, "figures")) stages[grepl("^(09|10)", stages)]
  else if (identical(args, "fast")) stages[!grepl("^(07|12)", stages)]
  else {
    pick <- stages[substr(stages, 1, 2) %in% args]
    if (!length(pick)) stop("no stage matches: ", paste(args, collapse = " "))
    pick
  }
}

source("R/00_config.R")
source("R/01_setup.R")
dir.create(cfg$out$logs, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(cfg$out$logs, format(Sys.time(), "run_%Y%m%d_%H%M%S.log"))
cat("log:", log_file, "\n")

t_all <- Sys.time()
for (s in selected) {
  cat("\n==================== ", s, " ====================\n", sep = "")
  .stage_t0 <- Sys.time()      # unusual name: stage scripts run in this environment and may define t0
  ok <- tryCatch({ source(file.path("R", paste0(s, ".R")), echo = FALSE); TRUE },
                 error = function(e) { message("ERROR in ", s, ": ", conditionMessage(e)); FALSE })
  elapsed <- round(as.numeric(difftime(Sys.time(), .stage_t0, units = "mins")), 2)
  cat(sprintf("%s %s in %.2f min\n", s, if (ok) "finished" else "FAILED", elapsed))
  cat(sprintf("%s | %s | %s | %.2f min\n", format(Sys.time()), s, if (ok) "ok" else "FAILED", elapsed),
      file = log_file, append = TRUE)
  if (!ok) stop("stopping: stage ", s, " failed")
}
cat(sprintf("\nAll selected stages done in %.1f min\n", as.numeric(difftime(Sys.time(), t_all, units = "mins"))))
