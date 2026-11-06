# =============================================================================
# 05_unit_level_models.R - fit the four unit-level models (admin-2 / admin-3,
# with / without the urban-rural covariate), draw joint posterior samples and
# post-stratify them to area-level MCV1 coverage.
#
# Outputs: output/estimates/sae_<level>_<model>.csv     (proportion scale)
#          output/estimates/sae_all_models.csv
#          output/estimates/sae_draws_<level>_<model>.rds  (n_area x n_samples)
#          output/tables/sae_hyperparameters.csv
# =============================================================================

source("R/sae_functions.R")
log_msg("05 unit-level models (INLA)")

fits <- list()
for (level in LEVELS) for (urban_model in c(TRUE, FALSE)) {
  model <- if (urban_model) "urban" else "nourban"
  log_msg("fitting %s / %s", level, model)
  fit <- fit_sae(level, urban_model)
  log_msg("  %s  (%.0f s)", fit$formula, fit$cpu_seconds)
  save_csv(fit$estimates, "estimates", sprintf("sae_%s_%s.csv", level, model))
  save_rds(fit$draws, "estimates", sprintf("sae_draws_%s_%s.rds", level, model))
  fits[[paste(level, model, sep = "_")]] <- fit
}

all_est <- bind_rows(lapply(fits, `[[`, "estimates"))
save_csv(all_est, "estimates", "sae_all_models.csv")

hyper <- bind_rows(lapply(names(fits), function(k) {
  f <- fits[[k]]
  bind_rows(as_tibble(f$fixed, rownames = "parameter"), as_tibble(f$hyper, rownames = "parameter")) %>%
    mutate(model = k, formula = f$formula, .before = 1)
}))
save_csv(hyper, "tables", "sae_hyperparameters.csv")

print(all_est %>% group_by(level, model) %>%
        summarise(n = n(), median_cov = round(100 * median(median), 1),
                  mean_sd_pp = round(100 * mean(sd), 2), min_sd_pp = round(100 * min(sd), 2),
                  max_sd_pp = round(100 * max(sd), 2)))
log_msg("05 unit-level models: done")
