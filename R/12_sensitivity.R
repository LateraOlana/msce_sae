# =============================================================================
# 12_sensitivity.R - sensitivity of the model-based estimates to modelling
# choices (reviewer 2, comments 2, 8, 10, 17, 18). Runs after stage 08 and before
# the figure stages (see run_all.R).
#
# For each level the model with the urban covariate is refitted under:
#   baseline refit          same specification (shows Monte Carlo noise); also
#                           returns the PIT values for the posterior predictive check
#   no urban-rural covariate
#   rook adjacency          shared border only, instead of queen
#   PC prior U = 0.5 / 2    scale of the PC prior on the BYM2 precision
#   extra IID term flipped  removed where the main model has it, added where not
#   survey-stratum effect   f(stratum, iid) over the 21 strata of the 2019 survey
#   areal urban fraction    post-stratification with the share of area (sensitivity)
#   uncalibrated fraction   population-weighted fraction without the survey-region calibration
# Each variant is compared with the main run (stage 05): correlation and mean /
# maximum absolute difference of the posterior medians, ratio of mean posterior
# SDs, and the units flagged as consistently under-performing.
#
# Outputs: output/tables/sensitivity_summary.csv
#          output/tables/pit_summary.csv, output/data/pit.rds   (Fig B6)
# =============================================================================

source("R/sae_functions.R")
if (!isTRUE(cfg$sensitivity$run)) { log_msg("12 sensitivity: skipped") } else {

log_msg("12 sensitivity analyses")
ns  <- cfg$sensitivity$n_post_samples
thr <- cfg$sensitivity$highlight_prob
uf  <- read_csv_out("data", "urban_fraction.csv")

hyp <- function(fit, pattern) { i <- grep(pattern, rownames(fit$hyper)); if (length(i)) round(fit$hyper$mean[i[1]], 4) else NA_real_ }

compare <- function(level, label, fit, base_est, base_hl) {
  e <- fit$estimates
  hl <- prob_bottom_quintile(fit$draws) >= thr
  tibble(level = level, variant = label, n_units = nrow(e),
         cor_median = round(cor(e$median, base_est$median), 3),
         mean_abs_diff_pp = round(100 * mean(abs(e$median - base_est$median)), 2),
         max_abs_diff_pp = round(100 * max(abs(e$median - base_est$median)), 2),
         mean_sd_pp = round(100 * mean(e$sd), 2),
         sd_ratio_vs_main = round(mean(e$sd) / mean(base_est$sd), 3),
         n_highlighted = sum(hl), n_highlighted_main = sum(base_hl), overlap_highlighted = sum(hl & base_hl),
         bym2_precision = hyp(fit, "Precision for area_id"), bym2_phi = hyp(fit, "Phi for area_id"),
         overdispersion = hyp(fit, "overdispersion"), iid_precision = hyp(fit, "Precision for area_iid"),
         stratum_precision = hyp(fit, "Precision for stratum_id"),
         waic = round(fit$waic, 1), formula = fit$formula)
}

rows <- list(); pits <- list()
for (level in LEVELS) {
  base_est <- read_csv_out("estimates", sprintf("sae_%s_urban.csv", level))
  base_hl  <- prob_bottom_quintile(read_out("estimates", sprintf("sae_draws_%s_urban.rds", level))) >= thr
  main_iid <- isTRUE(cfg$model$iid_term[[paste0(level, "_urban")]])
  run <- function(label, ...) {
    log_msg("  %s: %s", level, label)
    fit <- fit_sae(level, TRUE, n_samples = ns, ...)
    rows[[length(rows) + 1]] <<- compare(level, label, fit, base_est, base_hl)
    fit
  }
  f0 <- run("baseline refit (Monte Carlo noise)", cpo = TRUE)
  set.seed(cfg$seed)
  u <- runif(length(f0$pit))
  pits[[level]] <- tibble(level = level, cluster_row = seq_along(f0$pit),
                          y = f0$cluster_data$y, N = f0$cluster_data$N,
                          pit_nonrandomised = f0$pit, cpo = f0$cpo_value,
                          # randomised PIT for discrete outcomes: P(Y < y) + U P(Y = y)
                          pit = pmin(pmax(f0$pit - f0$cpo_value * (1 - u), 0), 1),
                          failure = f0$cpo_failure)
  f_nu <- fit_sae(level, FALSE, n_samples = ns)
  rows[[length(rows) + 1]] <- compare(level, "no urban-rural covariate", f_nu, base_est, base_hl)
  run("rook adjacency", adjacency = "rook")
  for (U in cfg$sensitivity$prec_U) run(sprintf("PC prior on BYM2 precision, U = %g", U), prec_U = U)
  run(if (main_iid) "extra IID term removed" else "extra IID term added", include_iid = !main_iid)
  run("survey-stratum random effect (21 strata)", stratum_effect = TRUE)
  ufl <- uf %>% filter(level == !!level) %>% arrange(area_id)
  run("areal urban fraction", urban_frac = ufl$urban_frac_area)
  run("population-weighted urban fraction, uncalibrated", urban_frac = ufl$urban_frac_pop)
}
sens <- bind_rows(rows)
save_csv(sens, "tables", "sensitivity_summary.csv")
print(as.data.frame(sens %>% select(level, variant, cor_median, mean_abs_diff_pp, max_abs_diff_pp, sd_ratio_vs_main,
                                    n_highlighted, overlap_highlighted, waic)))

pit_df <- bind_rows(pits)
save_rds(pit_df, "data", "pit.rds")
pit_summary <- pit_df %>% group_by(level) %>%
  summarise(n_clusters = n(), share_clusters_all_vaccinated = round(100 * mean(y == N), 1),
            share_clusters_none_vaccinated = round(100 * mean(y == 0), 1),
            mean_randomised_pit = round(mean(pit, na.rm = TRUE), 3),
            sd_randomised_pit = round(sd(pit, na.rm = TRUE), 3),            # uniform: 0.289
            share_pit_below_0.05 = round(100 * mean(pit < 0.05, na.rm = TRUE), 1),
            share_pit_above_0.95 = round(100 * mean(pit > 0.95, na.rm = TRUE), 1),
            ks_p_value_vs_uniform = round(suppressWarnings(ks.test(pit, "punif"))$p.value, 3),
            cpo_failures = sum(failure > 0, na.rm = TRUE))
save_csv(pit_summary, "tables", "pit_summary.csv")
print(as.data.frame(pit_summary))
log_msg("12 sensitivity: done")
}
