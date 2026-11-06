# =============================================================================
# 08_stability.R - how reliably can a unit be labelled "bottom quintile"?
#
# For every unit of the urban-covariate model we compute the probability that
# it falls in the lowest 20 % of coverage across all units, two ways:
#   beta_approx : moment-matched Beta(alpha, beta) per unit, independent draws
#                 (B = cfg$stability$n_draws, seed = cfg$seed)
#   posterior   : the joint posterior draws saved by stage 05 (what the Methods
#                 text describes)
# Units with probability >= cfg$stability$highlight_prob are the "consistently
# under-performing" units mapped in Fig B4. Exceedance probabilities
# P(coverage > 90 %) are added as an extra column.
#
# Outputs: output/estimates/stability_<level>.csv
#          output/tables/stability_summary.csv   (bands among point-estimate bottom-quintile units)
# =============================================================================

log_msg("08 bottom-quintile stability")

band_fine   <- function(p) cut(p, breaks = c(0, .3, .4, .5, .6, .7, .8, 1.0001),
                               labels = c("<30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", ">=80%"), right = FALSE)
band_coarse <- function(p) cut(p, breaks = c(0, .5, .6, .7, .8, 1.0001),
                               labels = c("<50%", "50-60%", "60-70%", "70-80%", ">=80%"), right = FALSE)

# prob_bottom_quintile() is defined in 01_setup.R

summary_rows <- list(); stab <- list()
for (level in LEVELS) {
  est   <- read_csv_out("estimates", sprintf("sae_%s_urban.csv", level))
  draws <- read_out("estimates", sprintf("sae_draws_%s_urban.rds", level))
  stopifnot(nrow(draws) == nrow(est))
  n <- nrow(est)

  # Beta approximation (moment matching on the posterior median and SD)
  set.seed(cfg$seed)
  n_eff <- pmax(est$median * (1 - est$median) / est$sd^2 - 1, 0.01)
  alpha <- 1 + est$median * n_eff
  beta  <- 1 + (1 - est$median) * n_eff
  M_beta <- replicate(cfg$stability$n_draws, rbeta(n, alpha, beta))

  thr <- cfg$stability$threshold_90
  res <- est %>%
    mutate(bottom_point = median <= quantile(median, 0.20, type = 8),
           prob_bottom_beta_approx = prob_bottom_quintile(M_beta),
           prob_bottom_posterior   = prob_bottom_quintile(draws),
           band_beta_approx = band_fine(prob_bottom_beta_approx),
           band_posterior   = band_fine(prob_bottom_posterior),
           p_gt90_beta_approx = 1 - pbeta(thr, alpha, beta),
           p_gt90_posterior   = rowMeans(draws > thr)) %>%
    mutate(prob_bottom = .data[[paste0("prob_bottom_", cfg$stability$method)]],
           highlight = prob_bottom >= cfg$stability$highlight_prob,
           stability_method = cfg$stability$method)
  stab[[level]] <- res
  save_csv(res, "estimates", sprintf("stability_%s.csv", level))

  for (method in c("beta_approx", "posterior")) {
    p <- res[[paste0("prob_bottom_", method)]]
    bp <- res %>% filter(bottom_point)
    pb <- bp[[paste0("prob_bottom_", method)]]
    summary_rows[[paste(level, method)]] <- bind_rows(
      tibble(level = level, method = method, group = "n units", band = "all", n = n, share = NA_real_),
      tibble(level = level, method = method, group = "n units", band = "bottom quintile by point estimate", n = nrow(bp), share = NA_real_),
      tibble(level = level, method = method, group = "n units", band = sprintf("prob >= %.0f%% (highlighted)", 100 * cfg$stability$highlight_prob),
             n = sum(p >= cfg$stability$highlight_prob), share = NA_real_),
      bp %>% mutate(band = band_coarse(pb)) %>% count(band, .drop = FALSE) %>%
        mutate(level = level, method = method, group = "point-estimate bottom quintile, by certainty band (coarse)",
               band = as.character(band), share = round(100 * n / nrow(bp), 1)),
      bp %>% mutate(band = band_fine(pb)) %>% count(band, .drop = FALSE) %>%
        mutate(level = level, method = method, group = "point-estimate bottom quintile, by certainty band (fine)",
               band = as.character(band), share = round(100 * n / nrow(bp), 1)),
      tibble(level = level, method = method, group = "exceedance", band = sprintf("mean P(coverage > %.0f%%)", 100 * thr),
             n = NA_integer_, share = round(100 * mean(res[[paste0("p_gt90_", method)]]), 2))
    ) %>% select(level, method, group, band, n, share)
  }
}
stability_summary <- bind_rows(summary_rows)
save_csv(stability_summary, "tables", "stability_summary.csv")

# --- sensitivity to the certainty threshold (reviewer 2, comment 11) ---------------------
thresholds <- bind_rows(lapply(LEVELS, function(level) {
  res <- stab[[level]]; bp <- res %>% filter(bottom_point)
  bind_rows(lapply(c(0.5, 0.6, 0.7, 0.8, 0.9), function(t)
    tibble(level = level, method = cfg$stability$method, threshold = t,
           n_units_at_or_above = sum(res$prob_bottom >= t),
           n_bottom_point_units = nrow(bp),
           n_bottom_point_stable = sum(bp$prob_bottom >= t),
           share_bottom_point_stable = round(100 * mean(bp$prob_bottom >= t), 1))))
}))
save_csv(thresholds, "tables", "stability_thresholds.csv")
print(as.data.frame(thresholds))
print(as.data.frame(stability_summary %>% filter(group != "exceedance", !grepl("fine", group))))
log_msg("08 stability: done")
