# =============================================================================
# 06_validation.R - put direct, model-based and MBG estimates side by side and
# compute the comparison statistics reported in the paper:
#   * posterior SD with vs without the urban-rural covariate (Fig 1a-b, Fig B1)
#   * direct SE vs model SD, and |direct SE - model SD| vs sample size (Fig 2, Fig B2)
#   * direct vs model-based coverage, Pearson and Spearman (Fig 3)
#   * model-based vs IHME MBG coverage and SE (Fig 4)
#
# Outputs: output/estimates/comparison_admin{2,3}.csv   one row per unit
#          output/tables/validation_correlations.csv
#          output/tables/sae_sd_summary.csv
# =============================================================================

log_msg("06 validation")
z975 <- qnorm(0.975)

read_sae <- function(level, model) read_csv_out("estimates", sprintf("sae_%s_%s.csv", level, model))

read_mbg <- function(level) {
  if (level == "admin2") {
    read_csv(cfg$paths$mbg_admin2, show_col_types = FALSE) %>%
      filter(vaccine == "mcv1", year == 2019) %>%
      transmute(code = ADM2_PCODE, mbg_mean = mean, mbg_lower = lower, mbg_upper = upper)
  } else {
    m <- read_csv(cfg$paths$mbg_admin3, show_col_types = FALSE, name_repair = "minimal")
    names(m)[1] <- "code"
    m %>% filter(vaccine == "mcv1", year == 2019) %>%
      transmute(code, mbg_mean = mean, mbg_lower = lower, mbg_upper = upper)
  }
}

comparison_table <- function(level) {
  key <- key_col(level); nm <- name_col(level)
  d <- read_csv_out("estimates", sprintf("direct_%s.csv", level)) %>%
    select(code = all_of(key), name = all_of(nm), ADM1_EN, area_id, n_clusters, n_resp, n_vacc,
           homogeneous, has_data, direct_est = est, direct_se = se, direct_ci_l = ci_l,
           direct_ci_u = ci_u, near_zero_se, reliable)
  u  <- read_sae(level, "urban") %>%
    select(code, sae_u_median = median, sae_u_mean = mean, sae_u_sd = sd,
           sae_u_lower = lower, sae_u_upper = upper, urban_frac)
  nu <- read_sae(level, "nourban") %>%
    select(code, sae_nu_median = median, sae_nu_mean = mean, sae_nu_sd = sd,
           sae_nu_lower = lower, sae_nu_upper = upper)
  d %>% left_join(u, by = "code") %>% left_join(nu, by = "code") %>%
    left_join(read_mbg(level), by = "code") %>%
    mutate(level = level, .before = 1,
           sae_u_se_ci = (sae_u_upper - sae_u_lower) / (2 * z975),   # CI-width based SE, as in Fig 4
           mbg_se_ci   = (mbg_upper - mbg_lower) / (2 * z975),
           abs_se_diff_urban   = abs(direct_se - sae_u_sd),
           abs_se_diff_nourban = abs(direct_se - sae_nu_sd))
}
comp <- lapply(setNames(LEVELS, LEVELS), comparison_table)
for (level in LEVELS) {
  stopifnot(!anyNA(comp[[level]]$sae_u_median), !anyNA(comp[[level]]$mbg_mean))
  save_csv(comp[[level]], "estimates", sprintf("comparison_%s.csv", level))
}

# --- correlations ---------------------------------------------------------------
cor_row <- function(level, what, x, y, subset_label = "all units") {
  ok <- !is.na(x) & !is.na(y)
  tibble(level = level, comparison = what, subset = subset_label, n = sum(ok),
         pearson = round(cor(x[ok], y[ok]), 3),
         spearman = round(cor(x[ok], y[ok], method = "spearman"), 3))
}
correlations <- bind_rows(lapply(LEVELS, function(level) {
  cc <- comp[[level]]; wd <- cc %>% filter(has_data); rel <- cc %>% filter(reliable)
  bind_rows(
    cor_row(level, "posterior SD: with urban vs without urban", cc$sae_u_sd, cc$sae_nu_sd),
    cor_row(level, "direct SE vs posterior SD (urban model)", wd$direct_se, wd$sae_u_sd, "units with data"),
    cor_row(level, "direct SE vs posterior SD (no-urban model)", wd$direct_se, wd$sae_nu_sd, "units with data"),
    cor_row(level, "direct vs model-based coverage (urban model)", wd$direct_est, wd$sae_u_median, "units with data"),
    cor_row(level, "direct vs model-based coverage (urban model)", rel$direct_est, rel$sae_u_median, "reliable direct estimates"),
    cor_row(level, "direct vs model-based coverage (no-urban model)", wd$direct_est, wd$sae_nu_median, "units with data"),
    cor_row(level, "model-based (urban) vs MBG coverage", cc$sae_u_median, cc$mbg_mean),
    cor_row(level, "model-based (urban) CI-based SE vs MBG CI-based SE", cc$sae_u_se_ci, cc$mbg_se_ci),
    cor_row(level, "model-based (urban) posterior SD vs MBG CI-based SE", cc$sae_u_sd, cc$mbg_se_ci)
  )
}))
save_csv(correlations, "tables", "validation_correlations.csv")
print(as.data.frame(correlations))

# --- in-sample agreement beyond correlation (reviewer 2, comment 6) -----------------------
agree_row <- function(level, cc, subset_label) {
  d <- cc$sae_u_median - cc$direct_est
  tibble(level = level, subset = subset_label, n = nrow(cc),
         bias_pp = round(100 * mean(d), 2), mae_pp = round(100 * mean(abs(d)), 2),
         rmse_pp = round(100 * sqrt(mean(d^2)), 2),
         pearson = round(cor(cc$direct_est, cc$sae_u_median), 3),
         spearman = round(cor(cc$direct_est, cc$sae_u_median, method = "spearman"), 3),
         cri95_coverage = round(100 * mean(cc$direct_est >= cc$sae_u_lower & cc$direct_est <= cc$sae_u_upper), 1))
}
agreement <- bind_rows(lapply(LEVELS, function(level) {
  cc <- comp[[level]]
  bind_rows(agree_row(level, cc %>% filter(has_data), "units with data"),
            agree_row(level, cc %>% filter(reliable), "reliable direct estimates"))
}))
save_csv(agreement, "tables", "agreement_in_sample.csv")
print(as.data.frame(agreement))

# --- benchmark: model estimates aggregated to the level the survey is designed for -----
# Population-weighted means of the unit medians per region (and nationally) compared
# with the design-based regional / national estimates (reviewer 2: an internal check
# with a design-unbiased reference).
direct1 <- read_csv_out("estimates", "direct_admin1.csv") %>% select(ADM1_EN, direct = est, direct_se = se)
uf <- read_csv_out("data", "urban_fraction.csv")
children <- read_out("data", "children.rds") %>% filter(!is.na(value))
nat_direct <- weighted.mean(children$value, children$weight)
bench_rows <- list(); bench_reg <- list()
for (level in LEVELS) {
  e <- read_sae(level, "urban") %>% left_join(uf %>% filter(level == !!level) %>% select(code, pop_total), by = "code")
  reg <- e %>% group_by(ADM1_EN) %>% summarise(model = weighted.mean(median, pop_total), pop = sum(pop_total)) %>%
    left_join(direct1, by = "ADM1_EN")
  bench_reg[[level]] <- reg %>% mutate(level = level, .before = 1, model = round(100 * model, 1), direct = round(100 * direct, 1),
                                       direct_se = round(100 * direct_se, 1), diff_pp = model - direct)
  bench_rows[[level]] <- tibble(level = level, n_regions = nrow(reg),
                                national_model_pct = round(100 * weighted.mean(reg$model, reg$pop), 1),
                                national_direct_pct = round(100 * nat_direct, 1),
                                region_pearson = round(cor(reg$model, reg$direct), 3),
                                region_bias_pp = round(100 * mean(reg$model - reg$direct), 1),
                                region_mae_pp = round(100 * mean(abs(reg$model - reg$direct)), 1),
                                region_rmse_pp = round(100 * sqrt(mean((reg$model - reg$direct)^2)), 1),
                                regions_within_direct_ci = sum(abs(reg$model - reg$direct) <= 1.96 * reg$direct_se))
}
save_csv(bind_rows(bench_rows), "tables", "benchmark_regional.csv")
save_csv(bind_rows(bench_reg), "tables", "benchmark_regional_detail.csv")
print(as.data.frame(bind_rows(bench_rows)))

# --- posterior SD summaries (Fig B1 reference lines, Fig 1 text) --------------------
sd_summary <- bind_rows(lapply(LEVELS, function(level) {
  cc <- comp[[level]]
  tibble(level = level,
         mean_sd_urban_pp = round(100 * mean(cc$sae_u_sd), 2),
         mean_sd_nourban_pp = round(100 * mean(cc$sae_nu_sd), 2),
         median_sd_urban_pp = round(100 * median(cc$sae_u_sd), 2),
         median_sd_nourban_pp = round(100 * median(cc$sae_nu_sd), 2),
         share_units_sd_above_10pp_urban = round(100 * mean(cc$sae_u_sd > 0.10), 1),
         share_units_urban_sd_smaller = round(100 * mean(cc$sae_u_sd < cc$sae_nu_sd), 1),
         mean_direct_se_pp_units_with_data = round(100 * mean(cc$direct_se[cc$has_data]), 2),
         mean_mbg_se_ci_pp = round(100 * mean(cc$mbg_se_ci), 2),
         mean_sae_se_ci_pp = round(100 * mean(cc$sae_u_se_ci), 2))
}))
save_csv(sd_summary, "tables", "sae_sd_summary.csv")
print(as.data.frame(sd_summary))
log_msg("06 validation: done")
