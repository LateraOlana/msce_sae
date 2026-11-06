# =============================================================================
# 11_summary.R - collect every number used in the manuscript text into one
# Markdown file.
#
# Output: output/results_summary.md
# =============================================================================

suppressPackageStartupMessages(library(knitr))
log_msg("11 summary")

md <- character()
add <- function(...) md <<- c(md, ..., "")
tbl <- function(df, digits = 3) add(kable(df, format = "pipe", digits = digits))
rd  <- function(dir, f) if (file.exists(out_path(dir, f))) read_csv_out(dir, f) else NULL

add("# MCV1 small-area estimation (Ethiopia, 2019 mini-EDHS): regenerated results",
    sprintf("Generated %s by `run_all.R` stage 11.", format(Sys.time(), "%Y-%m-%d %H:%M")),
    sprintf("Configuration: age_year = %s; direct-estimate strata = %s; iid term = %s; iid in prediction = %s; no-urban intercept = %s; priors = %s; posterior draws = %d; stability method = %s (B = %d).",
            if (is.null(cfg$age_year)) "all" else paste(cfg$age_year, collapse = ","),
            paste(names(cfg$direct$strata), cfg$direct$strata, sep = "=", collapse = ", "),
            paste(names(cfg$model$iid_term), cfg$model$iid_term, sep = "=", collapse = ", "),
            cfg$model$iid_in_prediction, cfg$model$nourban_intercept, cfg$model$priors,
            cfg$n_post_samples, cfg$stability$method, cfg$stability$n_draws))

# 1 data ---------------------------------------------------------------------------
add("## 1. Data description"); tbl(rd("tables", "data_description.csv"))

# 2 direct ---------------------------------------------------------------------------
add("## 2. Design-based (direct) estimates", "### Admin-1 (Fig C1)")
d1 <- rd("estimates", "direct_admin1.csv") %>%
  transmute(region = ADM1_EN, n_resp, n_clusters, `estimate (%)` = round(100 * est, 1),
            `95% CI` = sprintf("%.1f - %.1f", 100 * ci_l, 100 * ci_u), `SE (pp)` = round(100 * se, 1))
tbl(d1)
add("### Reliability at admin-2 / admin-3 (S3 Appendix)"); tbl(rd("tables", "direct_reliability.csv"))
add("### Lower-level SE compared with the regional SE"); tbl(rd("tables", "direct_vs_regional_se.csv"))
add("### Taylor linearisation vs bootstrap replicate weights (admin-2)"); tbl(rd("tables", "direct_taylor_vs_bootstrap.csv"))

# 3 model ---------------------------------------------------------------------------
add("## 3. Unit-level beta-binomial BYM2 model")
hy <- rd("tables", "sae_hyperparameters.csv") %>% select(model, parameter, mean, sd, `0.025quant`, `0.5quant`, `0.975quant`)
tbl(hy)
add("### Posterior SD summaries (Fig 1a-b, Fig B1)"); tbl(rd("tables", "sae_sd_summary.csv"))

# 4 validation ---------------------------------------------------------------------
add("## 4. Validation: correlations (Fig 1, Fig 3, Fig 4)"); tbl(rd("tables", "validation_correlations.csv"))
lo <- rd("tables", "loocv_summary.csv")
add("## 5. Leave-one-area-out cross-validation")
if (is.null(lo)) add("_not run yet (stage 07)_") else tbl(lo)

# 6 stability ----------------------------------------------------------------------
add("## 6. Bottom-quintile stability")
st <- rd("tables", "stability_summary.csv")
tbl(st %>% filter(!grepl("fine", group)) %>% select(level, method, group, band, n, share))
reg <- bind_rows(lapply(LEVELS, function(l) rd("estimates", sprintf("stability_%s.csv", l)) %>%
                          filter(highlight) %>% count(level, ADM1_EN, name = "highlighted_units")))
add(sprintf("### Units with >= %.0f%% probability of being in the bottom quintile (Fig B4), by region",
            100 * cfg$stability$highlight_prob)); tbl(reg)

# 6b manuscript tables ------------------------------------------------------------
unit2 <- rd("data", "unit_summary_admin2.csv"); unit3 <- rd("data", "unit_summary_admin3.csv")
d1raw <- rd("estimates", "direct_admin1.csv")
table1 <- bind_rows(
  tibble(level = "First (regions)", units = nrow(d1raw), with_data = nrow(d1raw), without_data = 0L,
         n_le_5 = 0L, homogeneous = 0L, clusters_median = median(d1raw$n_clusters), clusters_max = max(d1raw$n_clusters),
         responses_median = median(d1raw$n_resp), responses_max = max(d1raw$n_resp), responses_total = sum(d1raw$n_resp)),
  lapply(list(`Second (zones)` = unit2, `Third (woredas)` = unit3), function(u) {
    w <- u %>% filter(has_data)
    tibble(units = nrow(u), with_data = nrow(w), without_data = sum(!u$has_data),
           n_le_5 = sum(w$n_resp <= cfg$small_n), homogeneous = sum(w$homogeneous, na.rm = TRUE),
           clusters_median = median(w$n_clusters), clusters_max = max(w$n_clusters),
           responses_median = median(w$n_resp), responses_max = max(w$n_resp), responses_total = sum(w$n_resp))
  }) %>% bind_rows(.id = "level")
) %>% mutate(without_data_pct = round(100 * without_data / units, 1), .after = without_data)
# median uncertainty per level: design-based SE among units
# with data, posterior SD of the model with the urban covariate over all units
unc <- bind_rows(
  tibble(level = "First (regions)", median_direct_se_pp = round(100 * median(d1raw$se), 1), median_posterior_sd_pp = NA_real_),
  lapply(list(`Second (zones)` = "admin2", `Third (woredas)` = "admin3"), function(l) {
    cc <- rd("estimates", sprintf("comparison_%s.csv", l))
    tibble(median_direct_se_pp = round(100 * median(cc$direct_se[cc$has_data], na.rm = TRUE), 1),
           median_posterior_sd_pp = round(100 * median(cc$sae_u_sd), 1))
  }) %>% bind_rows(.id = "level"))
table1 <- table1 %>% left_join(unc, by = "level")
write_csv(table1, out_path("tables", "table1_sample_sizes.csv"))
add("## Table 1 (manuscript): sample sizes and units without data"); tbl(table1)

cmp2 <- rd("estimates", "comparison_admin2.csv"); cmp3 <- rd("estimates", "comparison_admin3.csv")
perf_row <- function(level, cc) {
  wd <- cc %>% filter(has_data); l <- if (is.null(lo)) NULL else lo %>% filter(level == !!level, model == "urban")
  d <- wd$sae_u_median - wd$direct_est
  g <- function(col) if (is.null(l) || !nrow(l)) NA else l[[col]]
  tibble(level = level,
         n_units_with_direct = nrow(wd),
         in_sample_pearson = round(cor(wd$direct_est, wd$sae_u_median), 2),
         in_sample_spearman = round(cor(wd$direct_est, wd$sae_u_median, method = "spearman"), 2),
         in_sample_bias_pp = round(100 * mean(d), 1),
         in_sample_mae_pp = round(100 * mean(abs(d)), 1),
         in_sample_rmse_pp = round(100 * sqrt(mean(d^2)), 1),
         in_sample_cri95_coverage = round(100 * mean(wd$direct_est >= wd$sae_u_lower & wd$direct_est <= wd$sae_u_upper), 1),
         loaocv_n = g("n_areas"), loaocv_pearson = g("pearson"), loaocv_spearman = g("spearman"),
         loaocv_bias_pp = g("mean_bias_pp"), loaocv_mae_pp = g("mae_pp"), loaocv_rmse_pp = g("rmse_pp"),
         loaocv_coverage_95 = g("coverage_95"),
         mbg_n = nrow(cc),
         mbg_pearson_coverage = round(cor(cc$sae_u_median, cc$mbg_mean), 2),
         mbg_spearman_coverage = round(cor(cc$sae_u_median, cc$mbg_mean, method = "spearman"), 2),
         mbg_bias_pp = round(100 * mean(cc$sae_u_median - cc$mbg_mean), 1),
         mbg_pearson_se = round(cor(cc$sae_u_se_ci, cc$mbg_se_ci), 2))
}
bench <- rd("tables", "benchmark_regional.csv")
table2 <- bind_rows(perf_row("admin2", cmp2), perf_row("admin3", cmp3)) %>%
  left_join(bench %>% transmute(level, bench_n_regions = n_regions, bench_pearson = region_pearson,
                                bench_bias_pp = region_bias_pp, bench_rmse_pp = region_rmse_pp,
                                bench_regions_within_ci = regions_within_direct_ci,
                                bench_national_model_pct = national_model_pct, bench_national_direct_pct = national_direct_pct),
            by = "level")
write_csv(table2, out_path("tables", "table2_performance.csv"))
add("## Table 2 (manuscript): model performance"); tbl(table2)

# --- additional tables ---------------------------------------------------------------------
extra_tables <- list(c("urban_calibration.csv", "Urban fraction: GHSL population shares vs survey shares by survey region"),
                     c("benchmark_regional.csv", "Benchmark: model estimates aggregated to regions vs design-based regional estimates"),
                     c("benchmark_regional_detail.csv", "Benchmark by region"),
                     c("flow_analytic_sample.csv", "Analytic sample flow and ascertainment"),
                     c("region_crosswalk_counts.csv", "Clusters by survey region and current region"),
                     c("weight_spread_by_stratum.csv", "Sampling-weight spread within the 21 survey strata"),
                     c("agreement_in_sample.csv", "In-sample agreement beyond correlation"),
                     c("stability_thresholds.csv", "Bottom-quintile stability by certainty threshold"),
                     c("sensitivity_summary.csv", "Sensitivity of the model-based estimates"),
                     c("pit_summary.csv", "Posterior predictive check: PIT of cluster counts"))
for (pair in extra_tables) {
  x <- rd("tables", pair[1])
  if (!is.null(x)) { add(paste0("## ", pair[2])); tbl(if ("formula" %in% names(x)) select(x, -formula) else x) }
}

