# =============================================================================
# 03_direct_estimates.R - design-based (Horvitz-Thompson) MCV1 estimates at
# admin-1, admin-2 and admin-3, with Taylor-linearised SEs, a bootstrap
# replicate-weight check at admin-2, and the reliability summaries reported in
# S3 Appendix.
#
# Outputs: output/estimates/direct_admin{1,2,3}.csv
#          output/tables/direct_reliability.csv
#          output/tables/direct_vs_regional_se.csv
#          output/tables/direct_taylor_vs_bootstrap.csv
# =============================================================================

log_msg("03 direct estimates")
children <- read_out("data", "children.rds") %>% filter(!is.na(value))
unit_admin2 <- read_csv_out("data", "unit_summary_admin2.csv")
unit_admin3 <- read_csv_out("data", "unit_summary_admin3.csv")

make_design <- function(strata_opt) {
  switch(strata_opt,
         none  = svydesign(ids = ~cluster + householdID, weights = ~weight, data = children, nest = TRUE),
         urban = svydesign(ids = ~cluster + householdID, weights = ~weight, strata = ~ur,
                           data = children, nest = TRUE),
         dhs   = svydesign(ids = ~cluster + householdID, weights = ~weight, strata = ~dhs_stratum,
                           data = children, nest = TRUE),
         stop("unknown strata option: ", strata_opt))
}

direct_level <- function(level) {
  key <- key_col(level)
  des <- make_design(cfg$direct$strata[[level]])
  est <- svyby(as.formula("~value"), as.formula(paste0("~", key)), des, svymean,
               vartype = c("se", "ci")) %>%
    as_tibble() %>%
    rename(est = value, ci_l = ci_l, ci_u = ci_u) %>%
    mutate(ci_u_capped = pmin(ci_u, 1))
  est
}

# admin-1 ------------------------------------------------------------------------
direct_admin1 <- direct_level("admin1") %>%
  left_join(children %>% group_by(ADM1_EN) %>%
              summarise(n_clusters = n_distinct(cluster), n_resp = n()), by = "ADM1_EN") %>%
  mutate(strata = cfg$direct$strata[["admin1"]]) %>%
  arrange(desc(est))
save_csv(direct_admin1, "estimates", "direct_admin1.csv")
cat("\nAdmin-1 direct estimates (%):\n")
print(direct_admin1 %>% transmute(ADM1_EN, est = round(100 * est, 1), ci = sprintf("%.1f-%.1f", 100 * ci_l, 100 * ci_u), n_resp))

# admin-2 / admin-3 ----------------------------------------------------------------
direct_lower <- function(level, unit_tbl) {
  key <- key_col(level)
  direct_level(level) %>%
    right_join(unit_tbl, by = key) %>%
    mutate(strata = cfg$direct$strata[[level]],
           near_zero_se = !is.na(se) & se < cfg$near_zero_se,
           reliable = !is.na(est) & !near_zero_se) %>%
    select(all_of(c(key, name_col(level))), ADM1_EN, area_id, n_clusters, n_resp, n_vacc,
           homogeneous, small_n, has_data, est, se, ci_l, ci_u, near_zero_se, reliable, strata) %>%
    arrange(area_id)
}
direct_admin2 <- direct_lower("admin2", unit_admin2)
direct_admin3 <- direct_lower("admin3", unit_admin3)

# Taylor vs bootstrap replicate weights (admin-2) ----------------------------------
set.seed(cfg$seed)
des2 <- make_design(cfg$direct$strata[["admin2"]])
rep2 <- as.svrepdesign(des2, type = cfg$direct$boot_type, replicates = cfg$direct$n_boot)
boot2 <- svyby(~value, ~ADM2_PCODE, rep2, svymean) %>% as_tibble() %>% select(ADM2_PCODE, se_boot = se)
direct_admin2 <- direct_admin2 %>% left_join(boot2, by = "ADM2_PCODE")

save_csv(direct_admin2, "estimates", "direct_admin2.csv")
save_csv(direct_admin3, "estimates", "direct_admin3.csv")

# --- reliability summaries (S3 Appendix) ----------------------------------------------
reliability <- function(d, level) {
  withdata <- d %>% filter(has_data)
  nz <- withdata %>% filter(near_zero_se)
  rel <- withdata %>% filter(reliable)
  tribble(
    ~level, ~quantity, ~value,
    level, "units with data", nrow(withdata),
    level, "near-zero SE (n)", nrow(nz),
    level, "near-zero SE (%)", round(100 * nrow(nz) / nrow(withdata), 1),
    level, "near-zero SE and n <= 5 (n)", sum(nz$small_n),
    level, "near-zero SE and homogeneous (n)", sum(nz$homogeneous),
    level, "flagged unreliable (%)", round(100 * (1 - nrow(rel) / nrow(withdata)), 1),
    level, "median SE among reliable (pp)", round(100 * median(rel$se), 1),
    level, "min SE among reliable (pp)", round(100 * min(rel$se), 1),
    level, "max SE among reliable (pp)", round(100 * max(rel$se), 1),
    level, "units with SE < 10 pp among reliable (%)", round(100 * mean(rel$se < 0.10), 1)
  )
}
direct_reliability <- bind_rows(reliability(direct_admin2, "admin2"), reliability(direct_admin3, "admin3"),
                                tibble(level = "admin1", quantity = "median SE (pp)",
                                       value = round(100 * median(direct_admin1$se), 1)))
save_csv(direct_reliability, "tables", "direct_reliability.csv")
print(as.data.frame(direct_reliability))

# lower-level SE compared with the regional SE ------------------------------------------
regional_se <- direct_admin1 %>% select(ADM1_EN, se_admin1 = se)
vs_regional <- function(d, level) {
  d %>% filter(has_data) %>%
    left_join(regional_se, by = "ADM1_EN") %>%
    mutate(diff_pp = 100 * (se - se_admin1)) %>%
    group_by(level = level, ADM1_EN) %>%
    summarise(n_units = n(),
              max_increase_pp = round(max(diff_pp), 1),
              median_diff_pp = round(median(diff_pp), 1),
              share_below_regional_se = round(100 * mean(se < se_admin1), 1),
              share_nonzero_se = round(100 * mean(!near_zero_se), 1),
              share_increase_gt10pp = round(100 * mean(diff_pp > 10), 1), .groups = "drop") %>%
    arrange(level, desc(max_increase_pp))
}
direct_vs_regional <- bind_rows(vs_regional(direct_admin3, "admin3"), vs_regional(direct_admin2, "admin2"))
save_csv(direct_vs_regional, "tables", "direct_vs_regional_se.csv")

# Taylor vs bootstrap --------------------------------------------------------------------
tb <- direct_admin2 %>% filter(has_data, !is.na(se_boot))
tb_rel <- tb %>% filter(reliable, se_boot > 0)
tb_row <- function(d, label) tibble(
  subset = label, n = nrow(d),
  correlation = round(cor(d$se, d$se_boot), 3),
  mean_taylor_pp = round(100 * mean(d$se), 2), mean_boot_pp = round(100 * mean(d$se_boot), 2),
  ratio_of_means = round(mean(d$se) / mean(d$se_boot), 3),
  mean_ratio = round(mean(d$se / d$se_boot), 3), median_ratio = round(median(d$se / d$se_boot), 3),
  pct_diff_ratio_of_means = round(100 * (mean(d$se) / mean(d$se_boot) - 1), 1))
taylor_boot <- bind_rows(tb_row(tb, "all units with data"),
                         tb_row(tb_rel, "reliable units (both SEs > 0)"))
save_csv(taylor_boot, "tables", "direct_taylor_vs_bootstrap.csv")
print(as.data.frame(taylor_boot))
log_msg("03 direct estimates: done")
