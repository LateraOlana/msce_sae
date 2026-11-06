# =============================================================================
# 07_loocv.R - leave-one-area-out cross-validation of the unit-level model.
#
# Every area with a usable direct estimate (non-missing, SE not near zero) is
# held out in turn: its clusters' counts are set to NA, the model is refitted,
# and the post-stratified prediction for the held-out area is compared with
# its direct estimate (Pearson / Spearman correlation, RMSE in percentage
# points, 95%-interval coverage).
#
# CACHED. Each held-out area's result is stored in
# output/estimates/loocv_<level>_<model>.csv together with a fingerprint of
# every setting that influences it (age filter, model switches, priors,
# number of draws). On a rerun only areas with no cached row, or a row made
# under different settings, are refitted; with nothing changed the stage
# takes seconds. Delete the CSV to force a full rerun.
#
# One INLA fit per held-out area (about 25-40 s). Results are written after
# every fit, so an interrupted run keeps what it did.
#
# Outputs: output/estimates/loocv_<level>_<model>.csv
#          output/tables/loocv_summary.csv
# =============================================================================

source("R/sae_functions.R")
if (!isTRUE(cfg$loocv$run)) { log_msg("07 loocv: skipped (cfg$loocv$run = FALSE)") } else {

log_msg("07 leave-one-area-out cross-validation")
models <- cfg$loocv$models %||% "urban"
# Settings that change a held-out prediction. The direct estimate it is
# compared with is joined fresh each time, so strata are not part of the key.
loocv_hash <- rlang::hash(list(age_year = cfg$age_year, model = cfg$model, urban = cfg$urban,
                               n_post_samples = cfg$loocv$n_post_samples, inla = as.character(packageVersion("INLA"))))
pred_cols <- c("pred_median", "pred_mean", "pred_sd", "pred_lower", "pred_upper")
summary_rows <- list()

for (level in cfg$loocv$levels) for (model in models) {
  urban_model <- model == "urban"
  direct <- read_csv_out("estimates", sprintf("direct_%s.csv", level)) %>%
    select(code = all_of(key_col(level)), name = all_of(name_col(level)), ADM1_EN, area_id,
           n_clusters, n_resp, direct_est = est, direct_se = se, reliable)
  hold <- direct %>% filter(reliable)
  out_file <- out_path("estimates", sprintf("loocv_%s_%s.csv", level, model))

  cached <- if (file.exists(out_file)) {
    old <- read_csv(out_file, show_col_types = FALSE)
    if ("config_hash" %in% names(old)) old %>% filter(config_hash == loocv_hash) %>% select(code, all_of(pred_cols), config_hash)
    else NULL
  } else NULL
  preds <- if (is.null(cached)) tibble(code = character(), config_hash = character()) else cached
  todo <- hold %>% filter(!code %in% preds$code)
  log_msg("%s / %s: %d hold-out areas, %d cached, %d to fit", level, model, nrow(hold), nrow(hold) - nrow(todo), nrow(todo))

  write_current <- function() {
    res <- hold %>% inner_join(preds, by = "code") %>% mutate(level = level, model = model, .before = 1)
    write_csv(res, out_file, na = "")
    res
  }
  for (i in seq_len(nrow(todo))) {
    t0 <- Sys.time()
    fit <- fit_sae(level, urban_model, holdout_area = todo$area_id[i], n_samples = cfg$loocv$n_post_samples)
    pred <- fit$estimates %>% filter(area_id == todo$area_id[i]) %>%
      transmute(code = todo$code[i], pred_median = median, pred_mean = mean, pred_sd = sd,
                pred_lower = lower, pred_upper = upper, config_hash = loocv_hash)
    preds <- bind_rows(preds, pred)
    write_current()
    log_msg("  %3d/%d %-35s direct %5.1f  pred %5.1f  (%.0f s)", i, nrow(todo), substr(todo$name[i], 1, 35),
            100 * todo$direct_est[i], 100 * pred$pred_median, as.numeric(difftime(Sys.time(), t0, units = "secs")))
  }
  res <- write_current()
  stopifnot(nrow(res) == nrow(hold))
  summary_rows[[paste(level, model)]] <- tibble(
    level = level, model = model, n_areas = nrow(res),
    pearson = round(cor(res$direct_est, res$pred_median), 3),
    spearman = round(cor(res$direct_est, res$pred_median, method = "spearman"), 3),
    rmse_pp = round(100 * sqrt(mean((res$direct_est - res$pred_median)^2)), 1),
    mae_pp = round(100 * mean(abs(res$direct_est - res$pred_median)), 1),
    mean_bias_pp = round(100 * mean(res$pred_median - res$direct_est), 1),
    coverage_95 = round(100 * mean(res$direct_est >= res$pred_lower & res$direct_est <= res$pred_upper), 1))
  print(as.data.frame(summary_rows[[paste(level, model)]]))
}
loocv_summary <- bind_rows(summary_rows)
save_csv(loocv_summary, "tables", "loocv_summary.csv")
log_msg("07 loocv: done")
}
