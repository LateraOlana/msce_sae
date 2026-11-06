# =============================================================================
# sae_functions.R - the unit-level beta-binomial BYM2 model (INLA) and its
# post-stratified area predictions. Used by 05_unit_level_models.R,
# 07_loocv.R and 12_sensitivity.R.
#
# fit_sae(level, urban_model, holdout_area = NULL, ...)
#   level          "admin2" or "admin3"
#   urban_model    TRUE  -> y ~ -1 + urban + rural + f(area, bym2) [+ f(area, iid)]
#                  FALSE -> y ~ -1 + f(area, bym2) [+ f(area, iid)]   (no covariate)
#   holdout_area   area_id(s) whose cluster counts are set to NA (cross-validation)
#   adjacency      "queen" (default, cfg$model$adjacency) or "rook"
#   prec_U         scale U of the PC prior on the BYM2 precision (1 = INLA default)
#   include_iid    NULL = cfg$model$iid_term for this model; TRUE/FALSE to override
#   stratum_effect TRUE adds f(stratum, iid) over the 21 survey strata (urban model only)
#   urban_frac     vector (one value per area) overriding the configured urban fraction
#   cpo            TRUE also returns the PIT values (posterior predictive check)
#
# Returns a list: estimates (one row per area, proportion scale), draws
# (matrix n_area x n_samples of post-stratified coverage), fixed, hyper,
# formula, waic, pit, cpo_failure, cpu_seconds.
# =============================================================================

sae_areas <- function(level) if (level == "admin2") read_out("data", "adm2.rds") else read_out("data", "adm3.rds")

# Contiguity adjacency (binary), cached as an adjacency matrix per definition.
sae_graph <- function(level, adjacency = cfg$model$adjacency %||% "queen") {
  f <- out_path("data", sprintf("adjacency_%s_%s.rds", level, adjacency))
  if (!file.exists(f)) {
    areas <- sae_areas(level)
    nb <- suppressWarnings(poly2nb(areas, queen = (adjacency == "queen"), row.names = as.character(areas$area_id)))
    if (any(card(nb) == 0)) {
      msg <- paste(level, adjacency, ": areas without neighbours:", paste(areas[[name_col(level)]][card(nb) == 0], collapse = ", "))
      if (adjacency == "queen") stop(msg) else warning(msg)
    }
    saveRDS(nb2mat(nb, style = "B", zero.policy = TRUE), f)
  }
  inla.read.graph(readRDS(f))
}

# Cluster-level binomial counts with the area index of the requested level.
sae_cluster_data <- function(level) {
  key <- key_col(level)
  areas <- st_drop_geometry(sae_areas(level)) %>% select(all_of(key), area_id)
  read_out("data", "clusters.rds") %>% st_drop_geometry() %>%
    filter(n_resp > 0) %>%
    transmute(cluster, code = .data[[key]], y = n_vacc, N = n_resp,
              urban = as.integer(urban), rural = 1L - as.integer(urban),
              stratum_id = as.integer(as.character(dhs_stratum))) %>%
    left_join(areas, by = c("code" = key)) %>%
    mutate(area_iid = area_id)
}

sae_urban_fraction <- function(level) {
  uf <- read_csv_out("data", "urban_fraction.csv") %>% filter(level == !!level) %>% arrange(area_id)
  stopifnot(all(uf$area_id == seq_len(nrow(uf))))
  uf$urban_frac
}

# Survey stratum (1..21) to use for the urban and the rural prediction of each area:
# the area's current region is mapped to the 2019 survey region, then to that region's
# urban / rural stratum (Addis Ababa has no rural stratum; its urban one is used).
sae_area_strata <- function(level) {
  areas <- st_drop_geometry(sae_areas(level))
  sr <- read_csv_out("data", "survey_region_lookup.csv")
  st <- read_csv_out("data", "stratum_lookup.csv")
  pick <- function(region, u) {
    s <- st$dhs_stratum[st$survey_region == region & st$urban == u]
    if (!length(s)) s <- st$dhs_stratum[st$survey_region == region][1]
    s
  }
  reg <- sr$survey_region[match(areas$ADM1_EN, sr$current_region)]
  stopifnot(!anyNA(reg))
  list(urban = vapply(reg, pick, 1, u = 1L), rural = vapply(reg, pick, 1, u = 0L))
}

fit_sae <- function(level, urban_model, holdout_area = NULL,
                    n_samples = cfg$n_post_samples, seed = cfg$seed, verbose = FALSE,
                    adjacency = cfg$model$adjacency %||% "queen", prec_U = cfg$model$bym2_prec_U %||% 1,
                    include_iid = NULL, stratum_effect = FALSE, urban_frac = NULL, cpo = FALSE) {
  areas  <- sae_areas(level); n_area <- nrow(areas)
  graph  <- sae_graph(level, adjacency)
  dat    <- sae_cluster_data(level)
  if (!is.null(holdout_area)) dat$y[dat$area_id %in% holdout_area] <- NA
  if (stratum_effect && !urban_model) stop("stratum_effect is only implemented for the model with the urban covariate")

  model_key <- paste0(level, "_", if (urban_model) "urban" else "nourban")
  use_iid   <- if (is.null(include_iid)) isTRUE(cfg$model$iid_term[[model_key]]) else isTRUE(include_iid)
  hyper_iid <- if (cfg$model$priors == "pc_as_described")
    list(prec = list(prior = "pc.prec", param = c(1, 0.01))) else NULL
  hyper_bym2 <- if (prec_U != 1)
    list(prec = list(prior = "pc.prec", param = c(prec_U, 0.01)), phi = list(prior = "pc", param = c(0.5, 0.5))) else NULL
  hyper_stratum <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))

  rhs <- c(if (urban_model) "-1 + urban + rural" else if (isTRUE(cfg$model$nourban_intercept)) "1" else "-1",
           "f(area_id, model = 'bym2', graph = graph, hyper = hyper_bym2)",
           if (use_iid) "f(area_iid, model = 'iid', hyper = hyper_iid)",
           if (stratum_effect) "f(stratum_id, model = 'iid', hyper = hyper_stratum)")
  form <- as.formula(paste("y ~", paste(rhs, collapse = " + ")))

  t0 <- Sys.time()
  res <- inla(form, family = "betabinomial", Ntrials = dat$N, data = dat,
              control.compute = list(config = TRUE, waic = TRUE, cpo = cpo),
              control.predictor = list(compute = FALSE), verbose = verbose)
  # NOTE: INLA 25.6.7 does not honour `seed` here (verified); draws differ on
  # every call. Stage 05 saves them so later stages are reproducible.
  samp <- inla.posterior.sample(n_samples, res, seed = seed, num.threads = "1:1")
  cpu <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  rn  <- rownames(samp[[1]]$latent)
  lat <- vapply(samp, function(s) s$latent[, 1], numeric(length(rn)))   # latent x draws
  pick <- function(names) { i <- match(names, rn); if (anyNA(i)) stop("latent rows missing: ", names[is.na(i)][1]); lat[i, , drop = FALSE] }

  # BYM2 rows 1..n are the total (structured + unstructured) area effect.
  eta <- pick(paste0("area_id:", seq_len(n_area)))
  if (use_iid && isTRUE(cfg$model$iid_in_prediction)) eta <- eta + pick(paste0("area_iid:", seq_len(n_area)))

  if (urban_model) {
    w   <- if (!is.null(urban_frac)) urban_frac else sae_urban_fraction(level)
    stopifnot(length(w) == n_area)
    eta_u <- sweep(eta, 2, pick("urban:1")[1, ], "+")
    eta_r <- sweep(eta, 2, pick("rural:1")[1, ], "+")
    if (stratum_effect) {
      strata <- sae_area_strata(level)
      s_all <- pick(paste0("stratum_id:", seq_len(max(dat$stratum_id))))
      eta_u <- eta_u + s_all[strata$urban, , drop = FALSE]
      eta_r <- eta_r + s_all[strata$rural, , drop = FALSE]
    }
    cov <- w * expit(eta_u) + (1 - w) * expit(eta_r)
  } else {
    if (isTRUE(cfg$model$nourban_intercept)) eta <- sweep(eta, 2, pick("(Intercept):1")[1, ], "+")
    w <- NA_real_
    cov <- expit(eta)
  }

  codes <- areas[[key_col(level)]]; names_ <- areas[[name_col(level)]]; lvl <- level
  estimates <- tibble(level = lvl, model = if (urban_model) "urban" else "nourban",
                      area_id = seq_len(n_area), code = codes,
                      name = names_, ADM1_EN = areas$ADM1_EN,
                      median = apply(cov, 1, median), mean = rowMeans(cov), sd = apply(cov, 1, sd),
                      lower = apply(cov, 1, quantile, 0.025), upper = apply(cov, 1, quantile, 0.975),
                      urban_frac = w)
  list(estimates = estimates, draws = cov, fixed = res$summary.fixed, hyper = res$summary.hyperpar,
       formula = paste(deparse(form), collapse = ""), waic = res$waic$waic,
       # INLA's pit is P(Y <= y_obs) (non-randomised); cpo is P(Y = y_obs). The randomised
       # PIT for discrete data, P(Y < y) + U * P(Y = y), is built from both in stage 12.
       pit = if (cpo) res$cpo$pit else NULL, cpo_value = if (cpo) res$cpo$cpo else NULL,
       cpo_failure = if (cpo) res$cpo$failure else NULL, cluster_data = dat,
       cpu_seconds = cpu)
}
