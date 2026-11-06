# =============================================================================
# 04_urban_fraction.R - urban fraction per admin-2 / admin-3 unit from the GHSL
# settlement model (SMOD) and population grid (GHS-POP), used to post-stratify the
# model-based estimates.
#
# SMOD codes: 10 water, 11-13 rural, 21-23 urban cluster classes, 30 urban centre.
# Three fractions are computed for every unit (see cfg$urban):
#   urban_frac_area  share of the unit's AREA in classes >= 21 (sensitivity)
#   urban_frac_pop   share of the unit's GHS-POP POPULATION in cfg$urban$urban_classes
#   urban_frac_cal   urban_frac_pop shifted on the logit scale within each 2019
#                    survey region so that the region's population-weighted urban
#                    share equals the survey's design-weighted urban share of
#                    children 12-23 months (Addis Ababa: 1; Sidama and South West
#                    Ethiopia use the former SNNPR share)
# urban_frac is the one selected by cfg$urban$weighting.
#
# Output: output/data/urban_fraction.csv
#         output/tables/urban_calibration.csv  (per survey region: GHSL vs survey share)
# =============================================================================

log_msg("04 urban fraction from GHSL (epoch %d, weighting = %s, classes {%s})",
        cfg$urban$epoch, cfg$urban$weighting, paste(cfg$urban$urban_classes, collapse = ","))
adm2 <- read_out("data", "adm2.rds"); adm3 <- read_out("data", "adm3.rds")

smod_files <- if (cfg$urban$epoch == 2020) cfg$paths$ghsl_smod_2020 else cfg$paths$ghsl_smod_2030
pop_files  <- if (cfg$urban$epoch == 2020) cfg$paths$ghsl_pop_2020 else NULL
stopifnot(all(file.exists(smod_files)))
if (cfg$urban$weighting != "area" && is.null(pop_files)) stop("population weighting needs the 2020 GHS-POP tiles")

mosaic_tiles <- function(files) do.call(mosaic, c(lapply(files, rast), fun = "max"))
smod <- mosaic_tiles(smod_files)
eth_moll <- vect(st_transform(adm3, crs(smod)))
smod_eth <- mask(crop(smod, eth_moll), eth_moll)
codes <- c(10, 11, 12, 13, 21, 22, 23, 30)
mask_area <- classify(smod_eth, cbind(codes, as.numeric(codes >= 21)), others = NA)
mask_pop  <- classify(smod_eth, cbind(codes, as.numeric(codes %in% cfg$urban$urban_classes)), others = NA)

pop_eth <- NULL
if (!is.null(pop_files)) {
  stopifnot(all(file.exists(pop_files)))
  pop_eth <- mask(crop(mosaic_tiles(pop_files), eth_moll), eth_moll)
  pop_eth <- resample(pop_eth, smod_eth, method = "near")     # same grid; guards against 1-cell offsets
}

# calibration targets: design-weighted urban share of the analytic children per survey region
children <- read_out("data", "children.rds") %>% filter(!is.na(value))
target <- children %>% group_by(survey_region) %>% summarise(survey_urban_share = weighted.mean(urban, weight))
lookup <- read_csv_out("data", "survey_region_lookup.csv")

calibrate <- function(u, popw, survey_region) {
  out <- u
  for (r in unique(survey_region)) {
    i <- survey_region == r
    t <- target$survey_urban_share[target$survey_region == r]
    if (!length(t)) next
    if (t >= 0.999) { out[i] <- 1; next }
    ui <- pmin(pmax(u[i], 0.001), 0.999); pw <- popw[i]
    f <- function(delta) sum(pw * expit(logit(ui) + delta)) / sum(pw) - t
    delta <- tryCatch(uniroot(f, c(-15, 15))$root, error = function(e) 0)
    out[i] <- expit(logit(ui) + delta)
  }
  out
}

frac_for <- function(areas, level) {
  areas_m <- st_transform(areas, crs(smod_eth))
  fix <- function(x, fallback) ifelse(is.finite(x), x, fallback)
  out <- tibble(level = level, code = areas[[key_col(level)]], area_id = areas$area_id, ADM1_EN = areas$ADM1_EN,
                urban_frac_area = fix(exact_extract(mask_area, areas_m, fun = "mean", progress = FALSE), 0))
  if (!is.null(pop_eth)) {
    out$pop_total      <- exact_extract(pop_eth, areas_m, fun = "sum", progress = FALSE)
    out$urban_frac_pop <- fix(exact_extract(mask_pop, areas_m, fun = "weighted_mean", weights = pop_eth, progress = FALSE),
                              out$urban_frac_area)
    out$survey_region  <- lookup$survey_region[match(out$ADM1_EN, lookup$current_region)]
    out$urban_frac_cal <- calibrate(out$urban_frac_pop, out$pop_total, out$survey_region)
  } else {
    out$pop_total <- NA_real_; out$urban_frac_pop <- NA_real_; out$survey_region <- NA_character_; out$urban_frac_cal <- NA_real_
  }
  out %>% mutate(urban_frac = switch(cfg$urban$weighting, area = urban_frac_area, population = urban_frac_pop,
                                     population_calibrated = urban_frac_cal))
}
urban_fraction <- bind_rows(frac_for(adm2, "admin2"), frac_for(adm3, "admin3"))
save_csv(urban_fraction, "data", "urban_fraction.csv")

# per survey region: GHSL population shares vs the survey share (reviewer 2 comment 10)
if (!is.null(pop_eth)) {
  calib <- urban_fraction %>% filter(level == "admin3") %>% group_by(survey_region) %>%
    summarise(ghsl_pop_share_area_classes = round(weighted.mean(urban_frac_area, pop_total), 3),
              ghsl_pop_share_urban_classes = round(weighted.mean(urban_frac_pop, pop_total), 3),
              calibrated_share = round(weighted.mean(urban_frac_cal, pop_total), 3)) %>%
    left_join(target %>% mutate(survey_urban_share = round(survey_urban_share, 3)), by = "survey_region")
  save_csv(calib, "tables", "urban_calibration.csv")
  print(as.data.frame(calib))
}
print(urban_fraction %>% group_by(level) %>%
        summarise(n = n(), mean_area = round(mean(urban_frac_area), 3), mean_pop = round(mean(urban_frac_pop), 3),
                  mean_cal = round(mean(urban_frac_cal), 3), mean_selected = round(mean(urban_frac), 3)))
log_msg("04 urban fraction: done")
