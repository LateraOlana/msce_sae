# =============================================================================
# 02_data_prep.R - build the analysis data set.
#
# Inputs : survey extract (children), DHS cluster GPS, CSA 2021 woreda polygons
# Outputs: output/data/children.rds      one row per child (12-23 months)
#          output/data/clusters.rds      one row per cluster (sf points) with counts
#          output/data/adm1.rds, adm2.rds, adm3.rds   polygons with area_id
#          output/data/unit_summary_admin{2,3}.csv    sample size per unit
#          output/tables/data_description.csv         numbers for the text
# =============================================================================

log_msg("02 data prep: reading inputs")

# --- children ------------------------------------------------------------------
raw <- read_csv(cfg$paths$survey_csv, show_col_types = FALSE)
children <- raw %>%
  { if (is.null(cfg$age_year)) . else filter(., age_year %in% cfg$age_year) } %>%
  transmute(cluster, householdID, weight,
            urban = as.integer(urban),
            ur = factor(ifelse(urban == 1, "urban", "rural"), levels = c("urban", "rural")),
            dhs_stratum = factor(strata),
            age_year,
            value = mcv1,
            # ascertainment (reviewer 2 comment 3)
            has_card = has_vacc_card,
            mcv1_card   = !is.na(mcv_dose_from_card)   & mcv_dose_from_card   >= 1,
            mcv1_recall = !is.na(mcv_dose_from_recall) & mcv_dose_from_recall >= 1)

# --- cluster GPS and point-in-polygon ------------------------------------------
gps <- st_read(cfg$paths$dhs_gps, quiet = TRUE) %>%
  st_drop_geometry() %>%
  select(cluster = DHSCLUST, DHSYEAR, URBAN_RURA, survey_region = ADM1NAME, LATNUM, LONGNUM) %>%
  st_as_sf(coords = c("LONGNUM", "LATNUM"), crs = 4326, remove = FALSE)

woreda <- st_read(cfg$paths$woreda_shp, quiet = TRUE) %>%
  st_make_valid() %>%
  select(ADM1_EN, ADM1_PCODE, ADM2_EN, ADM2_PCODE, ADM3_EN, ADM3_PCODE)
stopifnot(st_crs(woreda) == st_crs(gps))

adm3 <- woreda %>% arrange(ADM3_PCODE) %>% mutate(area_id = row_number())
adm2 <- suppressWarnings(woreda %>% group_by(ADM1_EN, ADM1_PCODE, ADM2_EN, ADM2_PCODE) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")) %>%
  arrange(ADM2_PCODE) %>% mutate(area_id = row_number())
adm1 <- suppressWarnings(woreda %>% group_by(ADM1_EN, ADM1_PCODE) %>%
  summarise(geometry = st_union(geometry), .groups = "drop"))

cluster_loc <- suppressWarnings(st_join(gps, woreda, join = st_within))
if (anyNA(cluster_loc$ADM3_PCODE)) {
  bad <- cluster_loc$cluster[is.na(cluster_loc$ADM3_PCODE)]
  stop("clusters not inside any woreda polygon: ", paste(bad, collapse = ", "))
}
stopifnot(!anyDuplicated(cluster_loc$cluster))

missing_gps <- setdiff(unique(children$cluster), cluster_loc$cluster)
if (length(missing_gps)) stop("clusters without GPS record: ", paste(missing_gps, collapse = ", "))

children <- children %>%
  left_join(st_drop_geometry(cluster_loc) %>%
              select(cluster, survey_region, ADM1_EN, ADM2_EN, ADM2_PCODE, ADM3_EN, ADM3_PCODE),
            by = "cluster")

# --- cluster-level counts ---------------------------------------------------------
clusters <- children %>%
  group_by(cluster) %>%
  summarise(urban = first(urban), dhs_stratum = first(dhs_stratum), n_children = n(),
            n_resp = sum(!is.na(value)), n_vacc = sum(value, na.rm = TRUE)) %>%
  left_join(cluster_loc %>% select(cluster, URBAN_RURA, survey_region, LATNUM, LONGNUM,
                                   ADM1_EN, ADM2_EN, ADM2_PCODE, ADM3_EN, ADM3_PCODE),
            by = "cluster") %>%
  st_as_sf()
stopifnot(all((clusters$urban == 1) == (clusters$URBAN_RURA == "U")))

# --- geographic crosswalk (reviewer 2 comment 4) ------------------------------------
# every cluster of the 2019 survey, its survey region (11 regions in place at the
# time) and the current (CSA 2021) units it falls in
crosswalk <- cluster_loc %>% st_drop_geometry() %>%
  transmute(cluster, urban_rural = URBAN_RURA, survey_region,
            current_region = ADM1_EN, current_admin2 = ADM2_EN, ADM2_PCODE, current_admin3 = ADM3_EN, ADM3_PCODE) %>%
  arrange(cluster)
save_csv(crosswalk, "tables", "cluster_region_crosswalk.csv")
crosswalk_counts <- crosswalk %>% count(survey_region, current_region, name = "clusters")
save_csv(crosswalk_counts, "tables", "region_crosswalk_counts.csv")
# current region -> survey region (used by the stratum-effect sensitivity, stage 12)
survey_region_lookup <- crosswalk %>% distinct(current_region, survey_region)
stopifnot(!anyDuplicated(survey_region_lookup$current_region))
save_csv(survey_region_lookup, "data", "survey_region_lookup.csv")
stratum_lookup <- clusters %>% st_drop_geometry() %>% distinct(dhs_stratum, survey_region, urban) %>% arrange(dhs_stratum)
stopifnot(!anyDuplicated(stratum_lookup$dhs_stratum))
save_csv(stratum_lookup, "data", "stratum_lookup.csv")

# --- sampling-weight spread within the 21 strata (reviewer 2 comment 2) ----------------
weight_spread <- children %>% filter(!is.na(value)) %>%
  group_by(dhs_stratum) %>%
  summarise(survey_region = first(survey_region), urban = first(urban), n_children = n(),
            n_clusters = n_distinct(cluster), mean_weight = mean(weight),
            cv_weight = sd(weight) / mean(weight), max_over_min = max(weight) / min(weight),
            deff_weighting = 1 + (sd(weight) / mean(weight))^2) %>%
  mutate(across(c(cv_weight, max_over_min, deff_weighting), ~ round(.x, 3)))
overall_cv <- with(children %>% filter(!is.na(value)), sd(weight) / mean(weight))
save_csv(weight_spread, "tables", "weight_spread_by_stratum.csv")

# --- analytic sample flow and ascertainment (reviewer 2 comment 3) --------------------
n_age <- nrow(children); n_known <- sum(!is.na(children$value)); vacc <- children %>% filter(value == 1)
flow <- kv_flow <- tibble(
  step = c("children aged 12-23 months in the extract", "MCV1 status missing (no card and no recall answer)",
           "children with a known MCV1 status (analytic sample)", "clusters with at least one child in the analytic sample",
           "vaccinated (MCV1 = 1)", "  documented on a vaccination card only", "  caregiver recall only",
           "  card and recall", "vaccination card seen by the interviewer (analytic sample)"),
  n = c(n_age, n_age - n_known, n_known, sum(clusters$n_resp > 0), nrow(vacc),
        sum(vacc$mcv1_card & !vacc$mcv1_recall), sum(!vacc$mcv1_card & vacc$mcv1_recall), sum(vacc$mcv1_card & vacc$mcv1_recall),
        sum(children$has_card == "seen" & !is.na(children$value), na.rm = TRUE)))
save_csv(flow, "tables", "flow_analytic_sample.csv")
print(as.data.frame(flow))

# --- per-unit sample summaries ------------------------------------------------------
unit_summary <- function(level) {
  key <- key_col(level); nm <- name_col(level)
  areas <- if (level == "admin2") adm2 else adm3
  resp <- children %>% filter(!is.na(value)) %>%
    group_by(.data[[key]]) %>%
    summarise(n_clusters = n_distinct(cluster), n_resp = n(), n_vacc = sum(value),
              homogeneous = n_distinct(value) == 1)
  any_cluster <- clusters %>% st_drop_geometry() %>% count(.data[[key]], name = "n_clusters_any")
  st_drop_geometry(areas) %>%
    select(all_of(c(key, nm)), ADM1_EN, area_id) %>%
    left_join(any_cluster, by = key) %>%
    left_join(resp, by = key) %>%
    mutate(across(c(n_clusters_any, n_clusters, n_resp, n_vacc), ~ replace_na(.x, 0L)),
           has_data = n_resp > 0,
           small_n = has_data & n_resp <= cfg$small_n)
}
unit_admin2 <- unit_summary("admin2")
unit_admin3 <- unit_summary("admin3")

# --- numbers for the data-description paragraph -------------------------------------
n_resp_total <- sum(!is.na(children$value))
kv <- function(...) { x <- list(...); tibble(quantity = names(x), value = as.character(unlist(x))) }
data_description <- kv(
  "age filter (age_year)"                     = if (is.null(cfg$age_year)) "all" else paste(cfg$age_year, collapse = ","),
  "children in extract"                        = nrow(children),
  "vaccination responses (non-missing MCV1)"   = n_resp_total,
  "vaccinated (n)"                             = sum(children$value == 1, na.rm = TRUE),
  "vaccinated (%)"                             = round(100 * mean(children$value, na.rm = TRUE), 1),
  "clusters with responses"                    = sum(clusters$n_resp > 0),
  "admin2 units (total)"                       = nrow(unit_admin2),
  "admin2 units with no cluster"               = sum(!unit_admin2$has_data),
  "admin2 units with data"                     = sum(unit_admin2$has_data),
  "admin2 units with n <= 5"                   = sum(unit_admin2$small_n),
  "admin2 units homogeneous (all 0 or all 1)"  = sum(unit_admin2$homogeneous, na.rm = TRUE),
  "admin2 max responses per unit"              = max(unit_admin2$n_resp),
  "admin3 units (total)"                       = nrow(unit_admin3),
  "admin3 units with no cluster"               = sum(!unit_admin3$has_data),
  "admin3 units with data"                     = sum(unit_admin3$has_data),
  "admin3 units with n <= 5"                   = sum(unit_admin3$small_n),
  "admin3 units homogeneous (all 0 or all 1)"  = sum(unit_admin3$homogeneous, na.rm = TRUE),
  "admin3 max responses per unit"              = max(unit_admin3$n_resp),
  "CV of sampling weights, overall"            = round(overall_cv, 3),
  "CV of sampling weights, within strata (median)" = round(median(weight_spread$cv_weight), 3),
  "CV of sampling weights, within strata (max)"    = round(max(weight_spread$cv_weight), 3)
)
print(as.data.frame(data_description))

# --- save ---------------------------------------------------------------------------
save_rds(children, "data", "children.rds")
save_rds(clusters, "data", "clusters.rds")
save_rds(adm1, "data", "adm1.rds"); save_rds(adm2, "data", "adm2.rds"); save_rds(adm3, "data", "adm3.rds")
save_csv(unit_admin2, "data", "unit_summary_admin2.csv")
save_csv(unit_admin3, "data", "unit_summary_admin3.csv")
save_csv(data_description, "tables", "data_description.csv")
log_msg("02 data prep: done (%d children, %d clusters)", nrow(children), nrow(clusters))
