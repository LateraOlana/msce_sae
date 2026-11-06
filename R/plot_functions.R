# =============================================================================
# plot_functions.R - panel builders shared by 09_figures_main.R and
# 10_figures_supp.R. Every function returns a ggplot; the figure scripts
# compose them with patchwork and save them.
# =============================================================================

# --- data access -------------------------------------------------------------------
load_fig_data <- function() {
  list(comp = list(admin2 = read_csv_out("estimates", "comparison_admin2.csv"),
                   admin3 = read_csv_out("estimates", "comparison_admin3.csv")),
       adm  = list(admin1 = read_out("data", "adm1.rds"), admin2 = read_out("data", "adm2.rds"),
                   admin3 = read_out("data", "adm3.rds")),
       clusters = read_out("data", "clusters.rds"))
}

# join a per-unit table (with a `code` column) onto the polygons of a level
geo_join <- function(fd, level, df) {
  fd$adm[[level]] %>% select(code = all_of(key_col(level)), geometry) %>% left_join(df, by = "code")
}

adm1_outline <- function(fd, colour = "#1A1A1A", linewidth = 0.9)
  geom_sf(data = fd$adm$admin1, fill = NA, colour = colour, linewidth = linewidth)

r_label <- function(size = 7, x = 0.05, y = 0.95)
  stat_cor(aes(label = after_stat(r.label)), method = "pearson", label.x.npc = x, label.y.npc = y,
           colour = "#34495E", size = size)

# text block with Pearson r, Spearman rho and RMSE (percentage points) for x vs y (both 0-100)
fit_stats_label <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  sprintf("r = %.2f\nrho = %.2f\nRMSE = %.1f pp", cor(x[ok], y[ok]), cor(x[ok], y[ok], method = "spearman"),
          sqrt(mean((x[ok] - y[ok])^2)))
}

# --- scatter panels -----------------------------------------------------------------
# Fig 1a-b: posterior SD with vs without the urban covariate
panel_sd_urban_vs_nourban <- function(fd, level) {
  d <- fd$comp[[level]] %>% transmute(x = 100 * sae_u_sd, y = 100 * sae_nu_sd)
  lim <- c(0, max(c(d$x, d$y)) * 1.05)
  ggplot(d, aes(x, y)) +
    geom_hex(bins = 40, show.legend = FALSE) +
    scale_fill_gradient(low = "#D6EAF8", high = "#2E86C1") +
    geom_point(alpha = 0.4, colour = "#2E86C1", size = 3) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "#E74C3C", se = FALSE,
                linewidth = 1.2, linetype = "dotdash") +
    geom_abline(intercept = 0, slope = 1, colour = "#27AE60", linewidth = 1, linetype = "longdash") +
    r_label(size = 6, x = 0.12, y = 0.93) +
    annotate("text", x = lim[2] * 0.08, y = lim[2] * 0.80, label = "Linear trend", colour = "#E74C3C",
             size = 5, fontface = "bold", hjust = 0) +
    annotate("text", x = lim[2] * 0.08, y = lim[2] * 0.72, label = "Identity line", colour = "#27AE60",
             size = 5, fontface = "bold", hjust = 0) +
    coord_equal(xlim = lim, ylim = lim) +
    labs(x = "Posterior SD of coverage (pp),\nmodel with urban-rural covariate",
         y = "Posterior SD of coverage (pp),\nmodel without covariate") +
    theme_pubclean(base_size = BASE) +
    theme(axis.title = element_text(size = BASE - 1), panel.grid.major = element_line(colour = "#EAEDED"),
          axis.line = element_line(colour = "#2C3E50"), plot.background = element_rect(fill = "white", colour = NA))
}

# Fig 2: |direct SE - model SD| against sample size
panel_shrinkage <- function(fd, level, max_n = cfg$figures$shrinkage_max_n) {
  d <- fd$comp[[level]] %>% filter(has_data, !is.na(direct_se)) %>%
    transmute(n = n_resp, diff = 100 * abs_se_diff_urban)
  if (!is.null(max_n)) d <- d %>% filter(n <= max_n)
  ggplot(d, aes(n, diff)) +
    geom_hline(yintercept = 0, colour = "#DDDDDD", linetype = "dashed", linewidth = 0.45) +
    geom_point(aes(size = diff), colour = "#C0392B", fill = "white", shape = 21, stroke = 1.1, alpha = 0.85) +
    scale_size_continuous(range = c(1.5, 8), guide = "none") +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, colour = "#C0392B", linewidth = 2.2) +
    scale_x_continuous(name = "Sample size", breaks = breaks_extended(7), expand = expansion(mult = c(0.02, 0.04))) +
    scale_y_continuous(name = "|direct SE - model SD|  (pp)",
                       limits = c(0, max(d$diff) * 1.05), expand = expansion(mult = c(0, 0.03))) +
    theme(plot.background = element_rect(fill = "white", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA),
          panel.grid.major.y = element_line(colour = "#EFEFEF", linewidth = 0.5),
          panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "#CCCCCC", linewidth = 0.55),
          axis.ticks = element_line(colour = "#CCCCCC"),
          axis.title = element_text(face = "bold", size = BASE, colour = "#2C2C2C"),
          axis.title.x = element_text(margin = margin(t = 10)), axis.title.y = element_text(margin = margin(r = 10)),
          axis.text = element_text(size = BASE - 2, colour = "#666666"), legend.position = "none",
          plot.margin = margin(20, 24, 16, 18))
}

# Fig 3: direct vs model-based coverage
panel_direct_vs_model <- function(fd, level) {
  d <- fd$comp[[level]] %>% filter(has_data) %>%
    transmute(x = 100 * direct_est, y = 100 * sae_u_median, s = 100 * sae_u_sd)
  ggplot(d, aes(x, y)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = COL_REF, linewidth = 1) +
    geom_point(aes(size = s), colour = COL_POINT, shape = 21, alpha = 0.85) +
    scale_size_continuous(range = c(1.5, 7)) +
    annotate("text", x = 2, y = 98, label = fit_stats_label(d$x, d$y), hjust = 0, vjust = 1,
             size = 6, colour = "#34495E", lineheight = 0.9) +
    coord_equal(xlim = c(0, 100), ylim = c(0, 100), expand = FALSE) +
    labs(x = "Direct (design-based) coverage estimate (%)", y = "Model-based coverage estimate (%)") +
    theme_scatter() + theme(plot.margin = margin(8, 14, 8, 8)) +
    annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf, colour = "grey40", linewidth = 0.8) +
    annotate("segment", y = -Inf, yend = Inf, x = -Inf, xend = -Inf, colour = "grey40", linewidth = 0.8)
}

# Fig 4: model-based vs MBG, coverage and SE
panel_mbg <- function(fd, level, what = c("coverage", "se")) {
  what <- match.arg(what)
  cc <- fd$comp[[level]]
  d <- if (what == "coverage") cc %>% transmute(x = 100 * sae_u_median, y = 100 * mbg_mean, s = 100 * sae_u_se_ci)
       else cc %>% transmute(x = 100 * sae_u_se_ci, y = 100 * mbg_se_ci, s = 1)
  labs_xy <- if (what == "coverage")
    labs(x = "Small-area model coverage estimate\n(urban-rural covariate) (%)",
         y = "MBG coverage estimate (%)")
  else labs(x = "Small-area model standard error\n(urban-rural covariate) (pp)",
            y = "MBG standard error (pp)")
  ggplot(d, aes(x, y)) +
    { if (what == "coverage") geom_point(aes(size = s), colour = COL_POINT, shape = 21, alpha = 0.85)
      else geom_point(colour = COL_POINT, shape = 21, alpha = 0.85, size = 2) } +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = COL_REF, linewidth = 1) +
    r_label(size = 8) + labs_xy + theme_scatter() +
    annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf, colour = "grey40", linewidth = 0.8) +
    annotate("segment", y = -Inf, yend = Inf, x = -Inf, xend = -Inf, colour = "grey40", linewidth = 0.8)
}

# --- maps --------------------------------------------------------------------------------
# Fig 1c-d, Fig B3: model-based coverage map
map_coverage <- function(fd, level, model = c("urban", "nourban"),
                         legend_name = "The first dose of measles (MCV1) coverage (%)") {
  model <- match.arg(model)
  col <- if (model == "urban") "sae_u_median" else "sae_nu_median"
  d <- geo_join(fd, level, fd$comp[[level]] %>% transmute(code, cov = 100 * .data[[col]]))
  ggplot() +
    geom_sf(data = d, aes(fill = cov), colour = "white", linewidth = 0.04) +
    scale_fill_coverage(legend_name) +
    adm1_outline(fd) +
    theme_map() + theme(plot.margin = margin(6, 6, 6, 6))
}

# Fig 5 a-b: posterior SD of coverage (urban model). `limits` shared across levels
# so patchwork can merge the two legends into one.
map_sd <- function(fd, level, limits = NULL) {
  d <- geo_join(fd, level, fd$comp[[level]] %>% transmute(code, sd = 100 * sae_u_sd))
  if (is.null(limits)) limits <- range(d$sd, na.rm = TRUE)
  ggplot() +
    geom_sf(data = d, aes(fill = sd), colour = "white", linewidth = 0.04) +
    scale_fill_viridis_c(option = "magma", direction = -1, limits = limits,
                         name = "Posterior SD of MCV1 coverage (pp)", guide = map_colourbar()) +
    adm1_outline(fd) +
    theme_map() + theme(plot.margin = margin(6, 6, 6, 6))
}

# Fig 5 c-d: probability of being in the bottom coverage quintile
map_prob_bottom <- function(fd, level, stab) {
  d <- geo_join(fd, level, stab %>% transmute(code, p = prob_bottom))
  ggplot() +
    geom_sf(data = d, aes(fill = p), colour = "white", linewidth = 0.04) +
    scale_fill_viridis_c(option = "rocket", direction = -1, limits = c(0, 1),
                         breaks = c(0, 0.25, 0.5, 0.7, 1), labels = c("0", "0.25", "0.5", "0.7", "1"),
                         name = "P(unit in bottom coverage quintile)", guide = map_colourbar()) +
    adm1_outline(fd) +
    theme_map() + theme(plot.margin = margin(6, 6, 6, 6))
}

# Fig C2 a-b: sample size per unit with cluster locations
map_samples <- function(fd, level) {
  d <- geo_join(fd, level, fd$comp[[level]] %>% transmute(code, n = ifelse(has_data, n_resp, NA_real_)))
  ggplot() +
    geom_sf(data = d, aes(fill = n), colour = "black", linewidth = 0.15) +
    adm1_outline(fd, colour = "gray34", linewidth = 1.1) +
    geom_sf(data = fd$clusters, colour = "red", alpha = 0.65, size = 1.5) +
    scale_fill_viridis_c(option = "mako", name = "Responses per unit", direction = -1, na.value = "#f8f9fa",
                         breaks = pretty_breaks(n = 6), guide = map_colourbar(100, 6)) +
    theme_map() + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
}

# Fig C2 e-f: design-based coverage map (red = near-zero SE, dark grey = no data)
map_direct <- function(fd, level) {
  d <- geo_join(fd, level, fd$comp[[level]] %>%
                  transmute(code, has_data, near_zero_se, cov = 100 * direct_est))
  ggplot() +
    geom_sf(data = d %>% filter(has_data, !near_zero_se), aes(fill = cov), colour = "white", linewidth = 0.05) +
    scale_fill_viridis(name = "Direct MCV1 coverage estimate (%)", option = "mako", begin = 0.15, end = 0.85,
                       direction = -1, limits = c(0, 100), na.value = "transparent", guide = map_colourbar(120, 6)) +
    new_scale_fill() +
    geom_sf(data = d %>% filter(has_data, near_zero_se), fill = COL_ZERO_SE, colour = "white", linewidth = 0.05) +
    geom_sf(data = d %>% filter(!has_data), fill = COL_MISSING, colour = "white", linewidth = 0.05) +
    adm1_outline(fd, colour = "black", linewidth = 1.1) +
    annotate("text", x = Inf, y = Inf, label = c("No data", "Near-zero standard error"),
             colour = c(COL_MISSING, COL_ZERO_SE), hjust = c(1.05, 1.05), vjust = c(1.5, 3), size = 6, fontface = "bold") +
    theme_map()
}

# Fig B4: consistently under-performing units (prob of bottom quintile >= threshold)
map_highlight <- function(fd, level, stab) {
  unit_word <- if (level == "admin2") "zones" else "districts"
  d <- geo_join(fd, level, stab %>% transmute(code, name, ADM1_EN, highlight, prob_bottom))
  region_stats <- d %>% st_drop_geometry() %>% group_by(ADM1_EN) %>%
    summarise(total = n(), high = sum(highlight), .groups = "drop") %>%
    filter(high > 0) %>% arrange(desc(high))
  lab <- d %>% filter(highlight) %>% st_point_on_surface() %>%
    mutate(X = st_coordinates(.)[, 1], Y = st_coordinates(.)[, 2]) %>% st_drop_geometry()
  summary_txt <- paste0("REGIONAL RISK SUMMARY\n",
                        paste0(region_stats$ADM1_EN, ": ", region_stats$high, " ", unit_word, collapse = "\n"))
  bb <- st_bbox(d)
  ggplot() +
    geom_sf(data = d, aes(fill = highlight), colour = "grey", linewidth = 0.25, alpha = 0.9) +
    scale_fill_manual(values = c(`FALSE` = "#F5F5F5", `TRUE` = "#1565C0"), guide = "none") +
    geom_text_repel(data = lab, aes(X, Y, label = paste0(name, "\n(", ADM1_EN, ")")), size = 3.4,
                    fontface = "bold", colour = "#1A1A1A", bg.color = alpha("white", 0.95), bg.r = 0.12,
                    box.padding = 0.5, point.padding = 0.4, segment.colour = "#1565C0", segment.alpha = 0.7,
                    segment.size = 0.6, max.overlaps = Inf, min.segment.length = 0.2, seed = 42, lineheight = 0.9) +
    annotate("text", x = bb["xmax"], y = bb["ymax"], label = "N ↑", hjust = 1, vjust = 1,
             size = 4, colour = "#999999", fontface = "bold") +
    annotate("label", x = bb["xmax"], y = bb["ymin"], label = summary_txt, hjust = 1, vjust = 0,
             fontface = "bold", size = 3.2, colour = "#1565C0", fill = alpha("white", 0.95),
             label.padding = unit(0.5, "lines"), label.r = unit(0.25, "lines")) +
    coord_sf(expand = FALSE, clip = "off") +
    theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA),
                         plot.margin = margin(20, 20, 20, 20))
}

# --- other supplementary panels --------------------------------------------------------
# Fig B1: posterior SD per unit, with and without the urban covariate
panel_sd_by_unit <- function(fd, level) {
  d <- fd$comp[[level]] %>% transmute(area_id, `Without urban` = 100 * sae_nu_sd, `With urban` = 100 * sae_u_sd) %>%
    pivot_longer(-area_id, names_to = "Model", values_to = "SD")
  means <- d %>% group_by(Model) %>% summarise(mean_SD = mean(SD))
  ggplot(d, aes(area_id, SD, colour = Model)) +
    geom_point(size = 3, shape = 21, alpha = 0.4) +
    geom_hline(data = means, aes(yintercept = mean_SD, colour = Model), linetype = "dashed", linewidth = 0.8) +
    scale_colour_manual(values = c(`Without urban` = "#0072B2", `With urban` = "#D55E00")) +
    labs(x = sprintf("%s administrative unit (index, ordered by P-code)", level_ord(level)),
         y = "Posterior standard deviation (pp)", colour = NULL) +
    theme_minimal(base_size = BASE) +
    theme(legend.position = "top", panel.grid.major.x = element_blank(),
          plot.background = element_rect(fill = "white", colour = NA))
}

# Fig B2: direct SE vs model SD for both model specifications
panel_direct_se_vs_model_sd <- function(fd, level) {
  d <- fd$comp[[level]] %>% filter(has_data) %>%
    transmute(Direct_SE = 100 * direct_se, `Urban_SE` = 100 * sae_u_sd, `NoUrban_SE` = 100 * sae_nu_sd) %>%
    pivot_longer(c(Urban_SE, NoUrban_SE), names_to = "Method", values_to = "Unit_SE")
  ggplot(d, aes(Direct_SE, Unit_SE, colour = Method, shape = Method)) +
    geom_point(alpha = 0.9, size = 3, shape = 21) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = c(Urban_SE = "#E69F00", NoUrban_SE = "#56B4E9"),
                        labels = c(Urban_SE = "Urban Unit-Level", NoUrban_SE = "Non-Urban Unit-Level")) +
    labs(x = "Direct estimate standard error (pp)", y = "Model posterior standard deviation (pp)") +
    theme_pubclean(base_size = BASE) +
    theme(legend.position = "top", legend.title = element_text(face = "bold"),
          panel.grid.major = element_line(colour = "#EAEDED"), axis.line = element_line(colour = "#2C3E50"),
          plot.background = element_rect(fill = "white", colour = NA)) +
    scale_x_continuous(expand = expansion(mult = 0.02)) + scale_y_continuous(expand = expansion(mult = 0.02))
}

# Fig B5: Bland-Altman style agreement, model-based (urban) vs direct
panel_bland_altman <- function(fd, level) {
  d <- fd$comp[[level]] %>% filter(has_data) %>%
    transmute(avg = 100 * (direct_est + sae_u_median) / 2, diff = 100 * (sae_u_median - direct_est),
              reliable = factor(ifelse(reliable, "reliable direct estimate", "near-zero design-based SE"),
                                levels = c("reliable direct estimate", "near-zero design-based SE")))
  m <- mean(d$diff); s <- sd(d$diff)
  ggplot(d, aes(avg, diff)) +
    geom_hline(yintercept = 0, colour = "grey60") +
    geom_hline(yintercept = m, colour = COL_POINT, linewidth = 1) +
    geom_hline(yintercept = c(m - 1.96 * s, m + 1.96 * s), colour = COL_POINT, linetype = "dashed") +
    geom_point(aes(shape = reliable), colour = "#2E86C1", size = 2.8, alpha = 0.85) +
    scale_shape_manual(values = c(16, 1), name = NULL) +
    annotate("text", x = 2, y = m, label = sprintf("mean difference %.1f pp", m), hjust = 0, vjust = -0.5,
             colour = COL_POINT, size = 4.5) +
    annotate("text", x = 2, y = m + 1.96 * s, label = sprintf("+1.96 SD = %.1f pp", m + 1.96 * s), hjust = 0, vjust = -0.5,
             colour = COL_POINT, size = 4) +
    annotate("text", x = 2, y = m - 1.96 * s, label = sprintf("-1.96 SD = %.1f pp", m - 1.96 * s), hjust = 0, vjust = 1.5,
             colour = COL_POINT, size = 4) +
    scale_x_continuous(limits = c(0, 100)) +
    labs(x = "Mean of direct and model-based estimate (%)", y = "Model-based minus direct estimate (pp)") +
    theme_scatter() + theme(legend.position = "top")
}

# Fig B6: probability integral transform of the cluster counts (posterior predictive check)
panel_pit <- function(pit_df, level) {
  d <- pit_df %>% filter(level == !!level, !is.na(pit))
  ggplot(d, aes(pit)) +
    geom_histogram(breaks = seq(0, 1, 0.1), fill = "#2E86C1", colour = "white") +
    geom_hline(yintercept = nrow(d) / 10, linetype = "dashed", colour = "grey40") +
    labs(x = "Randomised probability integral transform (PIT) of cluster counts", y = "Clusters") +
    theme_scatter()
}

# Fig C1: admin-1 direct estimates with the WHO 90 % line
panel_admin1_forest <- function(direct_admin1) {
  d <- direct_admin1 %>%
    mutate(est = 100 * est, lo = 100 * ci_l, hi = pmin(100 * ci_u, 100),
           status = ifelse(est > 90, "Above 90%", "Below 90%"),
           ADM1_EN = fct_reorder(ADM1_EN, est))
  ggplot(d, aes(x = est, y = ADM1_EN, colour = status)) +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.2, linewidth = 0.8, alpha = 0.7, orientation = "y") +
    geom_point(size = 3) +
    geom_vline(xintercept = 90, linetype = "dashed", colour = "#e63946", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.0f%%", est)), hjust = -0.4, vjust = -0.1, size = 3, show.legend = FALSE) +
    scale_colour_manual(values = c("Above 90%" = "#2a9d8f", "Below 90%" = "#264653")) +
    scale_x_continuous(breaks = c(seq(0, 100, 20), 90), labels = percent_format(scale = 1)) +
    labs(x = "Estimated MCV1 coverage (%)", y = NULL, colour = "Coverage status") +
    theme_minimal(base_size = BASE) +
    theme(panel.grid.major.y = element_line(colour = "gray90", linewidth = 0.2), panel.grid.minor = element_blank(),
          axis.text = element_text(colour = "gray30"), axis.title.x = element_text(margin = margin(t = 10), colour = "gray30"),
          legend.position = "top", legend.justification = "right", legend.title = element_text(size = 10),
          plot.background = element_rect(fill = "white", colour = NA), plot.margin = unit(c(1, 2, 1, 2), "cm"))
}

# Fig C2 c-d: design-based SE against sample size, split by response homogeneity
panel_n_vs_se <- function(fd, level) {
  d <- fd$comp[[level]] %>% filter(has_data) %>%
    transmute(n = n_resp, se = 100 * direct_se,
              homogeneous = factor(ifelse(homogeneous, "Homogeneous", "Heterogeneous"),
                                   levels = c("Heterogeneous", "Homogeneous")))
  ggplot(d, aes(n, se)) +
    geom_hex(bins = 30) +
    scale_fill_distiller(palette = "Spectral", direction = 1, name = "Count") +
    facet_wrap(~homogeneous, drop = FALSE) +
    labs(x = sprintf("Total samples per %s administrative unit", level_ord(level)),
         y = "Design-based standard error (pp)") +
    theme_bw(base_size = BASE) +
    theme(axis.title = element_text(face = "bold", size = BASE), legend.position = "bottom",
          axis.text.x = element_text(angle = 90, hjust = 1, colour = "gray40"), axis.text.y = element_text(colour = "gray40"),
          panel.grid.major.y = element_blank(), strip.background = element_rect(fill = "#F7F7F7", colour = NA),
          strip.text = element_text(face = "bold", size = 15, colour = "#2A5783"),
          plot.background = element_rect(fill = "white", colour = NA))
}
