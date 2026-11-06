# =============================================================================
# 10_figures_supp.R - supplementary figures.
#   S2 Appendix (model-based):
#     Fig B1 (a) admin-3, (b) admin-2: posterior SD per unit, with / without urban
#     Fig B2 (a) admin-3, (b) admin-2: direct SE vs model SD, both specifications
#     Fig B3 (a) admin-3, (b) admin-2: coverage map, model without urban
#     Fig B4 (a) admin-3, (b) admin-2: units with >= 70 % probability of bottom quintile
#   S3 Appendix (design-based):
#     Fig C1  admin-1 direct estimates with the WHO 90 % line
#     Fig C2 (a) admin-3 samples map, (b) admin-2 samples map,
#            (c) admin-2 n vs SE, (d) admin-3 n vs SE,
#            (e) admin-3 direct coverage map, (f) admin-2 direct coverage map
# Panels are also saved individually (FigB1a.png, ...) for assembling in Word.
# =============================================================================

source("R/plot_functions.R")
log_msg("10 supplementary figures")
fd <- load_fig_data()
stab <- list(admin2 = read_csv_out("estimates", "stability_admin2.csv"),
             admin3 = read_csv_out("estimates", "stability_admin3.csv"))
direct_admin1 <- read_csv_out("estimates", "direct_admin1.csv")

save_pair <- function(pa, pb, dir, name, width, height, widths = c(1, 1)) {
  save_fig(pa, dir, paste0(name, "a"), width = width / 2, height = height)
  save_fig(pb, dir, paste0(name, "b"), width = width / 2, height = height)
  save_fig(tag_panels((pa | pb) + plot_layout(widths = widths)), dir, name, width = width, height = height)
}

# --- S2 ------------------------------------------------------------------------------
save_pair(panel_sd_by_unit(fd, "admin3"), panel_sd_by_unit(fd, "admin2"), "fig_s2", "FigB1", 16, 6)
save_pair(panel_direct_se_vs_model_sd(fd, "admin3"), panel_direct_se_vs_model_sd(fd, "admin2"), "fig_s2", "FigB2", 16, 6)

b3 <- (map_coverage(fd, "admin3", "nourban") | map_coverage(fd, "admin2", "nourban")) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
save_fig(map_coverage(fd, "admin3", "nourban"), "fig_s2", "FigB3a", 9, 8)
save_fig(map_coverage(fd, "admin2", "nourban"), "fig_s2", "FigB3b", 9, 8)
save_fig(tag_panels(b3), "fig_s2", "FigB3", 16, 8)

b4a <- map_highlight(fd, "admin3", stab$admin3); b4b <- map_highlight(fd, "admin2", stab$admin2)
save_fig(b4a, "fig_s2", "FigB4a", 12, 9); save_fig(b4b, "fig_s2", "FigB4b", 12, 9)
save_fig(tag_panels(b4a | b4b), "fig_s2", "FigB4", 22, 9)

# revision (reviewer 2): agreement beyond correlation, posterior predictive check
save_pair(panel_bland_altman(fd, "admin3"), panel_bland_altman(fd, "admin2"), "fig_s2", "FigB5", 16, 6.5)
if (file.exists(out_path("data", "pit.rds"))) {
  pit_df <- read_out("data", "pit.rds")
  save_pair(panel_pit(pit_df, "admin3"), panel_pit(pit_df, "admin2"), "fig_s2", "FigB6", 16, 6)
}

# --- S3 ------------------------------------------------------------------------------
save_fig(panel_admin1_forest(direct_admin1), "fig_s3", "FigC1", 8, 5)

c2 <- list(a = map_samples(fd, "admin3"), b = map_samples(fd, "admin2"),
           c = panel_n_vs_se(fd, "admin2"),  d = panel_n_vs_se(fd, "admin3"),
           e = map_direct(fd, "admin3"),     f = map_direct(fd, "admin2"))
for (k in names(c2)) save_fig(c2[[k]], "fig_s3", paste0("FigC2", k), width = 8, height = if (k %in% c("c", "d")) 5.5 else 7)
c2_all <- (c2$a | c2$b) / (c2$c | c2$d) / (c2$e | c2$f) + plot_layout(heights = c(1.2, 1, 1.2))
save_fig(tag_panels(c2_all), "fig_s3", "FigC2", width = 16, height = 20, dpi = 250)

log_msg("10 supplementary figures: done")
