# =============================================================================
# 09_figures_main.R - Fig 1 to Fig 5 of the main manuscript.
#   Fig 1 (a) SD with vs without urban, admin-2; (b) same, admin-3;
#         (c) coverage map admin-2 (urban model); (d) coverage map admin-3
#   Fig 2 (a-b) |direct SE - model SD| vs sample size, admin-2 / admin-3
#   Fig 3 (a-b) direct vs model-based coverage, admin-2 / admin-3
#   Fig 4 (a) coverage vs MBG admin-2; (b) SE vs MBG admin-2; (c-d) admin-3
#   Fig 5 (a-b) posterior SD maps; (c-d) probability of bottom quintile   (revision, reviewer 12)
# Outputs: output/figures/main/Fig{1..5}.png (+ .tiff)
# =============================================================================

source("R/plot_functions.R")
log_msg("09 main figures")
fd <- load_fig_data()
stab <- list(admin2 = read_csv_out("estimates", "stability_admin2.csv"),
             admin3 = read_csv_out("estimates", "stability_admin3.csv"))
tiff <- isTRUE(cfg$figures$tiff)

# Fig 1 ------------------------------------------------------------------------------
fig1 <- (panel_sd_urban_vs_nourban(fd, "admin2") | panel_sd_urban_vs_nourban(fd, "admin3")) /
        (map_coverage(fd, "admin2") | map_coverage(fd, "admin3")) +
  plot_layout(heights = c(1, 1.7), guides = "collect") &
  theme(legend.position = "bottom")
save_fig(tag_panels(fig1), "fig_main", "Fig1", width = 16, height = 15, tiff = tiff)

# Fig 2 ------------------------------------------------------------------------------
fig2 <- panel_shrinkage(fd, "admin2") | panel_shrinkage(fd, "admin3")
save_fig(tag_panels(fig2), "fig_main", "Fig2", width = 16, height = 6, tiff = tiff)

# Fig 3 ------------------------------------------------------------------------------
fig3 <- panel_direct_vs_model(fd, "admin2") | panel_direct_vs_model(fd, "admin3")
save_fig(tag_panels(fig3), "fig_main", "Fig3", width = 16, height = 8, tiff = tiff)

# Fig 4 ------------------------------------------------------------------------------
fig4 <- (panel_mbg(fd, "admin2", "coverage") | panel_mbg(fd, "admin2", "se")) /
        (panel_mbg(fd, "admin3", "coverage") | panel_mbg(fd, "admin3", "se"))
save_fig(tag_panels(fig4), "fig_main", "Fig4", width = 16, height = 13, tiff = tiff)

# Fig 5 ------------------------------------------------------------------------------
sd_lim <- range(100 * c(fd$comp$admin2$sae_u_sd, fd$comp$admin3$sae_u_sd))
fig5 <- (map_sd(fd, "admin2", sd_lim) | map_sd(fd, "admin3", sd_lim)) /
        (map_prob_bottom(fd, "admin2", stab$admin2) | map_prob_bottom(fd, "admin3", stab$admin3)) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
save_fig(tag_panels(fig5), "fig_main", "Fig5", width = 16, height = 15, tiff = tiff)

log_msg("09 main figures: done")
