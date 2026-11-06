# =============================================================================
# 01_setup.R - packages, helper functions, plotting theme.
# Sourced by every stage after 00_config.R.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(purrr); library(stringr)
  library(tibble); library(forcats)
  library(sf); library(terra); library(exactextractr); library(spdep)
  library(survey); library(INLA)
  library(ggplot2); library(patchwork); library(ggpubr); library(ggnewscale)
  library(ggrepel); library(ggtext); library(viridis); library(scales); library(msm)
})
options(dplyr.summarise.inform = FALSE, survey.lonely.psu = "adjust", width = 120)
sf_use_s2(FALSE)                      # planar operations
set.seed(cfg$seed)

# --- small helpers -------------------------------------------------------------
expit <- function(x) 1 / (1 + exp(-x))
logit <- function(p) log(p / (1 - p))

ensure_dirs <- function() invisible(lapply(cfg$out, dir.create, recursive = TRUE, showWarnings = FALSE))
ensure_dirs()

log_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", sprintf(...), "\n")

out_path <- function(dir, file) file.path(cfg$out[[dir]], file)
save_rds <- function(x, dir, file) { saveRDS(x, out_path(dir, file)); invisible(x) }
read_out <- function(dir, file) readRDS(out_path(dir, file))
save_csv <- function(x, dir, file) { write_csv(x, out_path(dir, file), na = ""); invisible(x) }
read_csv_out <- function(dir, file) read_csv(out_path(dir, file), show_col_types = FALSE)

# Probability that each unit is in the bottom quintile, from a draws matrix (units x draws)
prob_bottom_quintile <- function(M) {
  n <- nrow(M); k <- floor(0.2 * n)
  in_bottom <- apply(M, 2, function(sim) { out <- logical(n); out[order(sim)[seq_len(k)]] <- TRUE; out })
  rowMeans(in_bottom)
}

# Administrative levels ---------------------------------------------------------
LEVELS <- c("admin2", "admin3")
key_col  <- function(level) switch(level, admin1 = "ADM1_EN", admin2 = "ADM2_PCODE", admin3 = "ADM3_PCODE")
name_col <- function(level) switch(level, admin1 = "ADM1_EN", admin2 = "ADM2_EN",    admin3 = "ADM3_EN")
level_label <- function(level) switch(level, admin1 = "first", admin2 = "second", admin3 = "third")
level_ord   <- function(level) switch(level, admin1 = "1st", admin2 = "2nd", admin3 = "3rd")

# --- figure helpers -----------------------------------------------------------
save_fig <- function(p, dir, name, width, height, dpi = cfg$figures$dpi, tiff = FALSE) {
  png_file <- out_path(dir, paste0(name, ".png"))
  ggsave(png_file, p, width = width, height = height, dpi = dpi, bg = "white")
  if (tiff) ggsave(out_path(dir, paste0(name, ".tiff")), p, width = width, height = height,
                   dpi = cfg$figures$tiff_dpi %||% 300, bg = "white", compression = "lzw")
  log_msg("figure -> %s", png_file)
  invisible(png_file)
}

BASE <- cfg$figures$base_size %||% 14      # base font size for all panels

tag_panels <- function(p) p + plot_annotation(tag_levels = "a", tag_suffix = ".") &
  theme(plot.tag = element_text(size = BASE + 6, face = "bold"))

# Colours used across the paper
coverage_pal <- c("#7F0000", "#C0392B", "#E67E22", "#F1C40F", "#A8D5A2", "#27AE60", "#145A32")
COL_MISSING  <- "#495057"    # dark grey: no data
COL_ZERO_SE  <- "#E63946"    # red: near-zero design-based SE
COL_POINT    <- "#D35400"
COL_REF      <- "#3a86ff"

theme_scatter <- function(base_size = BASE) {
  theme_pubclean(base_size = base_size) +
    theme(axis.title = element_text(face = "bold", size = base_size),
          axis.title.x = element_text(margin = margin(t = 8)),
          axis.title.y = element_text(margin = margin(r = 8)),
          axis.text = element_text(colour = "grey30", size = base_size - 2),
          panel.grid.major = element_line(colour = "grey92", linewidth = 0.4),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "white", colour = NA),
          legend.position = "none")
}

theme_map <- function(base_size = BASE) {
  theme_void(base_size = base_size) +
    theme(legend.position = "bottom",
          legend.title = element_text(face = "bold", size = base_size),
          legend.text = element_text(size = base_size - 2),
          plot.background = element_rect(fill = "white", colour = NA))
}

# Horizontal colour bar used by every map legend
map_colourbar <- function(barwidth = 150, barheight = 7)
  guide_colourbar(direction = "horizontal", barheight = unit(barheight, "mm"),
                  barwidth = unit(barwidth, "mm"), title.position = "top", title.hjust = 0.5,
                  frame.colour = "#555555", ticks.colour = "#555555")

# Coverage fill scale shared by the model-based maps (0-100 %)
scale_fill_coverage <- function(name = "The first dose of measles (MCV1) coverage (%)") {
  scale_fill_gradientn(name = name, colours = coverage_pal, limits = c(0, 100),
                       breaks = c(0, 25, 50, 75, 100), labels = c("0", "25", "50", "75", "100%"),
                       na.value = "transparent", guide = map_colourbar())
}
