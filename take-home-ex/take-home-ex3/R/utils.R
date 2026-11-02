# =============================================================================
# R/utils.R  -- shared helpers for the technical report and the analysis scripts
# Take-home Exercise 3 | ISSS626
#
# One typographic and colour system for every graphic in the report, plus the
# caching layer that keeps the bootstrap and cross-validation runs reproducible
# without recomputing them on every render.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2)
})

# ---- Reproducible caching ---------------------------------------------------
# Repeated calibration is the expensive part of this analysis. `cache_rds()`
# evaluates an expression once, stores the result, and returns the stored
# object on later runs. Deleting data/cache/ forces a full recomputation; every
# script sets the same seed, so the recomputed values are identical.
CACHE_DIR <- "data/cache"
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

cache_rds <- function(name, expr, refresh = getOption("th3.refresh", FALSE)) {
  f <- file.path(CACHE_DIR, paste0(name, ".rds"))
  if (file.exists(f) && !refresh) return(readRDS(f))
  t0  <- Sys.time()
  val <- force(expr)
  saveRDS(val, f)
  message(sprintf("[cache] computed %-30s %.1fs", name,
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  val
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- Colour ------------------------------------------------------------------
PAL <- list(
  rail   = "#0f766e", rail_l = "#99f6e4",
  new    = "#c2410c", new_l  = "#fed7aa",
  lose   = "#b2182b", gain   = "#2166ac",
  grey   = "#57534e", grey_l = "#d6d3d1", grid = "#e7e5e4",
  ink    = "#1c1917", paper  = "#fcfcfb"
)

# Official line colours, so that a reader who knows the network recognises them.
PAL_LINE <- c(NS = "#d42e12", EW = "#009645", CG = "#009645", NE = "#9900aa",
              CC = "#fa9e0d", CE = "#fa9e0d", DT = "#005ec4", TE = "#9d5b25",
              BP = "#748477", SE = "#748477", SW = "#748477", PE = "#748477",
              PW = "#748477", STC = "#748477", PTC = "#748477",
              JS = "#0099aa", JE = "#0099aa", JW = "#0099aa", CR = "#97c616")

# Sequential ramp for demand and ridership; diverging ramp for change.
RAMP_SEQ <- c("#f7fcf5", "#d9f0d3", "#a6dba0", "#5aae61", "#1b7837", "#00441b")
RAMP_DIV <- c("#b2182b", "#ef8a62", "#fddbc7", "#f7f7f7", "#d1e5f0", "#67a9cf", "#2166ac")

# ---- Presentation -----------------------------------------------------------
theme_th <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = rel(1.05),
                                      margin = margin(b = 2)),
      plot.subtitle    = element_text(colour = PAL$grey, size = rel(0.88),
                                      margin = margin(b = 9)),
      plot.caption     = element_text(colour = PAL$grey, size = rel(0.74),
                                      hjust = 0, margin = margin(t = 9)),
      axis.title       = element_text(size = rel(0.86)),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.28, colour = PAL$grid),
      strip.text       = element_text(face = "bold", size = rel(0.85)),
      legend.position  = "bottom",
      legend.key.height = unit(0.4, "cm"),
      legend.title     = element_text(size = rel(0.82)),
      legend.text      = element_text(size = rel(0.78)),
      plot.title.position   = "plot",
      plot.caption.position = "plot"
    )
}

theme_map <- function(base_size = 11) {
  theme_th(base_size) +
    theme(
      axis.text        = element_blank(),
      axis.title       = element_blank(),
      axis.ticks       = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_rect(fill = PAL$paper, colour = NA),
      panel.border     = element_rect(fill = NA, colour = PAL$grid, linewidth = 0.4),
      legend.position  = "right",
      legend.key.width = unit(0.32, "cm"),
      legend.key.height = unit(1.0, "cm")
    )
}

# ---- Formatters --------------------------------------------------------------
n_fmt <- function(x, d = 0) formatC(x, format = "f", big.mark = ",", digits = d)
pct   <- function(x, d = 1) paste0(formatC(x, format = "f", digits = d), "%")

pv <- function(p) {
  vapply(p, function(x) {
    if (is.na(x)) return("NA")
    if (x < 0.001) return("< 0.001")
    formatC(x, format = "f", digits = 3)
  }, character(1), USE.NAMES = FALSE)
}
