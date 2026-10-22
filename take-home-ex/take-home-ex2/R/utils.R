# =============================================================================
# R/utils.R  -- shared helpers for the technical report and the analysis scripts
# Take-home Exercise 2 | ISSS626
#
# One typographic and colour system for every graphic in the report, plus the
# caching layer that keeps the permutation-based statistics reproducible
# without recomputing them on every render.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2)
})

# ---- Reproducible caching ---------------------------------------------------
# The conditional permutation tests behind local Moran's I, local Geary's C and
# Gi* are the expensive part of this analysis, and the space-time cube runs
# them once per period. `cache_rds()` evaluates an expression once, stores the
# result, and returns the stored object on later renders. The code that
# produced it stays visible in the report, so caching removes recomputation and
# not transparency. Deleting data/cache/ forces a full recomputation; every
# script sets the same seed, so the recomputed values are identical.
CACHE_DIR <- "data/cache"
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

cache_rds <- function(name, expr, refresh = getOption("th2.refresh", FALSE)) {
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
# Three palettes, each tied to one statistic, so that a colour means the same
# thing everywhere in the report.
PAL <- list(
  conflict = "#9f1239", conflict_l = "#fb7185", conflict_xl = "#ffe4e6",
  cool     = "#1d4ed8", cool_l     = "#93c5fd",
  grey     = "#57534e", grey_l     = "#d6d3d1", grid = "#e7e5e4",
  ink      = "#1c1917", paper      = "#fcfcfb"
)

# Sequential ramp for counts and rates. Every intensity surface in the report
# uses it, so two choropleths are comparable at a glance.
RAMP_SEQ <- c("#fff7f3", "#fde0dd", "#fcc5c0", "#fa9fb5",
              "#f768a1", "#dd3497", "#ae017e", "#7a0177")

# LISA quadrants. High-High and Low-Low carry the saturated hues because they
# are the clusters; the two spatial-outlier quadrants are deliberately paler so
# that a reader does not confuse an outlier with a cluster.
PAL_LISA <- c(
  "High-High"     = "#b2182b",
  "Low-Low"       = "#2166ac",
  "High-Low"      = "#ef8a62",
  "Low-High"      = "#67a9cf",
  "Insignificant" = "#f0efee"
)

# Gi* hot and cold spots, diverging around zero.
PAL_GI <- c(
  "Hot spot"      = "#b2182b",
  "Cold spot"     = "#2166ac",
  "Insignificant" = "#f0efee"
)

# The eight emerging hot spot categories, in Esri's own ordering so that a
# reader who knows the tool recognises the legend. Hot categories run warm,
# cold categories run cool, and anything the test cannot classify stays grey.
PAL_EHSA <- c(
  "New hot spot"           = "#fcae91",
  "Consecutive hot spot"   = "#fb6a4a",
  "Intensifying hot spot"  = "#de2d26",
  "Persistent hot spot"    = "#a50f15",
  "Diminishing hot spot"   = "#fdd0a2",
  "Sporadic hot spot"      = "#fdae6b",
  "Oscillating hot spot"   = "#e6550d",
  "Historical hot spot"    = "#8c2d04",
  "New cold spot"          = "#c6dbef",
  "Consecutive cold spot"  = "#6baed6",
  "Intensifying cold spot" = "#2171b5",
  "Persistent cold spot"   = "#08306b",
  "Diminishing cold spot"  = "#dadaeb",
  "Sporadic cold spot"     = "#9e9ac8",
  "Oscillating cold spot"  = "#6a51a3",
  "Historical cold spot"   = "#3f007d",
  "No pattern detected"    = "#f0efee"
)

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
      panel.grid.major = element_line(linewidth = 0.2, colour = "#f0efee"),
      panel.background = element_rect(fill = PAL$paper, colour = NA),
      panel.border     = element_rect(fill = NA, colour = PAL$grid,
                                      linewidth = 0.4),
      legend.position  = "right",
      legend.key.width = unit(0.32, "cm"),
      legend.key.height = unit(1.0, "cm")
    )
}

# ---- Formatters --------------------------------------------------------------
n_fmt <- function(x, d = 0) formatC(x, format = "f", big.mark = ",", digits = d)

pct <- function(x, d = 1) paste0(formatC(x, format = "f", digits = d), "%")

# Vectorised: several tables format a whole column of p-values at once.
pv <- function(p) {
  vapply(p, function(x) {
    if (is.na(x)) return("NA")
    if (x < 0.001) return("< 0.001")
    formatC(x, format = "f", digits = 3)
  }, character(1), USE.NAMES = FALSE)
}

# Kendall's tau carries a sign that matters, so it is never shown unsigned.
tau_fmt <- function(x, d = 3) {
  sprintf("%s%s", ifelse(x >= 0, "+", "−"),
          formatC(abs(x), format = "f", digits = d))
}

# ---- Significance ------------------------------------------------------------
# Every local statistic in this report is tested at each of 330 townships at
# once. Reporting the raw 0.05 threshold alone would put roughly 17 false
# positives on each map, so both the raw and the false-discovery-rate adjusted
# verdicts are carried through and the maps say which one they are drawing.
sig_flag <- function(p, alpha = 0.05, method = c("none", "fdr", "bonferroni")) {
  method <- match.arg(method)
  if (method == "none") return(p < alpha)
  p.adjust(p, method = method) < alpha
}

# Expected number of false positives at a given threshold, used in the prose so
# that a count of significant townships is never quoted without its benchmark.
expected_fp <- function(n, alpha = 0.05) n * alpha

# ---- LISA classification -----------------------------------------------------
# The quadrant label is built from the sign of the centred value and the sign of
# its centred spatial lag, not from any quadrant field a package happens to
# return: that keeps the label meaning the same regardless of which
# implementation produced the statistic.
lisa_quadrant <- function(x, lag_x, p, alpha = 0.05,
                          method = c("none", "fdr", "bonferroni")) {
  sig <- sig_flag(p, alpha, match.arg(method))
  xc  <- x - mean(x, na.rm = TRUE)
  lc  <- lag_x - mean(lag_x, na.rm = TRUE)
  out <- ifelse(xc > 0 & lc > 0, "High-High",
         ifelse(xc < 0 & lc < 0, "Low-Low",
         ifelse(xc > 0 & lc < 0, "High-Low", "Low-High")))
  out[!sig | is.na(sig)] <- "Insignificant"
  factor(out, levels = names(PAL_LISA))
}

# Hot and cold spots are read off the SIGN OF Gi*, which is what the statistic
# measures: whether the neighbourhood sum around a township is large or small
# relative to the whole study area.
gi_class <- function(gi_star, p, alpha = 0.05,
                     method = c("none", "fdr", "bonferroni")) {
  sig <- sig_flag(p, alpha, match.arg(method))
  out <- ifelse(gi_star > 0, "Hot spot", "Cold spot")
  out[!sig | is.na(sig)] <- "Insignificant"
  factor(out, levels = names(PAL_GI))
}
