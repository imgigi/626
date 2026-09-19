# =============================================================================
# R/utils.R  -- shared helpers for the technical report and the analysis scripts
# Take-home Exercise 1 | ISSS626
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2)
})

# ---- Reproducible caching ---------------------------------------------------
# Monte-Carlo envelopes and point-process fits are the expensive part of this
# analysis. `cache_rds()` evaluates an expression once, stores the result, and
# returns the stored object on subsequent renders. The *code that ran remains
# visible in the report*, so the cache only removes recomputation, not
# transparency. Delete data/cache/ to force a full recomputation.
CACHE_DIR <- "data/cache"
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

cache_rds <- function(name, expr, refresh = getOption("th1.refresh", FALSE)) {
  f <- file.path(CACHE_DIR, paste0(name, ".rds"))
  if (file.exists(f) && !refresh) return(readRDS(f))
  t0  <- Sys.time()
  val <- force(expr)
  saveRDS(val, f)
  message(sprintf("[cache] computed %-28s %.1fs", name,
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  val
}

# ---- Presentation -----------------------------------------------------------
# One typographic system for every statistical graphic, so the report reads as
# a single document rather than a pile of default plots.
PAL <- list(
  fire   = "#b91c1c", fire_l = "#f97316", fire_xl = "#fed7aa",
  water  = "#0e7490", grey   = "#57534e", grey_l = "#d6d3d1",
  env    = "#cbd5e1", theo   = "#1e293b", ok = "#15803d"
)

theme_th <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title      = element_text(face = "bold", size = rel(1.05),
                                     margin = margin(b = 2)),
      plot.subtitle   = element_text(colour = PAL$grey, size = rel(0.88),
                                     margin = margin(b = 9)),
      plot.caption    = element_text(colour = PAL$grey, size = rel(0.74),
                                     hjust = 0, margin = margin(t = 9)),
      axis.title      = element_text(size = rel(0.86)),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.28, colour = "#e7e5e4"),
      strip.text      = element_text(face = "bold", size = rel(0.85)),
      legend.position = "bottom",
      legend.key.height = unit(0.4, "cm"),
      legend.title    = element_text(size = rel(0.82)),
      legend.text     = element_text(size = rel(0.78)),
      plot.title.position = "plot",
      plot.caption.position = "plot"
    )
}

# Envelope plots recur ~10 times; one function keeps them identical and keeps
# the axis labelling honest (distance in km, function value in km^2 for K).
plot_env <- function(e, title, subtitle = NULL, caption = NULL,
                     xlab = "Distance r (km)", ylab = NULL,
                     scale_x = 1000, scale_y = 1, ylab_default = "Summary function",
                     obs_lab = "Observed", theo_lab = "Null expectation") {
  d <- as.data.frame(e)
  d$r <- d$r / scale_x
  yy  <- c("obs", "mmean", "theo", "lo", "hi")
  for (v in intersect(yy, names(d))) d[[v]] <- d[[v]] / scale_y
  ref <- if ("mmean" %in% names(d)) "mmean" else "theo"
  ggplot(d, aes(r)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = "env"), alpha = 0.85) +
    geom_line(aes(y = .data[[ref]], colour = "theo"), linewidth = 0.5,
              linetype = "22") +
    geom_line(aes(y = obs, colour = "obs"), linewidth = 0.85) +
    scale_fill_manual(NULL, values = c(env = PAL$env),
                      labels = "Global envelope (α = 0.01)") +
    scale_colour_manual(NULL, values = c(obs = PAL$fire, theo = PAL$theo),
                        labels = c(obs = obs_lab, theo = theo_lab),
                        breaks = c("obs", "theo")) +
    labs(title = title, subtitle = subtitle, caption = caption,
         x = xlab, y = ylab %||% ylab_default) +
    theme_th()
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Envelope plots for L-type functions are far easier to read when CENTRED on
# the reference curve: deviations then appear as departures from a horizontal
# zero line instead of as two nearly-parallel rising curves. The y axis becomes
# "excess L", in metres, which for an L function has a direct reading -- the
# extra radius over which the observed pattern packs the same number of
# neighbours as the reference process.
env_centred <- function(e, scale_x = 1000) {
  d <- as.data.frame(e)
  ref <- if ("mmean" %in% names(d)) d$mmean else d$theo
  data.frame(r = d$r / scale_x, ref = ref,
             obs = d$obs - ref, lo = d$lo - ref, hi = d$hi - ref)
}

plot_env_centred <- function(e, title, subtitle = NULL, caption = NULL,
                             ylab = "Observed L(r) − reference (m)",
                             xlab = "Distance r (km)",
                             ref_lab = "Reference process",
                             obs_lab = "Observed",
                             band_lab = "Global envelope (α = 0.01)",
                             annotate_exceed = TRUE) {
  d <- env_centred(e)
  out <- which(d$obs > d$hi | d$obs < d$lo)
  p <- ggplot(d, aes(r)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = "band"), alpha = 0.9) +
    geom_hline(yintercept = 0, linetype = "22", colour = PAL$theo,
               linewidth = 0.45) +
    geom_line(aes(y = obs, colour = "obs"), linewidth = 0.9) +
    scale_fill_manual(NULL, values = c(band = PAL$env), labels = band_lab) +
    scale_colour_manual(NULL, values = c(obs = PAL$fire), labels = obs_lab) +
    labs(title = title, subtitle = subtitle, caption = caption,
         x = xlab, y = ylab) +
    theme_th()
  if (annotate_exceed && length(out)) {
    p <- p + annotate("rect", xmin = d$r[min(out)], xmax = d$r[max(out)],
                      ymin = -Inf, ymax = Inf, fill = PAL$fire, alpha = 0.05)
  }
  p
}

# How far outside the envelope the observed curve strays, as a compact record
# that can be tabulated across bandwidths, units and windows.
env_summary <- function(e) {
  d <- env_centred(e)
  out <- which(d$obs > d$hi | d$obs < d$lo)
  above <- which(d$obs > d$hi); below <- which(d$obs < d$lo)
  list(n_out = length(out),
       pct_out = 100 * length(out) / nrow(d),
       r_first = if (length(out)) d$r[min(out)] else NA_real_,
       r_last  = if (length(out)) d$r[max(out)] else NA_real_,
       direction = if (length(above) > length(below)) "above (clustered)"
                   else if (length(below) > 0) "below (inhibited)" else "inside",
       max_excess = if (length(out)) max(d$obs[out] - d$hi[out]) else NA_real_)
}

# ---- Cartography ------------------------------------------------------------
# All maps in this report are drawn with ggplot2 + ggspatial rather than a
# mixture of plotting systems, so that typography, colour and legend treatment
# are identical across statistical graphics and maps.
theme_map <- function(base_size = 11) {
  theme_th(base_size) +
    theme(
      axis.text        = element_blank(),
      axis.title       = element_blank(),
      axis.ticks       = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2, colour = "#f0efee"),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "#fcfcfb", colour = NA),
      panel.border     = element_rect(fill = NA, colour = "#e7e5e4",
                                      linewidth = 0.4),
      legend.position  = "right",
      legend.key.width = unit(0.32, "cm"),
      legend.key.height= unit(1.0, "cm")
    )
}

# A spatstat `im` carries its values in a [ny, nx] matrix whose first row is the
# LOWEST y. as.data.frame.im already resolves that correctly, so it is used
# rather than hand-rolled matrix flipping; NA cells (outside the window) are
# dropped so geom_raster leaves them blank.
im_df <- function(im, value_name = "value", scale = 1) {
  d <- as.data.frame(im)
  names(d) <- c("x", "y", value_name)
  d[[value_name]] <- d[[value_name]] * scale
  d[!is.na(d[[value_name]]), ]
}

# Sequential ramp used for every intensity surface, so that two density maps in
# this report are always comparable at a glance.
scale_fill_fire <- function(name = expression(paste("per ", km^2)), ...) {
  scale_fill_gradientn(
    name   = name,
    colours = c("#fffbf5", "#fee9c8", "#fdbb6f", "#f97c31",
                "#dc4a1f", "#a81c0c", "#6b0d05"),
    na.value = NA, ...)
}

# Compact number formatting used throughout the prose.
n_fmt <- function(x, d = 0) formatC(x, format = "f", big.mark = ",", digits = d)
km    <- function(x, d = 1) paste0(formatC(x / 1000, format = "f", digits = d), " km")
# Vectorised, because several tables format a whole column of p-values at once.
pv <- function(p) {
  vapply(p, function(x) {
    if (is.na(x)) return("NA")
    if (x < 0.001) return("< 0.001")
    formatC(x, format = "f", digits = 3)
  }, character(1), USE.NAMES = FALSE)
}

# spatstat's envelope-based tests return `statistic` as a one-row DATA FRAME
# (the statistic plus its rank among the simulations), not the named numeric an
# htest normally carries. This extracts the statistic itself in either case.
stat_of <- function(test) {
  st <- test$statistic
  if (is.data.frame(st)) return(unname(st[[1]]))
  unname(st[1])
}

# A single place to describe a Monte-Carlo test result in words, so that no
# p-value in the report is reported without its null model.
mc_line <- function(test, null_label) {
  sprintf("%s = %s, p %s (%d simulations of %s)",
          if (grepl("Diggle", test$method)) "DCLF u" else "MAD",
          formatC(unname(test$statistic), format = "g", digits = 4),
          pv(test$p.value),
          test$parameter[["nsim"]] %||% NA_integer_, null_label)
}
