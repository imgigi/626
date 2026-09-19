# =============================================================================
# 05b_stkde.R -- spatio-temporal first-order analysis
# Take-home Exercise 1 | ISSS626
#
# Two complementary spatio-temporal intensity estimates:
#
#  (1) DAILY KERNEL SURFACES. One fixed-bandwidth KDE per local day, using the
#      same sigma as the pooled surface so the panels are directly comparable
#      with it and with each other. Units are detections per km2 per day --
#      the most directly interpretable form, and no temporal smoothing is
#      applied, so each panel shows only that day's evidence.
#
#  (2) SPATIO-TEMPORAL KDE (sparr::spattemp.density). A separable space-time
#      kernel that borrows strength across adjacent days. This is the estimator
#      covered in the course material; it trades the literal interpretation of
#      (1) for a less noisy picture of how the field evolves. The CONDITIONAL
#      surfaces z.cond are used, i.e. the spatial density given the day, which
#      answers "given that fire occurred on this day, where was it?"
#
#      tlim is set to c(0.5, 8.5) rather than c(1, 8) so that the tres = 8 grid
#      cell centres land exactly on integer day indices 1..8. With c(1, 8) the
#      grid falls at 1.4375, 2.3125, ... and no panel corresponds to a real day.
# =============================================================================

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(spatstat.geom); library(spatstat.explore)
  library(sparr)
})
source("R/utils.R")
sf_use_s2(FALSE); set.seed(1234)

DER <- "data/derived"
fire_px <- readRDS(file.path(DER, "fire_pxday_sf.rds"))
win     <- readRDS(file.path(DER, "win.rds"))
SIGMA_MAIN <- readRDS(file.path(DER, "fo_kde.rds"))$sigma_main

days <- sort(unique(fire_px$date_local))
t0   <- min(days)
tt   <- as.numeric(fire_px$date_local - t0) + 1

# (1) One KDE per day. `edge = TRUE` applies Diggle's correction; `positive`
# guards against the tiny negative values that correction can introduce.
xy <- st_coordinates(fire_px)
dens_daily <- lapply(seq_along(days), function(k) {
  sel <- fire_px$date_local == days[k]
  pk  <- ppp(xy[sel, 1], xy[sel, 2], window = win, check = FALSE)
  density(pk, sigma = SIGMA_MAIN, edge = TRUE, positive = TRUE,
          dimyx = c(512, 256))
})
names(dens_daily) <- as.character(days)
cat("Daily KDE surfaces (per km2 per day), sigma =", SIGMA_MAIN, "m\n")
print(data.frame(day = as.character(days),
                 n = as.integer(table(factor(fire_px$date_local, levels = days))),
                 max_per_km2 = round(sapply(dens_daily, max) * 1e6, 2),
                 integral = round(sapply(dens_daily, integral), 1)),
      row.names = FALSE)

# (2) Separable spatio-temporal KDE.
ppp_marked <- ppp(xy[, 1], xy[, 2], window = win, marks = tt, check = FALSE)
LAMBDA_T <- 1.0
stk <- cache_rds("stkde_v2", spattemp.density(
  ppp_marked, h = SIGMA_MAIN, lambda = LAMBDA_T,
  tlim = c(0.5, length(days) + 0.5), sres = 256, tres = length(days),
  verbose = FALSE), refresh = TRUE)

cat("\nSTKDE tgrid (should be integer day indices):",
    paste(round(stk$tgrid, 4), collapse = ", "), "\n")
cat("Conditional slices integrate to:",
    paste(round(sapply(stk$z.cond, integral), 3), collapse = ", "), "\n")

saveRDS(list(stkde = stk, dens_daily = dens_daily, days = days, t0 = t0,
             tt = tt, lambda_t = LAMBDA_T, sigma = SIGMA_MAIN),
        file.path(DER, "fo_stkde.rds"))
cat("\n================ 05b_stkde.R complete ================\n")
