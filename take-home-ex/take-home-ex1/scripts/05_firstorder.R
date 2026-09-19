# =============================================================================
# 05_firstorder.R -- first-order (intensity) analysis
# Take-home Exercise 1 | ISSS626
#
# Computes and caches every first-order result used by the technical report,
# and prints a console digest so the numbers can be checked independently of
# the rendered document.
# =============================================================================

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(tidyr); library(lubridate)
  library(spatstat.geom); library(spatstat.explore); library(spatstat.random)
  library(sparr)
})
source("R/utils.R")
sf_use_s2(FALSE)
set.seed(1234)

DER <- "data/derived"
P         <- readRDS(file.path(DER, "params.rds"))
fire_sf   <- readRDS(file.path(DER, "fire_sf.rds"))
fire_px   <- readRDS(file.path(DER, "fire_pxday_sf.rds"))
fire_ev   <- readRDS(file.path(DER, "fire_event.rds"))
fire_all  <- readRDS(file.path(DER, "fire_all_sf.rds"))
win       <- readRDS(file.path(DER, "win.rds"))
bnd       <- readRDS(file.path(DER, "bnd_proj.rds"))
covars    <- readRDS(file.path(DER, "covars.rds"))
ppp_det   <- readRDS(file.path(DER, "ppp_det.rds"))
ppp_px    <- readRDS(file.path(DER, "ppp_pxday.rds"))
ppp_ev    <- readRDS(file.path(DER, "ppp_event.rds"))

AREA_KM2 <- area(win) / 1e6
cat(sprintf("\nWindow area %.1f km2 | detections %d | pixel-days %d | events %d\n",
            AREA_KM2, npoints(ppp_det), npoints(ppp_px), npoints(ppp_ev)))
cat(sprintf("Mean intensity (pixel-day unit) = %.4f per km2 = 1 per %.1f km2\n\n",
            npoints(ppp_px) / AREA_KM2, AREA_KM2 / npoints(ppp_px)))

# =============================================================================
# 1. TEMPORAL FIRST-ORDER STRUCTURE
# =============================================================================
# Question: does the rate of detection vary over the observation period, and
# how much of any variation is attributable to the observation process (which
# overpass, day or night, cloud) rather than to fire activity?

daily <- fire_px |> st_drop_geometry() |>
  count(date_local, name = "n_pxday") |>
  full_join(fire_sf |> st_drop_geometry() |> count(date_local, name = "n_det"),
            by = "date_local") |>
  full_join(fire_sf |> st_drop_geometry() |>
              count(date_local, overpass) |>
              pivot_wider(names_from = overpass, values_from = n,
                          names_prefix = "n_", values_fill = 0),
            by = "date_local") |>
  arrange(date_local) |>
  mutate(across(starts_with("n_"), ~replace_na(.x, 0L)))

# Event onsets: an event is counted on the day it is first detected, which is
# a closer proxy for *new ignitions* than the detection count.
onsets <- fire_ev |> count(start, name = "n_new_events") |> rename(date_local = start)
daily  <- daily |> left_join(onsets, by = "date_local") |>
  mutate(n_new_events = replace_na(n_new_events, 0L))

# Multi-sensor comparison: the four products have different overpass times and
# footprints, so agreement between them is evidence that a day-to-day change
# is real rather than an artefact of one satellite's sampling.
daily_sensor <- fire_all |> st_drop_geometry() |>
  count(date_local, firms_source) |>
  pivot_wider(names_from = firms_source, values_from = n, values_fill = 0)

hourly <- fire_sf |> st_drop_geometry() |>
  mutate(hr = floor(hour_local)) |> count(hr, overpass)

cat("---- Daily counts ----\n"); print(as.data.frame(daily), row.names = FALSE)
cat("\n---- Daily counts by product ----\n"); print(as.data.frame(daily_sensor), row.names = FALSE)
cat("\n---- Local hour of acquisition ----\n"); print(as.data.frame(hourly), row.names = FALSE)

cv <- sd(daily$n_pxday) / mean(daily$n_pxday)
cat(sprintf("\nDay-to-day variation in pixel-day counts: mean %.0f, sd %.0f, CV %.2f, range %d-%d\n",
            mean(daily$n_pxday), sd(daily$n_pxday), cv,
            min(daily$n_pxday), max(daily$n_pxday)))
# Poisson dispersion check: under a constant-rate Poisson process the daily
# counts would have variance equal to their mean.
disp <- var(daily$n_pxday) / mean(daily$n_pxday)
pois_test <- sum((daily$n_pxday - mean(daily$n_pxday))^2) / mean(daily$n_pxday)
cat(sprintf("Index of dispersion = %.1f (=1 under constant-rate Poisson); chi2_%d = %.1f, p = %s\n",
            disp, nrow(daily) - 1, pois_test,
            pv(pchisq(pois_test, nrow(daily) - 1, lower.tail = FALSE))))

saveRDS(list(daily = daily, daily_sensor = daily_sensor, hourly = hourly,
             cv = cv, disp = disp,
             pois_chi2 = pois_test,
             pois_p = pchisq(pois_test, nrow(daily) - 1, lower.tail = FALSE)),
        file.path(DER, "fo_temporal.rds"))

# =============================================================================
# 2. QUADRAT ANALYSIS
# =============================================================================
# A deliberately coarse first look. The chi-squared quadrat test asks only
# whether counts are consistent with a *homogeneous* Poisson process; it says
# nothing about interaction, and its answer depends on quadrat size. Both
# limitations are demonstrated rather than asserted, by running a ladder of
# quadrat sizes.
quad_ladder <- lapply(list(c(3, 8), c(5, 12), c(8, 20), c(12, 30)), function(nxy) {
  qt <- quadrat.test(ppp_px, nx = nxy[1], ny = nxy[2], method = "MonteCarlo",
                     nsim = 999)
  qc <- quadratcount(ppp_px, nx = nxy[1], ny = nxy[2])
  v  <- as.vector(as.table(qc))
  tibble(nx = nxy[1], ny = nxy[2],
         n_quadrats_used = length(v),
         mean_count = mean(v), var_count = var(v),
         vmr = var(v) / mean(v),
         statistic = unname(qt$statistic), p_value = qt$p.value)
}) |> bind_rows()

cat("\n---- Quadrat tests (Monte Carlo, 999 sim, null = CSR) ----\n")
print(as.data.frame(quad_ladder), row.names = FALSE)
saveRDS(quad_ladder, file.path(DER, "fo_quadrat.rds"))

# =============================================================================
# 3. KERNEL DENSITY ESTIMATION: BANDWIDTH IS THE DECISION THAT MATTERS
# =============================================================================
bw_tab <- tibble(
  selector = c("bw.diggle", "bw.ppl", "bw.CvL", "bw.scott (x)", "bw.scott (y)"),
  sigma_m  = c(as.numeric(bw.diggle(ppp_px)), as.numeric(bw.ppl(ppp_px)),
               as.numeric(bw.CvL(ppp_px)),    as.numeric(bw.scott(ppp_px))[1],
               as.numeric(bw.scott(ppp_px))[2]),
  criterion = c("minimises MSE for a Cox process",
                "maximises likelihood cross-validation",
                "Cronie--van Lieshout: matches expected to observed area",
                "normal reference rule, x direction",
                "normal reference rule, y direction")
)
cat("\n---- Automatic bandwidth selectors ----\n")
print(as.data.frame(bw_tab |> mutate(sigma_m = round(sigma_m))), row.names = FALSE)
cat(sprintf("\nSensor resolution floor: VIIRS I-band pixel = %d m.\n", P$pixel_grid_m))
cat("Selectors below that floor describe the detection grid, not the fire field.\n")

SIGMA_LADDER <- c(500, 1000, 2000, 5000, 10000)
SIGMA_MAIN   <- 2000

dens <- lapply(SIGMA_LADDER, function(s)
  density(ppp_px, sigma = s, edge = TRUE, positive = TRUE, dimyx = c(512, 256)))
names(dens) <- paste0("sigma_", SIGMA_LADDER)

# Diagnostic: a density surface should integrate back to the point count.
# The shortfall quantifies the combined discretisation + edge-correction error.
int_check <- tibble(sigma = SIGMA_LADDER,
                    integral = vapply(dens, integral, numeric(1)),
                    n = npoints(ppp_px)) |>
  mutate(pct_error = 100 * (integral - n) / n,
         max_per_km2 = vapply(dens, function(d) max(d) * 1e6, numeric(1)),
         ratio_max_to_mean = max_per_km2 / (npoints(ppp_px) / AREA_KM2))
cat("\n---- Density surface diagnostics ----\n")
print(as.data.frame(int_check |> mutate(across(where(is.numeric), ~round(.x, 3)))),
      row.names = FALSE)

# Adaptive (variable) bandwidth: the intensity gradient across this window is
# very steep, and a single global bandwidth must either oversmooth the peaks
# or undersmooth the sparse north. Abramson's rule lets the bandwidth shrink
# where data are dense.
dens_adapt <- cache_rds("dens_adapt", {
  bw <- bw.abram(ppp_px, h0 = SIGMA_MAIN)
  densityAdaptiveKernel(ppp_px, bw = bw, edge = TRUE, dimyx = c(512, 256))
})
cat(sprintf("\nAdaptive KDE (Abramson, h0 = %d m): max %.3f per km2, integral %.0f\n",
            SIGMA_MAIN, max(dens_adapt) * 1e6, integral(dens_adapt)))

saveRDS(list(bw_tab = bw_tab, sigma_ladder = SIGMA_LADDER, sigma_main = SIGMA_MAIN,
             dens = dens, dens_adapt = dens_adapt, int_check = int_check),
        file.path(DER, "fo_kde.rds"))

# =============================================================================
# 4. ADMINISTRATIVE AGGREGATION (the unit planners actually act on)
# =============================================================================
kec <- readRDS(file.path(DER, "pp_adm3_cod.rds")) |> st_transform(P$crs_proj) |>
  st_make_valid()
kec$area_km2 <- as.numeric(st_area(kec)) / 1e6
kec_counts <- fire_px |>
  st_join(kec |> select(adm3_name), join = st_within) |>
  st_drop_geometry() |> count(adm3_name, name = "n_pxday")
kec <- kec |> left_join(kec_counts, by = "adm3_name") |>
  mutate(n_pxday = tidyr::replace_na(n_pxday, 0L),
         rate_per_1000km2 = n_pxday / area_km2 * 1000,
         share_pct = 100 * n_pxday / sum(n_pxday))
cat("\n---- Detections by kecamatan (district) ----\n")
print(as.data.frame(kec |> st_drop_geometry() |>
        select(adm3_name, area_km2, n_pxday, rate_per_1000km2, share_pct) |>
        mutate(across(where(is.numeric), ~round(.x, 1))) |>
        arrange(desc(rate_per_1000km2))), row.names = FALSE)
saveRDS(kec, file.path(DER, "fo_kecamatan.rds"))

# =============================================================================
# 5. INTENSITY AS A FUNCTION OF COVARIATES (non-parametric)
# =============================================================================
# rhohat() estimates lambda(Z) -- intensity as a smooth function of a covariate
# -- without assuming a functional form. It answers "where in covariate space
# is detection intensity high?", which is the first step towards explaining the
# first-order inhomogeneity rather than merely mapping it.
rho <- cache_rds("rhohat_all", lapply(names(covars), function(nm)
  rhohat(ppp_px, covars[[nm]], confidence = 0.95)))
names(rho) <- names(covars)

cat("\n---- rhohat summary: intensity per km2 at selected covariate values ----\n")
rho_digest <- lapply(names(rho), function(nm) {
  r <- rho[[nm]]; d <- as.data.frame(r)
  zz <- quantile(d[[1]], c(0, .1, .25, .5, .75, .9, 1), na.rm = TRUE)
  tibble(covariate = nm,
         q = names(zz), z_m = round(as.numeric(zz)),
         rho_per_km2 = round(approx(d[[1]], d$rho, as.numeric(zz))$y * 1e6, 3))
}) |> bind_rows()
print(as.data.frame(rho_digest), row.names = FALSE)
saveRDS(rho, file.path(DER, "fo_rhohat.rds"))

# Covariate values at the detections vs across the window: a simple, readable
# summary of the same signal that does not depend on smoothing choices.
cov_at_pts <- sapply(covars, function(im) lookup.im(im, ppp_px$x, ppp_px$y, naok = TRUE))
cov_in_win <- sapply(covars, function(im) as.vector(im[win, drop = TRUE]))
cov_compare <- tibble(
  covariate = names(covars),
  median_at_detections = round(apply(cov_at_pts, 2, median, na.rm = TRUE)),
  median_in_window     = round(apply(cov_in_win, 2, median, na.rm = TRUE))
) |> mutate(ratio = round(median_at_detections / median_in_window, 2))
cat("\n---- Covariate values at detections vs across the window (m) ----\n")
print(as.data.frame(cov_compare), row.names = FALSE)
# Formal tests that intensity depends on each covariate.
#   cdf.test(): Kolmogorov--Smirnov comparison of the covariate distribution at
#     the data points against its distribution over the window, under the null
#     of a homogeneous Poisson process.
#   berman.test(Z1): the standardised mean covariate value at the points, same
#     null; its SIGN tells the direction of the effect.
# CAVEAT, stated here and again in the report: both tests assume independent
# points. The pattern is strongly clustered, so the effective sample size is far
# smaller than n and these p-values are anti-conservative. They are reported as
# evidence of direction and rough strength, not as calibrated probabilities.
cov_tests <- lapply(names(covars), function(nm) {
  ks <- cdf.test(ppp_px, covars[[nm]], test = "ks")
  bm <- berman.test(ppp_px, covars[[nm]], which = "Z1")
  tibble(covariate = nm,
         ks_D = unname(ks$statistic), ks_p = ks$p.value,
         Z1 = unname(bm$statistic),  Z1_p = bm$p.value,
         direction = if_else(unname(bm$statistic) < 0,
                             "closer than expected", "further than expected"))
}) |> bind_rows()
cat("\n---- Formal covariate-dependence tests (null: homogeneous Poisson) ----\n")
print(as.data.frame(cov_tests |> mutate(across(where(is.numeric), ~signif(.x, 4)))),
      row.names = FALSE)

saveRDS(list(cov_at_pts = cov_at_pts, cov_compare = cov_compare,
             cov_tests = cov_tests),
        file.path(DER, "fo_covsummary.rds"))

# =============================================================================
# 6. SPATIO-TEMPORAL KDE -> see scripts/05b_stkde.R
# =============================================================================
# The spatio-temporal estimates live in their own script, because two
# complementary estimators are compared there (independent daily kernels, and a
# separable space-time kernel) and because the temporal grid needs care: with
# tlim = c(1, 8) and tres = 8, sparr places its grid at cell centres 1.4375,
# 2.3125, ... so no panel corresponds to an actual day. 05b_stkde.R uses
# tlim = c(0.5, 8.5) instead, putting the grid on integer day indices.
#
# 05b_stkde.R writes data/derived/fo_stkde.rds and must be run AFTER this one.

cat("\n================ 05_firstorder.R complete ================\n")
