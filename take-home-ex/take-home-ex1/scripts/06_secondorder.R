# =============================================================================
# 06_secondorder.R -- second-order (interaction) analysis
# Take-home Exercise 1 | ISSS626
#
# The organising argument of this script is that the choice of NULL MODEL, not
# the choice of summary function, determines what a second-order analysis can
# legitimately conclude here. Three nested null models are used:
#
#   N1  Complete spatial randomness (homogeneous Poisson).
#       Tests "is the pattern random?" -- expected to be rejected trivially,
#       because the first-order analysis already established strong intensity
#       variation. Included to make that confounding explicit.
#
#   N2  Inhomogeneous Poisson with intensity FIXED at a kernel estimate.
#       Tests "once the large-scale intensity gradient is accounted for, do
#       points still attract one another?" This is the substantive question.
#
#   N3  Random labelling of dates over the observed locations.
#       Tests "is there space-time interaction?", holding both the spatial
#       configuration and the daily totals fixed.
#
# TWO IMPLEMENTATION DECISIONS THAT MATERIALLY AFFECT THE ANSWER
#
#   (a) VARIANCE STABILISATION. Global envelopes have constant width. Applied
#       to K, whose variance grows like r^4, a constant-width band is set by
#       the largest distances and is blind at small r. All envelopes are
#       therefore built on L = sqrt(K/pi), which is variance-stabilised.
#
#   (b) A SINGLE FIXED INTENSITY. Kinhom/Linhom will happily re-estimate
#       lambda separately for every simulated pattern (via a leave-one-out
#       kernel estimate). Doing so was tried and rejected: it makes the
#       observed and simulated statistics incommensurable, and on this pattern
#       it is numerically degenerate (the simulated L collapsed to zero). The
#       null model here is "an inhomogeneous Poisson process with intensity
#       lambda-hat", so lambda-hat is estimated ONCE from the data and the same
#       image is used to compute L for the observed pattern and for every
#       simulation.
#
#       This has a known and welcome bias direction: lambda-hat is fitted to
#       the observed clustering, which DEFLATES the observed Linhom, whereas
#       the simulated patterns are drawn from lambda-hat and so match it by
#       construction. The test is therefore CONSERVATIVE with respect to
#       detecting clustering -- a positive finding is not an artefact of the
#       procedure.
# =============================================================================

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(tidyr)
  library(spatstat.geom); library(spatstat.explore); library(spatstat.random)
})
source("R/utils.R")
sf_use_s2(FALSE)
set.seed(1234)

DER <- "data/derived"
P       <- readRDS(file.path(DER, "params.rds"))
win     <- readRDS(file.path(DER, "win.rds"))
ppp_px  <- readRDS(file.path(DER, "ppp_pxday.rds"))
ppp_det <- readRDS(file.path(DER, "ppp_det.rds"))
ppp_ev  <- readRDS(file.path(DER, "ppp_event.rds"))
fire_px <- readRDS(file.path(DER, "fire_pxday_sf.rds"))

NSIM   <- 199         # global envelope level = 2/(199+1) = 1%
CORR   <- "translate" # Ohser--Stoyan translation correction: unbiased, and ~20x
                      # faster than isotropic on this 4,178-vertex polygon,
                      # which is decisive at 199 simulations x 4 bandwidths.
RMAX_K <- 15000       # ~ 1/6 of the window's shorter dimension
RMAX_G <- 3000
r_K    <- seq(0, RMAX_K, length.out = 201)
r_G    <- seq(0, RMAX_G, length.out = 201)
r_g    <- seq(0, 8000,  length.out = 161)
SIG_LADDER <- c(1000, 2000, 5000, 10000)
SIGMA_MAIN <- 2000

cat(sprintf("n(pixel-day) = %d | window %.0f km2 | nsim = %d | correction = %s\n",
            npoints(ppp_px), area(win) / 1e6, NSIM, CORR))
cat(sprintf("Window bbox: %.0f x %.0f km; rmax(L) = %.0f km, rmax(G) = %.1f km\n\n",
            diff(win$xrange) / 1000, diff(win$yrange) / 1000,
            RMAX_K / 1000, RMAX_G / 1000))

# Helper: the KDE intensity image, restricted to the polygonal window so that
# simulated patterns inherit exactly the observed window.
lam_hat <- function(sig, pp = ppp_px) {
  density(pp, sigma = sig, edge = TRUE, positive = TRUE,
          dimyx = c(512, 256))[Window(pp), drop = FALSE]
}
# Helper: simulator for N2. The final `[w]` clips to the true polygon, because
# the intensity image's window is a pixel mask that is not a strict subset of it.
sim_ipp <- function(lam, w) function(X) rpoispp(lam)[w]

# =============================================================================
# N1a. NEAREST-NEIGHBOUR AND EMPTY-SPACE FUNCTIONS AGAINST CSR
# =============================================================================
# G = nearest-neighbour distance distribution (clustering -> G above Poisson).
# F = empty-space function (clustering -> F below Poisson).
# J = (1-G)/(1-F), equal to 1 under CSR.
# Kaplan--Meier edge correction, the recommended choice for these distance
# distributions in an irregular window.
env_G <- cache_rds("env_G_csr", envelope(ppp_px, Gest, r = r_G, nsim = NSIM,
                                         global = TRUE, correction = "km",
                                         savefuns = TRUE, verbose = FALSE))
env_F <- cache_rds("env_F_csr", envelope(ppp_px, Fest, r = r_G, nsim = NSIM,
                                         global = TRUE, correction = "km",
                                         savefuns = TRUE, verbose = FALSE))
env_J <- cache_rds("env_J_csr", envelope(ppp_px, Jest,
                                         r = seq(0, 1500, length.out = 151),
                                         nsim = NSIM, global = TRUE,
                                         correction = "km", savefuns = TRUE,
                                         verbose = FALSE))

nnd_obs <- mean(nndist(ppp_px))
nnd_csr <- cache_rds("nnd_csr_sim", replicate(NSIM,
             mean(nndist(runifpoint(npoints(ppp_px), win = win)))))
cat(sprintf("Mean nearest-neighbour distance: observed %.0f m; CSR %.0f m (2.5-97.5%%: %.0f-%.0f); ratio %.3f\n",
            nnd_obs, mean(nnd_csr), quantile(nnd_csr, .025),
            quantile(nnd_csr, .975), nnd_obs / mean(nnd_csr)))

# =============================================================================
# N1b. L AGAINST CSR -- the deliberately inadequate analysis
# =============================================================================
env_L_csr <- cache_rds("env_L_csr",
  envelope(ppp_px, Lest, r = r_K, nsim = NSIM, global = TRUE,
           correction = CORR, savefuns = TRUE, verbose = FALSE))
test_L_csr <- cache_rds("test_L_csr", dclf.test(env_L_csr))
mad_L_csr  <- cache_rds("mad_L_csr",  mad.test(env_L_csr))

cat("\n---- N1: L against CSR ----\n")
print(test_L_csr); print(mad_L_csr)

# =============================================================================
# N2. INHOMOGENEOUS L AGAINST AN INHOMOGENEOUS POISSON PROCESS
# =============================================================================
# The key modelling assumption is SCALE SEPARATION: intensity varies at scales
# >= sigma_lambda, interaction operates below it. Too small a sigma_lambda and
# the kernel estimate absorbs the clustering; too large and residual first-order
# variation is misread as interaction. This cannot be settled from a single
# realisation, so the analysis is run across a ladder of sigma_lambda and the
# conclusion is reported as a function of it.
#
# use.theory = FALSE: the reference curve is the mean of independent
# simulations, not the theoretical L = r. This matters because the estimator is
# slightly biased under the null (lambda-hat is not the true intensity), and the
# simulated mean captures that bias while the theoretical value does not.
# Each bandwidth is cached separately as well as in aggregate. A single
# bandwidth costs ~400 simulations, so per-bandwidth caching means an
# interrupted run resumes rather than restarting.
fit_linhom <- function(sig) cache_rds(sprintf("linhom_sigma%05d", sig), {
  lam <- lam_hat(sig)
  e <- envelope(ppp_px, Linhom, r = r_K, nsim = NSIM, nsim2 = NSIM,
                global = TRUE, use.theory = FALSE,
                simulate = sim_ipp(lam, win),
                funargs = list(lambda = lam, correction = CORR, normpower = 2),
                savefuns = TRUE, verbose = FALSE)
  list(sigma = sig, env = e, test = dclf.test(e), mad = mad.test(e))
})
linhom_ladder <- cache_rds("linhom_ladder", lapply(SIG_LADDER, fit_linhom))
names(linhom_ladder) <- paste0("sigma_", SIG_LADDER)

digest_env <- function(e) {
  d <- as.data.frame(e)
  ref <- if ("mmean" %in% names(d)) d$mmean else d$theo
  out <- which(d$obs > d$hi | d$obs < d$lo)
  list(d = d, ref = ref, out = out)
}
linhom_digest <- lapply(linhom_ladder, function(x) {
  g <- digest_env(x$env)
  tibble(sigma_lambda_km = x$sigma / 1000,
         dclf_p = x$test$p.value, mad_p = x$mad$p.value,
         r_first_above_km = if (length(g$out)) g$d$r[min(g$out)] / 1000 else NA_real_,
         r_last_above_km  = if (length(g$out)) g$d$r[max(g$out)] / 1000 else NA_real_,
         pct_r_outside    = 100 * length(g$out) / nrow(g$d))
}) |> bind_rows()

cat("\n---- N2: Linhom against inhomogeneous Poisson, by intensity bandwidth ----\n")
print(as.data.frame(linhom_digest |> mutate(across(where(is.numeric), ~round(.x, 4)))),
      row.names = FALSE)

# Sensitivity to the renormalisation power, which changed default in
# spatstat.explore 3.8-1 and is therefore worth showing explicitly.
np_check <- cache_rds("np_check", {
  lam <- lam_hat(SIGMA_MAIN)
  bind_rows(lapply(1:2, function(np) {
    L <- Linhom(ppp_px, lambda = lam, r = r_K, correction = CORR, normpower = np)
    d <- as.data.frame(L)
    tibble(normpower = np,
           L_at_500  = approx(d$r, d$trans, 500)$y,
           L_at_2000 = approx(d$r, d$trans, 2000)$y,
           L_at_10000 = approx(d$r, d$trans, 10000)$y)
  }))
})
cat("\n---- Sensitivity of Linhom to normpower (observed pattern) ----\n")
print(as.data.frame(np_check |> mutate(across(where(is.numeric), ~round(.x, 1)))),
      row.names = FALSE)

# =============================================================================
# INHOMOGENEOUS PAIR CORRELATION: naming a characteristic scale
# =============================================================================
# K and L are cumulative, so a large value at r = 5 km may be entirely due to
# structure at 500 m. g(r) isolates interaction AT distance r, which is what is
# needed to name a cluster size.
env_g <- cache_rds("env_pcfinhom", {
  lam <- lam_hat(SIGMA_MAIN)
  envelope(ppp_px, pcfinhom, r = r_g, nsim = NSIM, nsim2 = NSIM,
           global = TRUE, use.theory = FALSE,
           simulate = sim_ipp(lam, win),
           funargs = list(lambda = lam, correction = CORR, normpower = 2,
                          divisor = "d"),
           savefuns = TRUE, verbose = FALSE)
})
gd <- as.data.frame(env_g)
inside_band <- which(gd$r > 500 & gd$obs <= gd$hi)
cat(sprintf("\nInhomogeneous pair correlation at sigma_lambda = %d m:\n", SIGMA_MAIN))
cat(sprintf("  g(500 m)  = %.2f   (envelope %.2f - %.2f)\n",
            approx(gd$r, gd$obs, 500)$y, approx(gd$r, gd$lo, 500)$y,
            approx(gd$r, gd$hi, 500)$y))
cat(sprintf("  g(1 km)   = %.2f   g(2 km) = %.2f   g(5 km) = %.2f\n",
            approx(gd$r, gd$obs, 1000)$y, approx(gd$r, gd$obs, 2000)$y,
            approx(gd$r, gd$obs, 5000)$y))
cat(sprintf("  observed re-enters the envelope at r = %s m\n",
            if (length(inside_band)) n_fmt(gd$r[min(inside_band)]) else "never"))

# =============================================================================
# ROBUSTNESS TO THE DEFINITION OF A POINT EVENT
# =============================================================================
# If the conclusion depends on the choice of analysis unit it is not a
# conclusion about fire, so the test is repeated on all three units.
unit_check <- cache_rds("unit_check", lapply(
  list(detections = ppp_det, `pixel-day` = ppp_px, events = ppp_ev),
  function(pp) {
    lam <- lam_hat(SIGMA_MAIN, pp)
    e <- envelope(pp, Linhom, r = r_K, nsim = 99, nsim2 = 99, global = TRUE,
                  use.theory = FALSE, simulate = sim_ipp(lam, Window(pp)),
                  funargs = list(lambda = lam, correction = CORR, normpower = 2),
                  savefuns = TRUE, verbose = FALSE)
    list(env = e, test = dclf.test(e), n = npoints(pp))
  }))

unit_digest <- bind_rows(lapply(names(unit_check), function(nm) {
  x <- unit_check[[nm]]; g <- digest_env(x$env)
  tibble(unit = nm, n = x$n, dclf_p = x$test$p.value,
         pct_r_outside = 100 * length(g$out) / nrow(g$d))
}))
cat("\n---- Robustness: Linhom by definition of a point event ----\n")
print(as.data.frame(unit_digest |> mutate(across(where(is.numeric), ~round(.x, 4)))),
      row.names = FALSE)

# Robustness to the study window: shrinking the window by an internal buffer
# removes the boundary zone entirely, which is the strongest possible check on
# whether edge correction is doing the work.
window_check <- cache_rds("window_check", {
  bnd <- readRDS(file.path(DER, "bnd_proj.rds"))
  bind_rows(lapply(c(0, 5000, 10000), function(b) {
    w2 <- if (b == 0) win else as.owin(st_geometry(st_buffer(bnd, -b)))
    pp2 <- ppp_px[w2]
    lam2 <- lam_hat(SIGMA_MAIN, pp2)
    e <- envelope(pp2, Linhom, r = r_K, nsim = 99, nsim2 = 99, global = TRUE,
                  use.theory = FALSE, simulate = sim_ipp(lam2, w2),
                  funargs = list(lambda = lam2, correction = CORR, normpower = 2),
                  savefuns = TRUE, verbose = FALSE)
    g <- digest_env(e)
    tibble(inner_buffer_km = b / 1000, n = npoints(pp2),
           area_km2 = area(w2) / 1e6,
           dclf_p = dclf.test(e)$p.value,
           pct_r_outside = 100 * length(g$out) / nrow(g$d))
  }))
})
cat("\n---- Robustness: shrinking the study window ----\n")
print(as.data.frame(window_check |> mutate(across(where(is.numeric), ~round(.x, 3)))),
      row.names = FALSE)

# =============================================================================
# N3. SPACE-TIME INTERACTION: MONTE-CARLO KNOX TEST
# =============================================================================
# Research question: are detections close in space also close in time, beyond
# what the separate spatial and temporal patterns imply?
#
# The null keeps the observed locations and the observed multiset of dates
# exactly as they are, and randomly re-assigns dates to locations. An excess of
# space-time-close pairs therefore cannot be explained by the spatial pattern
# alone, nor by the daily totals alone -- only by dependence between where and
# when. This random-labelling null needs no model for either marginal pattern,
# which is why it is preferred here to a space-time Poisson null.
xy <- cbind(ppp_px$x, ppp_px$y)
tt <- as.numeric(fire_px$date_local - min(fire_px$date_local))

knox_grid <- expand.grid(space_m = c(1000, 2000, 5000),
                         time_d  = c(0, 1, 2)) |> as_tibble()

knox_res <- cache_rds("knox_res", {
  # Spatial neighbour lists depend only on space_m, so they are built once per
  # radius and reused across all 999 date permutations.
  nbs <- lapply(unique(knox_grid$space_m), function(s) {
    nb <- dbscan::frNN(xy, eps = s)
    i  <- rep(seq_along(nb$id), lengths(nb$id)); j <- unlist(nb$id)
    keep <- i < j                      # count each unordered pair once
    list(i = i[keep], j = j[keep])
  })
  names(nbs) <- as.character(unique(knox_grid$space_m))
  bind_rows(lapply(seq_len(nrow(knox_grid)), function(k) {
    s <- knox_grid$space_m[k]; td <- knox_grid$time_d[k]
    ij <- nbs[[as.character(s)]]; i <- ij$i; j <- ij$j
    obs <- sum(abs(tt[i] - tt[j]) <= td)
    sim <- replicate(999, { tp <- sample(tt); sum(abs(tp[i] - tp[j]) <= td) })
    tibble(space_m = s, time_d = td, n_pairs_space = length(i),
           observed = obs, sim_mean = mean(sim), sim_sd = sd(sim),
           ratio = obs / mean(sim),
           p_value = (1 + sum(sim >= obs)) / (1 + length(sim)))
  }))
})
cat("\n---- N3: Monte-Carlo Knox test (999 permutations of the date labels) ----\n")
print(as.data.frame(knox_res |> mutate(across(where(is.numeric), ~round(.x, 3)))),
      row.names = FALSE)

# Per-lag view: does the excess of nearby pairs decay as the time lag grows?
st_by_lag <- cache_rds("st_by_lag", {
  nb <- dbscan::frNN(xy, eps = 5000)
  i  <- rep(seq_along(nb$id), lengths(nb$id)); j <- unlist(nb$id)
  keep <- i < j; i <- i[keep]; j <- j[keep]
  dt <- abs(tt[i] - tt[j])
  bind_rows(lapply(0:5, function(u) {
    obs <- sum(dt == u)
    sim <- replicate(499, { tp <- sample(tt); sum(abs(tp[i] - tp[j]) == u) })
    tibble(lag_days = u, observed = obs, expected = mean(sim),
           sd = sd(sim), ratio = obs / mean(sim),
           p_value = (1 + sum(sim >= obs)) / (1 + length(sim)))
  }))
})
cat("\n---- Space-time pair excess by temporal lag (pairs within 5 km) ----\n")
print(as.data.frame(st_by_lag |> mutate(across(where(is.numeric), ~round(.x, 3)))),
      row.names = FALSE)

# Same test at a tighter spatial radius, to see whether interaction is a
# short-range or a landscape-scale phenomenon.
st_by_lag_1km <- cache_rds("st_by_lag_1km", {
  nb <- dbscan::frNN(xy, eps = 1000)
  i  <- rep(seq_along(nb$id), lengths(nb$id)); j <- unlist(nb$id)
  keep <- i < j; i <- i[keep]; j <- j[keep]
  dt <- abs(tt[i] - tt[j])
  bind_rows(lapply(0:5, function(u) {
    obs <- sum(dt == u)
    sim <- replicate(499, { tp <- sample(tt); sum(abs(tp[i] - tp[j]) == u) })
    tibble(lag_days = u, observed = obs, expected = mean(sim),
           ratio = obs / mean(sim),
           p_value = (1 + sum(sim >= obs)) / (1 + length(sim)))
  }))
})
cat("\n---- Same, pairs within 1 km ----\n")
print(as.data.frame(st_by_lag_1km |> mutate(across(where(is.numeric), ~round(.x, 3)))),
      row.names = FALSE)

saveRDS(list(env_G = env_G, env_F = env_F, env_J = env_J,
             env_L_csr = env_L_csr, test_L_csr = test_L_csr, mad_L_csr = mad_L_csr,
             linhom_ladder = linhom_ladder, linhom_digest = linhom_digest,
             np_check = np_check,
             env_g = env_g, gd = gd,
             unit_check = unit_check, unit_digest = unit_digest,
             window_check = window_check,
             knox_res = knox_res, st_by_lag = st_by_lag,
             st_by_lag_1km = st_by_lag_1km,
             nnd_obs = nnd_obs, nnd_csr = nnd_csr,
             nsim = NSIM, corr = CORR, r_K = r_K, r_G = r_G, r_g = r_g,
             sig_ladder = SIG_LADDER, sigma_main = SIGMA_MAIN),
        file.path(DER, "so_results.rds"))

cat("\n================ 06_secondorder.R complete ================\n")
