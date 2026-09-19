# =============================================================================
# 07_valueadded.R -- value-added analysis: parametric point-process modelling
# Take-home Exercise 1 | ISSS626
#
# WHY GO BEYOND THE CORE METHODS
#   The core first- and second-order analyses can establish THAT intensity
#   varies and THAT points cluster, but neither can separate the two, because a
#   kernel intensity estimate and an interaction effect are confounded by
#   construction: any clustering can be re-described as intensity variation at a
#   smaller bandwidth, and vice versa. Fitting explicit point-process models
#   breaks that circularity by giving the two components different functional
#   roles and letting the data arbitrate:
#
#     M0   homogeneous Poisson                      (no trend, no interaction)
#     M1c  inhomogeneous Poisson, canals only       (one covariate)
#     M1r  inhomogeneous Poisson, rivers only       (one covariate)
#     M1   inhomogeneous Poisson, all four, log-linear
#     M2   inhomogeneous Poisson, all four, quadratic
#     M3   Thomas cluster process with M1's trend   (trend AND interaction)
#
#   Three things are gained that the core analysis cannot provide:
#
#   (a) A MULTIVARIATE answer on covariates. rhohat gives marginal curves, and
#       the four distance covariates are mutually correlated, so a marginal
#       canal effect might simply be an accessibility effect wearing a disguise.
#       ppm holds the others constant.
#
#   (b) A TEST OF WHETHER GEOGRAPHY EXPLAINS THE CLUSTERING. Simulating from the
#       fitted M1/M2 gives null model N4 -- an inhomogeneous Poisson process
#       whose intensity is driven by *measured* drainage and access geography
#       rather than by a kernel smooth of the data itself. This avoids the
#       circularity of N2 entirely, because the trend cannot absorb small-scale
#       clustering: it has only four smooth covariates to work with.
#
#   (c) A GENERATIVE MODEL WITH INTERPRETABLE PARAMETERS. M3's parameters are
#       quantities a planner can use -- how many distinct clusters, how large
#       each is in metres, how many detections each contains, and how likely two
#       nearby detections are to belong to the same cluster. No summary function
#       provides these.
# =============================================================================

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(tidyr); library(tibble)
  library(spatstat.geom); library(spatstat.explore); library(spatstat.model)
  library(spatstat.random)
})
source("R/utils.R")
sf_use_s2(FALSE)
set.seed(1234)

DER <- "data/derived"
win    <- readRDS(file.path(DER, "win.rds"))
ppp_px <- readRDS(file.path(DER, "ppp_pxday.rds"))
covars <- readRDS(file.path(DER, "covars.rds"))

# Covariates are rescaled from metres to kilometres. This changes nothing
# statistically but makes every coefficient a "per kilometre" effect, which is
# the unit the interpretation needs.
Z <- lapply(covars, function(im) eval.im(im / 1000))
names(Z) <- names(covars)

r_K      <- seq(0, 15000, length.out = 201)
NSIM_GOF <- 199
CORR     <- "translate"
# ppm fits by maximum pseudolikelihood on a quadrature scheme. The scheme is
# refined beyond the default because the covariates vary over ~250 m while the
# window is ~215 km long; too coarse a scheme would alias the canal signal.
QUAD_ND  <- 256

# =============================================================================
# 1. MODEL FITTING AND COMPARISON
# =============================================================================
fits <- cache_rds("va_fits", list(
  M0  = ppm(ppp_px ~ 1, nd = QUAD_ND),
  M1c = ppm(ppp_px ~ d_canal, covariates = Z, nd = QUAD_ND),
  M1r = ppm(ppp_px ~ d_river, covariates = Z, nd = QUAD_ND),
  M1  = ppm(ppp_px ~ d_canal + d_river + d_road + d_settle,
            covariates = Z, nd = QUAD_ND),
  # M2 lets each effect bend. A log-linear distance effect forces intensity to
  # decay at a constant proportional rate for ever, which is implausible; a
  # quadratic allows the effect to flatten out or reverse.
  M2  = ppm(ppp_px ~ polynom(d_canal, 2) + polynom(d_river, 2) +
                     polynom(d_road, 2) + polynom(d_settle, 2),
            covariates = Z, nd = QUAD_ND)
))
M0 <- fits$M0; M1c <- fits$M1c; M1r <- fits$M1r; M1 <- fits$M1; M2 <- fits$M2

aic_tab <- tibble(
  model = c("M0  homogeneous Poisson",
            "M1c trend: distance to canal only",
            "M1r trend: distance to river only",
            "M1  trend: all four, log-linear",
            "M2  trend: all four, quadratic"),
  df  = vapply(fits, function(m) length(coef(m)), integer(1)),
  AIC = vapply(fits, AIC, numeric(1))
) |> mutate(dAIC_vs_M0 = AIC[1] - AIC, dAIC_vs_best = AIC - min(AIC))

cat("---- Poisson model comparison ----\n")
print(as.data.frame(aic_tab |> mutate(across(where(is.numeric), ~round(.x, 1)))),
      row.names = FALSE)

co <- as.data.frame(summary(M1)$coefs.SE.CI)
cat("\n---- M1 coefficients (log intensity per km) ----\n")
print(co |> mutate(across(where(is.numeric), ~round(.x, 4))))

# The same numbers read multiplicatively: the factor by which intensity changes
# per additional kilometre, and the distance over which it halves.
eff <- co |>
  rownames_to_column("term") |>
  filter(term != "(Intercept)") |>
  transmute(term,
            estimate = Estimate,
            per_km_multiplier = exp(Estimate),
            ci_lo = exp(CI95.lo), ci_hi = exp(CI95.hi),
            halving_distance_km = log(0.5) / Estimate,
            z = Zval)
cat("\n---- M1 effects as multipliers per additional km ----\n")
print(as.data.frame(eff |> mutate(across(where(is.numeric), ~round(.x, 3)))),
      row.names = FALSE)

# =============================================================================
# 2. N4: DOES MEASURED GEOGRAPHY EXPLAIN THE CLUSTERING?
# =============================================================================
# envelope() applied to a fitted ppm simulates from THAT model, so the null is
# "an inhomogeneous Poisson process with this covariate-driven trend". If the
# observed L stays outside the envelope, four smooth covariates do not account
# for the clustering and a model with genuine interaction is required.
gof_M1 <- cache_rds("gof_M1", envelope(M1, Lest, r = r_K, nsim = NSIM_GOF,
                                       global = TRUE, correction = CORR,
                                       savefuns = TRUE, verbose = FALSE))
gof_M1_test <- cache_rds("gof_M1_test", dclf.test(gof_M1))
cat("\n---- N4: goodness of fit of M1 (null = fitted inhomogeneous Poisson) ----\n")
print(gof_M1_test)

gof_M2 <- cache_rds("gof_M2", envelope(M2, Lest, r = r_K, nsim = NSIM_GOF,
                                       global = TRUE, correction = CORR,
                                       savefuns = TRUE, verbose = FALSE))
gof_M2_test <- cache_rds("gof_M2_test", dclf.test(gof_M2))
cat("\n---- N4b: goodness of fit of M2 (flexible trend) ----\n")
print(gof_M2_test)

# =============================================================================
# 3. M3: THOMAS CLUSTER PROCESS WITH THE SAME TREND
# =============================================================================
# A Thomas process supposes unobserved "parent" locations at intensity kappa,
# each producing offspring scattered isotropically with standard deviation
# `scale`. Read substantively: parents are ignition foci, offspring are the
# detections they generate, and `scale` is the effective size of a burning
# patch. Fitted by minimum contrast on the inhomogeneous K, so the trend comes
# from M1 and the residual second-order structure is attributed to clustering.
M3 <- cache_rds("M3_thomas",
  kppm(ppp_px ~ d_canal + d_river + d_road + d_settle,
       clusters = "Thomas", covariates = Z, nd = QUAD_ND,
       statistic = "K", statargs = list(correction = CORR)))

cat("\n---- M3 Thomas cluster process ----\n")
print(M3)

cp <- M3$clustpar
clust_tab <- tibble(
  quantity = c("kappa: parent intensity",
               "implied number of clusters in the window",
               "scale: cluster standard deviation",
               "~68% of offspring within",
               "~95% of offspring within",
               "mu: mean offspring per parent",
               "phi: cluster strength",
               "psib: sibling probability"),
  value = c(sprintf("%.3g per m2  (%.2f per 1,000 km2)",
                    cp[["kappa"]], cp[["kappa"]] * 1e9),
            sprintf("%.0f", cp[["kappa"]] * area(win)),
            sprintf("%.0f m", cp[["scale"]]),
            sprintf("%.1f km of the parent", cp[["scale"]] / 1000),
            sprintf("%.1f km of the parent", 2 * cp[["scale"]] / 1000),
            sprintf("%.0f detections", mean(M3$mu)),
            sprintf("%.2f", tryCatch(unname(clusterstrength(M3)),
                                     error = function(e) NA_real_)),
            sprintf("%.3f", tryCatch(unname(psib(M3)),
                                     error = function(e) NA_real_)))
)
cat("\n---- M3 parameters, read substantively ----\n")
print(as.data.frame(clust_tab), row.names = FALSE)

# N5: is even trend-plus-clustering adequate?
gof_M3 <- cache_rds("gof_M3", envelope(M3, Lest, r = r_K, nsim = NSIM_GOF,
                                       global = TRUE, correction = CORR,
                                       savefuns = TRUE, verbose = FALSE))
gof_M3_test <- cache_rds("gof_M3_test", dclf.test(gof_M3))
cat("\n---- N5: goodness of fit of M3 (null = fitted Thomas process) ----\n")
print(gof_M3_test)

# =============================================================================
# 4. RESIDUAL DIAGNOSTICS
# =============================================================================
# Smoothed Pearson residuals show WHERE the trend model is wrong, which is more
# informative than a single goodness-of-fit p-value: a covariate model can be
# rejected overall and still be adequate over most of the window.
res_sm <- cache_rds("res_sm", Smooth(residuals(M1, type = "pearson"),
                                     sigma = 5000))
cat(sprintf("\nSmoothed Pearson residuals of M1: min %.4f, max %.4f\n",
            min(res_sm), max(res_sm)))

# Partial residuals: the non-parametric shape of one covariate effect after the
# others are accounted for, against the fitted log-linear form.
prd <- cache_rds("parres_all", lapply(names(Z), function(nm)
  tryCatch(parres(M1, nm), error = function(e) NULL)))
names(prd) <- names(Z)

# =============================================================================
# 5. WHAT THE MODELS IMPLY SPATIALLY
# =============================================================================
trend_M1 <- predict(M1, type = "trend", ngrid = c(512, 256))
trend_M2 <- predict(M2, type = "trend", ngrid = c(512, 256))
cat(sprintf("\nM1 fitted trend: %.4f to %.4f per km2 (ratio %.0f:1)\n",
            min(trend_M1) * 1e6, max(trend_M1) * 1e6,
            max(trend_M1) / min(trend_M1)))
cat(sprintf("M2 fitted trend: %.4f to %.4f per km2 (ratio %.0f:1)\n",
            min(trend_M2) * 1e6, max(trend_M2) * 1e6,
            max(trend_M2) / min(trend_M2)))

saveRDS(list(fits = fits, M0 = M0, M1 = M1, M2 = M2, M3 = M3,
             aic_tab = aic_tab, coefs = co, eff = eff,
             gof_M1 = gof_M1, gof_M1_test = gof_M1_test,
             gof_M2 = gof_M2, gof_M2_test = gof_M2_test,
             gof_M3 = gof_M3, gof_M3_test = gof_M3_test,
             res_sm = res_sm, parres = prd,
             trend_M1 = trend_M1, trend_M2 = trend_M2,
             clustpar = cp, mu = M3$mu, clust_tab = clust_tab,
             nsim = NSIM_GOF, quad_nd = QUAD_ND),
        file.path(DER, "va_results.rds"))

cat("\n================ 07_valueadded.R complete ================\n")
