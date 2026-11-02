# =============================================================================
# R/huff.R -- the probabilistic spatial interaction model, in matrix form
# Take-home Exercise 3 | ISSS626
#
# For origin i and station j,
#
#   P_ij = A_j f(d_ij) / (A_0 + sum_k A_k f(d_ik)),   E_j = sum_i M_i P_ij
#
# with attractiveness A_j = exp(x_j' beta + r_j), distance decay f, and a local
# market M_i that is a weighted sum of one or more market columns. A_0 is an optional
# outside option, the choice of not using the rail network at all; with A_0 = 0 this
# is the Huff model. MCI::huff.shares() computes the same shares, as 05_calibrate.R
# confirms, but it works through a data frame and cannot be called thousands of
# times inside a calibration, so the calibration uses these matrix forms.
# =============================================================================

decay_g <- function(D, type = c("pow", "exp")) {
  type <- match.arg(type)
  if (type == "pow") log(D) else D
}

# lA: log attractiveness (length J); lf: log decay (I x J). Everything stays on the
# log scale, because A_j f(d_ij) spans many orders of magnitude across the island.
huff_prob <- function(lA, lf, lA0 = NULL) {
  u  <- sweep(lf, 2, lA, "+")
  mx <- apply(u, 1, max)
  if (!is.null(lA0)) mx <- pmax(mx, lA0)
  # A hexagon with no station inside the maximum access distance sends no one.
  mx[!is.finite(mx)] <- 0
  v  <- exp(u - mx)
  den <- rowSums(v)
  if (!is.null(lA0)) den <- den + exp(lA0 - mx)
  den[den == 0] <- 1
  v / den
}

# par = (lambda, beta..., log kappa..., [log A0])
unpack <- function(par, spec) {
  p <- ncol(spec$X); m <- ncol(spec$M)
  list(lambda = par[1],
       beta   = if (p > 0) par[1 + seq_len(p)] else numeric(0),
       kappa  = exp(par[1 + p + seq_len(m)]),
       lA0    = if (isTRUE(spec$use_A0)) par[2 + p + m] else NULL)
}

model_parts <- function(par, spec, D = spec$D, X = spec$X, M = spec$M, r = NULL) {
  q  <- unpack(par, spec)
  lA <- if (ncol(X) > 0) as.vector(X %*% q$beta) else rep(0, ncol(D))
  if (!is.null(r)) lA <- lA + r
  g  <- decay_g(D, spec$type)
  lf <- q$lambda * g
  # The choice set stops at a maximum access distance, because no one reaches their
  # boarding station from 20 km away, however small the probability a decay gives it.
  if (!is.null(spec$dmax)) lf[D > spec$dmax] <- -Inf
  P  <- huff_prob(lA, lf, q$lA0)
  mk <- as.vector(M %*% q$kappa)
  list(q = q, lA = lA, g = g, P = P, mk = mk, E = as.vector(crossprod(mk, P)))
}

expected <- function(par, spec, ...) model_parts(par, spec, ...)$E

# Squared error on the log scale. Ridership runs from under 200 to over 70,000 a day,
# and the stations this exercise cares about most are at the bottom of that range; on
# the raw scale the fit would be decided by a dozen interchanges.
loss <- function(par, spec) {
  E <- expected(par, spec)
  w <- if (is.null(spec$w)) 1 else spec$w
  sum(w * (log(spec$T) - log(pmax(E, 1e-6)))^2)
}

# Exact gradient of the loss. The validation and bootstrap run hundreds of
# calibrations, so finite differences, at one model run per parameter, are avoided.
grad_loss <- function(par, spec) {
  z <- model_parts(par, spec)
  P <- z$P; g <- z$g; mk <- z$mk; E <- pmax(z$E, 1e-6)
  w <- -2 * (if (is.null(spec$w)) 1 else spec$w) * (log(spec$T) - log(E)) / E
  gbar <- rowSums(P * g)
  out <- sum(w * (colSums(mk * P * g) - as.vector(crossprod(mk * gbar, P))))
  X <- spec$X
  if (ncol(X) > 0) {
    xbar <- P %*% X
    cs <- colSums(mk * P)
    out <- c(out, vapply(seq_len(ncol(X)), function(m) {
      sum(w * (X[, m] * cs - as.vector(crossprod(mk * xbar[, m], P))))
    }, numeric(1)))
  }
  out <- c(out, vapply(seq_len(ncol(spec$M)), function(m) {
    sum(w * z$q$kappa[m] * as.vector(crossprod(spec$M[, m], P)))
  }, numeric(1)))
  if (isTRUE(spec$use_A0)) {
    P0 <- 1 - rowSums(P)
    out <- c(out, sum(w * (-as.vector(crossprod(mk * P0, P)))))
  }
  out
}

start_par <- function(spec, lambda0) {
  c(lambda0, rep(0, ncol(spec$X)), rep(-1, ncol(spec$M)),
    if (isTRUE(spec$use_A0)) -2)
}

fit_huff <- function(spec, lambda0 = if (spec$type == "pow") -1.5 else -0.3,
                     par0 = NULL, maxit = 300, fix_lambda = NULL) {
  if (is.null(par0)) par0 <- start_par(spec, lambda0)
  if (!is.null(fix_lambda)) {
    # Holds lambda at a given value, for the profile and for the literature's -2.
    fn <- function(th) loss(c(fix_lambda, th), spec)
    gr <- function(th) grad_loss(c(fix_lambda, th), spec)[-1]
    fit <- optim(par0[-1], fn, gr, method = "BFGS", control = list(maxit = maxit, reltol = 1e-8))
    fit$par <- c(fix_lambda, fit$par)
  } else {
    fit <- optim(par0, loss, grad_loss, spec = spec, method = "BFGS",
                 control = list(maxit = maxit, reltol = 1e-8))
  }
  fit$E <- expected(fit$par, spec)
  fit$spec <- spec
  fit
}

fit_stats <- function(T, E) {
  e <- log(T) - log(E)
  tibble::tibble(rmsle = sqrt(mean(e^2)),
                 r2_log = 1 - sum(e^2) / sum((log(T) - mean(log(T)))^2),
                 mape = median(abs(E / T - 1)),
                 rho = cor(T, E, method = "spearman"))
}

# Station-specific attractiveness that reproduces observed ridership exactly, given
# lambda and the market. The multiplicative update is the one huff.fit() in MCI
# applies; here it runs on the whole matrix at once. A closed model allocates all of
# its market, so its total must equal total ridership before the update can converge;
# a log-scale fit leaves it a little short, so the market is rescaled first. Without
# an outside option r is defined only up to a constant, which is fixed by giving r a
# mean of zero, so that r = 0 is what a typical station has and a new station is
# given. `r0` allows a warm start.
closure_par <- function(par, spec) {
  if (isTRUE(spec$use_A0)) return(par)
  mk <- as.vector(spec$M %*% unpack(par, spec)$kappa)
  p <- ncol(spec$X); m <- ncol(spec$M)
  par[1 + p + seq_len(m)] <- par[1 + p + seq_len(m)] + log(sum(spec$T) / sum(mk))
  par
}

calibrate_r <- function(par, spec, r0 = NULL, tol = 1e-4, maxit = 3000) {
  par <- closure_par(par, spec)
  r <- if (is.null(r0)) rep(0, ncol(spec$D)) else r0
  for (it in seq_len(maxit)) {
    E <- expected(par, spec, r = r)
    ratio <- spec$T / pmax(E, 1e-9)
    gap <- max(abs(ratio - 1))
    if (gap < tol) break
    r <- r + log(ratio)
  }
  if (!isTRUE(spec$use_A0)) r <- r - mean(r)
  list(par = par, r = r, iterations = it, max_gap = gap)
}
