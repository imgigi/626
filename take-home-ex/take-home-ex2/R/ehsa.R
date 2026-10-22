# =============================================================================
# R/ehsa.R -- emerging hot spot analysis with a period-conditional null
# Take-home Exercise 2 | ISSS626
#
# The space-time Gi* below is the statistic sfdep::emerging_hotspot_analysis()
# computes, written out so that it can be checked against sfdep value for value
# (scripts/04_ehsa.R does the check) and so that the one part that is changed
# is visible. That part is the permutation null.
#
# sfdep draws the pseudo-neighbours of a bin from every bin in the cube, at any
# location and in any period, and replaces the bin's own value as well. When
# the national level of the variable trends, as conflict in Myanmar does, the
# draws for an early, quiet quarter come mostly from later, busier quarters, so
# almost every township looks unusually low; in the busiest quarters the reverse
# happens. The significance of a bin then tracks the national trend rather than
# the township's standing among its contemporaries.
#
# Here the township's own values are held fixed and its neighbours are replaced
# by a random set of other townships, read in the same periods as the real
# neighbours. Each bin is therefore judged only against the other townships of
# its own quarter, which is the comparison the Gi* standardisation already makes.
# =============================================================================

# Observed space-time Gi*.
#   X       n_locs x n_times matrix, rows in the order of `nb_self`
#   nb_self self-inclusive neighbour list (spdep / sfdep `nb`)
#   k       number of earlier periods pooled into each bin's neighbourhood
# Weights are row-standardised within each period copy, as in sfdep, so a bin
# with k earlier periods available has total weight k + 1.
gi_spt <- function(X, nb_self, k = 1L) {
  n  <- nrow(X); TT <- ncol(X)
  w0 <- 1 / lengths(nb_self)
  G  <- matrix(NA_real_, n, TT)
  for (t in seq_len(TT)) {
    sl  <- max(1L, t - k):t
    S   <- rowSums(X[, sl, drop = FALSE])
    lx  <- w0 * vapply(nb_self, function(j) sum(S[j]), numeric(1))
    m   <- length(sl)
    x   <- X[, t]
    Wi  <- m
    S1i <- m * w0
    si2 <- sum((x - mean(x))^2) / n
    VG  <- si2 * ((n * S1i - Wi^2) / (n - 1))
    G[, t] <- (lx - Wi * mean(x)) / sqrt(VG)
  }
  G
}

# Period-conditional permutation p-values for the same statistic. Because the
# expectation and variance of Gi* depend only on the weights and on the
# period's own mean and variance, ranking the neighbourhood sum is equivalent
# to ranking Gi* itself, and is cheaper.
# Returns both the folded p-value that sfdep thresholds (min of the two tails)
# and the two-sided p-value, which is twice it.
gi_spt_perm <- function(X, nb_self, k = 1L, nsim = 499L, seed = 1234L) {
  set.seed(seed)
  n  <- nrow(X); TT <- ncol(X)
  w0 <- 1 / lengths(nb_self)
  cards <- lengths(nb_self) - 1L
  self_in <- vapply(seq_len(n), function(i) i %in% nb_self[[i]], logical(1))
  stopifnot(all(self_in))
  p_fold <- matrix(NA_real_, n, TT)
  for (t in seq_len(TT)) {
    sl <- max(1L, t - k):t
    S  <- rowSums(X[, sl, drop = FALSE])
    for (i in seq_len(n)) {
      obs <- sum(S[nb_self[[i]]])
      ci  <- cards[i]
      r   <- replicate(nsim, sample.int(n - 1L, ci))
      r   <- matrix(r, nrow = ci)
      idx <- r + (r >= i)                      # skip the township itself
      sim <- S[i] + colSums(matrix(S[idx], nrow = ci))
      l <- (sum(obs >= sim) + 1) / (nsim + 1)
      g <- (sum(obs <= sim) + 1) / (nsim + 1)
      p_fold[i, t] <- min(l, g)
    }
  }
  # pmin() takes its attributes from its first argument, so the matrix is
  # filled in place rather than rebuilt; otherwise the dimensions are lost.
  p_two <- p_fold
  p_two[] <- pmin(1, 2 * p_fold)
  list(p_fold = p_fold, p_two = p_two)
}

# The class rules. sfdep's own list is used, evaluated in sfdep's order, with one
# correction. In sfdep 0.2.5 the "intensifying hotspot" rule counts significant
# quarters of EITHER sign and never checks that the last one is hot, and it is
# evaluated before "diminishing coldspot". A township that is a significant cold
# spot in 90% of quarters and whose Gi* is rising therefore comes out as an
# intensifying HOT spot. The corrected rule counts only significant hot quarters
# and requires the last one to be hot, which is how sfdep already writes the
# cold version and how Esri defines the class.
ehsa_rules <- function(corrected = TRUE) {
  fx <- sfdep:::fxs
  if (corrected) {
    fx[["intensifying hotspot"]] <- function(gs, sigs, n, tau, tau_p, threshold, ...) {
      (sum(sigs & gs > 0) / n >= 0.9) && sigs[n] && gs[n] > 0 &&
        (tau > 0) && (tau_p < threshold)
    }
  }
  fx
}

# Mann-Kendall on each township's Gi* series, then the class rules. With
# corrected = FALSE this reproduces sfdep:::classify_hotspot() exactly, which is
# checked in scripts/04_ehsa.R, so a difference in the result can come only from
# the p-values fed in or from the one rule that is corrected.
# sfdep spells "oscilating" with one l; it is corrected here, once, for display.
ehsa_classify <- function(G, P, threshold = 0.05, locs, corrected = TRUE) {
  fx <- ehsa_rules(corrected)
  out <- lapply(seq_len(nrow(G)), function(i) {
    gs <- G[i, ]; sigs <- P[i, ] <= threshold; n <- length(gs)
    mk <- Kendall::MannKendall(gs)
    tau <- unname(mk$tau); tau_p <- as.numeric(mk$sl)
    hit <- vapply(fx, function(f) isTRUE(f(gs = gs, sigs = sigs, n = n, tau = tau,
                                            tau_p = tau_p, threshold = threshold)),
                  logical(1))
    data.frame(location = locs[i], tau = tau, p_value = tau_p,
               classification = names(fx)[hit][1])
  })
  res <- do.call(rbind, out)
  res$classification <- sub("oscilating", "oscillating", res$classification)
  res
}

# Esri-style display labels, in the order used by PAL_EHSA in R/utils.R.
ehsa_label <- function(x) {
  x <- sub("hotspot$", "hot spot", sub("coldspot$", "cold spot", x))
  x <- paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  factor(x, levels = names(PAL_EHSA))
}

# Panel (long, ordered by period then location) to an n_locs x n_times matrix.
panel_matrix <- function(panel, var, locs) {
  periods <- sort(unique(panel$period))
  M <- matrix(NA_real_, length(locs), length(periods),
              dimnames = list(locs, as.character(periods)))
  M[cbind(match(panel$TS_PCODE, locs), match(panel$period, periods))] <- panel[[var]]
  stopifnot(!anyNA(M))
  M
}
