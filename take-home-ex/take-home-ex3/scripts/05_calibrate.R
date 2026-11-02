# =============================================================================
# 05_calibrate.R -- calibrating the model, and testing it against simpler ones
# Take-home Exercise 3 | ISSS626
#
# A ladder of models from the nearest-station rule upwards. Each is scored in sample
# and by leaving out one line at a time, which is the closest available imitation of
# the real task: predicting a line the model has never seen, with its stations absent
# from the network during calibration. The selected model is then profiled, tested
# against a permutation null, bootstrapped, and checked against MCI.
# Run from take-home-ex/take-home-ex3/. About 15 minutes on 8 workers.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(spdep); library(sf); library(parallel); library(MCI)
})
source("R/utils.R"); source("R/huff.R")
DER <- "data/derived"
set.seed(1234)
N_WORKERS <- 8
N_BOOT <- 300
DMAX <- 5   # km; the choice set stops here (the profile below shows why)

origins <- readRDS(file.path(DER, "origins.rds"))
dest    <- readRDS(file.path(DER, "dest_existing.rds"))
dP      <- readRDS(file.path(DER, "dest_proposed.rds"))
codes   <- readRDS(file.path(DER, "station_codes.rds"))
D       <- readRDS(file.path(DER, "D_existing.rds"))
Dn      <- readRDS(file.path(DER, "D_proposed.rds"))
AF      <- readRDS(file.path(DER, "audit_frames.rds"))
X_VARS  <- AF$X_VARS
X_BASE  <- setdiff(X_VARS, "l_cbd")
idx_old <- match(dest$destination, dP$destination)
idx_new <- which(dP$role == "new access point")

# Each market column is scaled to the same total as residents, so that its kappa says
# how many tap-ins a resident-equivalent unit of it generates relative to the others.
TOT <- sum(origins$residents)
nz  <- function(v) v / sum(v) * TOT
MH  <- cbind(residents = origins$residents, dorm = nz(origins$dorm_m2))
MW  <- cbind(work = nz(origins$work_m2), poi = nz(origins$poi))
M3  <- cbind(residents = origins$residents, work = nz(origins$work_m2), dorm = nz(origins$dorm_m2))
M1  <- cbind(residents = origins$residents)

mkX <- function(vars, d = dest) {
  if (!length(vars)) return(matrix(0, nrow(d), 0))
  sweep(as.matrix(d[, vars, drop = FALSE]), 2, AF$centre[vars], "-")
}

# ---- Folds ------------------------------------------------------------------------
# A station is held out with a line only if it serves no other line; interchanges stay.
folds <- codes |>
  filter(status == "open", station_id %in% dest$destination) |>
  mutate(g = case_when(line %in% c("EW", "CG") ~ "EW",
                       line %in% c("CC", "CE") ~ "CC",
                       line == "BP" ~ "LRT Bukit Panjang",
                       line %in% c("SE", "SW", "STC") ~ "LRT Sengkang",
                       line %in% c("PE", "PW", "PTC") ~ "LRT Punggol",
                       TRUE ~ line)) |>
  group_by(station_id) |>
  summarise(n = n_distinct(g), fold = first(g), .groups = "drop") |>
  filter(n == 1) |>
  select(destination = station_id, fold)
dest <- left_join(dest, folds, by = "destination")
FOLDS <- sort(setdiff(unique(dest$fold), NA))

make_spec <- function(Tv, M, vars = X_VARS, type = "pow", dmax = DMAX, use_A0 = FALSE) {
  list(D = D, T = Tv, M = M, X = mkX(vars), type = type, dmax = dmax, use_A0 = use_A0)
}

lolo_part <- function(spec, fix_lambda = NULL) {
  pred <- rep(NA_real_, ncol(D))
  for (f in FOLDS) {
    te <- which(dest$fold == f); tr <- setdiff(seq_len(ncol(D)), te)
    s2 <- modifyList(spec, list(D = D[, tr, drop = FALSE], X = spec$X[tr, , drop = FALSE], T = spec$T[tr]))
    fit <- fit_huff(s2, fix_lambda = fix_lambda)
    pred[te] <- expected(fit$par, s2, D = D, X = spec$X)[te]
  }
  pred
}

stats_row <- function(name, group, E, pred, lambda, n_par) {
  ok <- !is.na(pred); T <- dest$daily_in
  tibble(model = name, group = group, lambda = lambda, n_par = n_par,
         insample = fit_stats(T, E)$r2_log, insample_rho = fit_stats(T, E)$rho,
         lolo_r2 = fit_stats(T[ok], pred[ok])$r2_log, lolo_rmsle = fit_stats(T[ok], pred[ok])$rmsle,
         lolo_mape = fit_stats(T[ok], pred[ok])$mape, lolo_rho = fit_stats(T[ok], pred[ok])$rho)
}

# Single-purpose: one model of all tap-ins.
score_one <- function(d) {
  sp <- make_spec(dest$daily_in, d$M, d$vars, d$type %||% "pow", d$dmax %||% Inf, d$use_A0 %||% FALSE)
  fit <- fit_huff(sp, fix_lambda = d$fix)
  pred <- lolo_part(sp, d$fix)
  list(name = d$name, group = d$group, fit = fit, lolo = pred,
       row = stats_row(d$name, d$group, fit$E, pred, fit$par[1], length(fit$par)))
}
# Two-purpose: the home end and the activity end, calibrated separately and summed.
score_two <- function(d) {
  sh <- make_spec(dest$T_H, d$MH %||% MH, d$vars, d$type %||% "pow", d$dmax %||% DMAX, d$use_A0 %||% FALSE)
  sw <- make_spec(dest$T_W, d$MW %||% MW, d$vars, d$type %||% "pow", d$dmax %||% DMAX, FALSE)
  fh <- fit_huff(sh, fix_lambda = d$fix); fw <- fit_huff(sw, fix_lambda = d$fix)
  ph <- lolo_part(sh, d$fix); pw <- lolo_part(sw, d$fix)
  list(name = d$name, group = d$group, H = fh, W = fw, lolo_H = ph, lolo_W = pw, lolo = ph + pw,
       row = stats_row(d$name, d$group, fh$E + fw$E, ph + pw, NA_real_, length(fh$par) + length(fw$par)) |>
         mutate(lambda_H = fh$par[1], lambda_W = fw$par[1]))
}

in_cluster <- function(items, fun) {
  cl <- makeCluster(N_WORKERS)
  on.exit(stopCluster(cl))
  clusterEvalQ(cl, {
    suppressPackageStartupMessages({library(dplyr); library(tibble)})
    source("R/huff.R")
    `%||%` <- function(a, b) if (is.null(a)) b else a
  })
  clusterExport(cl, c("dest", "D", "Dn", "dP", "AF", "MH", "MW", "M1", "M3", "DMAX", "FOLDS", "X_VARS", "X_BASE",
                      "idx_old", "idx_new", "mkX", "make_spec", "lolo_part", "stats_row", "score_one", "score_two"),
                envir = globalenv())
  parLapplyLB(cl, items, fun)
}

# ---- The ladder -------------------------------------------------------------------
ladder <- cache_rds("ladder", {
  entries <- list(
    list(kind = "one", name = "M1 Huff, equal A, lambda = -2", group = "Reference", M = M1, vars = character(0), fix = -2),
    list(kind = "one", name = "M2 Huff, equal A, lambda fitted", group = "Reference", M = M1, vars = character(0)),
    list(kind = "one", name = "M3 + LRT and interchange", group = "Reference", M = M1, vars = c("lrt", "ic")),
    list(kind = "one", name = "M4 + bus stops, retail, food", group = "Attraction", M = M1, vars = X_BASE),
    list(kind = "one", name = "M5 + work and dormitory markets", group = "Market", M = M3, vars = X_BASE),
    list(kind = "two", name = "M6 Two purposes", group = "Two purposes", vars = X_BASE, dmax = Inf),
    list(kind = "two", name = "M7 + rail time to the CBD", group = "Two purposes", vars = X_VARS, dmax = Inf),
    list(kind = "two", name = "M8 + 5 km choice set (selected)", group = "Selected", vars = X_VARS, dmax = DMAX),
    list(kind = "two", name = "M8e Exponential decay", group = "Variant", vars = X_VARS, dmax = DMAX, type = "exp"),
    list(kind = "two", name = "M8f lambda = -2 (Zhu et al.)", group = "Variant", vars = X_VARS, dmax = DMAX, fix = -2),
    list(kind = "two", name = "M9 + outside option at the home end", group = "Open model", vars = X_VARS, dmax = DMAX, use_A0 = TRUE)
  )
  out <- in_cluster(entries, function(d) if (d$kind == "one") score_one(d) else score_two(d))
  setNames(out, c("m1", "m2", "m3", "m4", "m5", "m6", "m7", "m8", "m8e", "m8f", "m9"))
})

# M0: everyone uses the nearest station, scaled to total ridership. It has nothing to
# learn beyond that scale, so its hold-out score is its score on the same stations.
nearest_E <- local({
  nn <- apply(D, 1, which.min)
  e <- tapply(M1[, 1], factor(nn, levels = seq_len(ncol(D))), sum)
  e[is.na(e)] <- 0
  as.numeric(e) * sum(dest$daily_in) / sum(e, na.rm = TRUE)
})
row0 <- local({
  E <- pmax(nearest_E, 1); ok <- !is.na(dest$fold); T <- dest$daily_in
  tibble(model = "M0 Nearest station", group = "Reference", lambda = NA_real_, n_par = 1L,
         insample = fit_stats(T, E)$r2_log, insample_rho = fit_stats(T, E)$rho,
         lolo_r2 = fit_stats(T[ok], E[ok])$r2_log, lolo_rmsle = fit_stats(T[ok], E[ok])$rmsle,
         lolo_mape = fit_stats(T[ok], E[ok])$mape, lolo_rho = fit_stats(T[ok], E[ok])$rho)
})

# The transport planner's usual alternative: a direct ridership regression on the same
# variables plus residents within 800 m, with no competition between stations.
drm <- local({
  df <- dest |> mutate(y = log(daily_in))
  fm <- reformulate(c("l_res8", X_VARS), "y")
  pr <- rep(NA_real_, nrow(df))
  for (f in FOLDS) { te <- which(df$fold == f); pr[te] <- predict(lm(fm, df[-te, ]), df[te, ]) }
  m <- lm(fm, df); ok <- !is.na(pr); T <- df$daily_in
  list(model = m, row = tibble(model = "DRM Direct regression", group = "Reference", lambda = NA_real_,
       n_par = length(coef(m)), insample = fit_stats(T, exp(fitted(m)))$r2_log,
       insample_rho = fit_stats(T, exp(fitted(m)))$rho,
       lolo_r2 = fit_stats(T[ok], exp(pr[ok]))$r2_log, lolo_rmsle = fit_stats(T[ok], exp(pr[ok]))$rmsle,
       lolo_mape = fit_stats(T[ok], exp(pr[ok]))$mape, lolo_rho = fit_stats(T[ok], exp(pr[ok]))$rho))
})

ladder_tbl <- bind_rows(row0, drm$row, map_dfr(ladder, "row"))
print(as.data.frame(ladder_tbl |> mutate(across(where(is.numeric), \(v) round(v, 3)))))

sel  <- ladder$m8
open <- ladder$m9
SP_H <- sel$H$spec; SP_W <- sel$W$spec; SP_O <- open$H$spec
message("selected: lambda H ", round(sel$H$par[1], 3), ", W ", round(sel$W$par[1], 3))

# ---- Fold-level hold-out errors -----------------------------------------------------
fold_tbl <- dest |>
  mutate(pred = sel$lolo, err = log(pred / daily_in)) |>
  filter(!is.na(fold)) |>
  group_by(fold) |>
  summarise(n = n(), median_error = median(err), rmsle = sqrt(mean(err^2)),
            actual = sum(daily_in), predicted = sum(pred),
            rho = cor(daily_in, pred, method = "spearman"), .groups = "drop") |>
  mutate(total_ratio = predicted / actual)
print(as.data.frame(mutate(fold_tbl, across(where(is.numeric), \(v) round(v, 2)))))

# Can the model pick out the lowest-ridership stations before they exist? Scored on
# MRT stations held out with their line. AUC is the chance that a station in the
# lowest fifth is ranked below one that is not.
mrt_held <- dest |> mutate(pred = sel$lolo) |> filter(!is.na(fold), lrt == 0)
low_cut  <- quantile(mrt_held$daily_in, 0.2)
auc_low <- local({
  y <- mrt_held$daily_in <= low_cut
  r <- rank(-mrt_held$pred)
  (sum(r[y]) - sum(y) * (sum(y) + 1) / 2) / (sum(y) * sum(!y))
})
hold_bias <- with(mrt_held, median(log(pred / daily_in)))

# ---- Far-field test -----------------------------------------------------------------
# A new station in Jurong should not take riders from Marina Bay. The change at
# stations more than 8 km from any new access point is a plausibility test that the
# hold-out score cannot supply, because a model fits better by letting CBD stations
# draw on residents across the island.
new_xy <- dP[idx_new, c("x", "y")]
km_new <- vapply(seq_len(nrow(dest)), \(i) min(sqrt((dest$x[i] - new_xy$x)^2 + (dest$y[i] - new_xy$y)^2)) / 1000, 0)

pivot_change <- function(fh, sh, fw, sw, X_n = mkX(X_VARS, dP)) {
  bh <- expected(fh$par, sh); nh <- expected(fh$par, sh, D = Dn, X = X_n)
  bw <- expected(fw$par, sw); nw <- expected(fw$par, sw, D = Dn, X = X_n)
  (dest$T_H * nh[idx_old] / pmax(bh, 1e-6) + dest$T_W * nw[idx_old] / pmax(bw, 1e-6)) / dest$daily_in - 1
}

# ---- Profiles -----------------------------------------------------------------------
prof_dmax <- cache_rds("prof_dmax", {
  items <- as.list(c(3, 4, 5, 6, 8, 12, Inf))
  in_cluster_res <- in_cluster(items, function(dm) {
    sh <- make_spec(dest$T_H, MH, X_VARS, "pow", dm); sw <- make_spec(dest$T_W, MW, X_VARS, "pow", dm)
    fh <- fit_huff(sh); fw <- fit_huff(sw)
    ph <- lolo_part(sh); pw <- lolo_part(sw); p <- ph + pw; ok <- !is.na(p)
    list(dmax = dm, fh = fh, fw = fw, sh = sh, sw = sw,
         insample = fit_stats(dest$daily_in, fh$E + fw$E)$r2_log,
         lolo = fit_stats(dest$daily_in[ok], p[ok])$r2_log, rho = fit_stats(dest$daily_in[ok], p[ok])$rho)
  })
  map_dfr(in_cluster_res, function(z) {
    chg <- pivot_change(z$fh, z$sh, z$fw, z$sw)
    tibble(dmax = z$dmax, lambda_H = z$fh$par[1], lambda_W = z$fw$par[1], insample = z$insample, lolo = z$lolo, rho = z$rho,
           far_median = median(chg[km_new > 8]), far_worst = min(chg[km_new > 8]),
           near_median = median(chg[km_new <= 2]), worst = min(chg))
  })
})
print(as.data.frame(mutate(prof_dmax, across(where(is.numeric), \(v) round(v, 3)))))

prof_lambda <- cache_rds("prof_lambda", {
  grid <- expand.grid(part = c("H", "W"), lambda = c(-0.5, -1, -1.5, -2, -2.5, -3, -3.5), stringsAsFactors = FALSE)
  res <- in_cluster(split(grid, seq_len(nrow(grid))), function(g) {
    sp <- make_spec(if (g$part == "H") dest$T_H else dest$T_W, if (g$part == "H") MH else MW)
    fit <- fit_huff(sp, fix_lambda = g$lambda)
    pred <- lolo_part(sp, g$lambda); ok <- !is.na(pred)
    tibble(part = g$part, lambda = g$lambda, insample = fit_stats(sp$T, fit$E)$r2_log,
           lolo = fit_stats(sp$T[ok], pred[ok])$r2_log)
  })
  bind_rows(res)
})

fit_fixed_A0 <- function(la0, sp, par0) {
  fn <- function(th) loss(c(th, la0), sp)
  gr <- function(th) { g <- grad_loss(c(th, la0), sp); g[-length(g)] }
  o <- optim(par0, fn, gr, method = "BFGS", control = list(maxit = 300, reltol = 1e-8))
  o$par <- c(o$par, la0); o
}
prof_a0 <- cache_rds("prof_a0", {
  grid <- c(-2, -1, 0, 0.5, 1, 1.5, 2, 3)
  X <- mkX(X_VARS)
  map_dfr(grid, function(la0) {
    sp <- modifyList(SP_H, list(use_A0 = TRUE))
    o <- fit_fixed_A0(la0, sp, sel$H$par)
    pred <- rep(NA_real_, ncol(D))
    for (f in FOLDS) {
      te <- which(dest$fold == f); tr <- setdiff(seq_len(ncol(D)), te)
      s2 <- modifyList(sp, list(D = D[, tr, drop = FALSE], X = X[tr, , drop = FALSE], T = dest$T_H[tr]))
      o2 <- fit_fixed_A0(la0, s2, sel$H$par)
      pred[te] <- expected(o2$par, s2, D = D, X = X)[te]
    }
    ok <- !is.na(pred)
    z <- model_parts(o$par, sp)
    Xn <- mkX(X_VARS, dP)
    z2 <- model_parts(o$par, sp, D = Dn, X = Xn)
    tibble(lA0 = la0, A0 = exp(la0), lambda = o$par[1], insample = fit_stats(dest$T_H, z$E)$r2_log,
           lolo = fit_stats(dest$T_H[ok], pred[ok])$r2_log,
           share_outside = weighted.mean(1 - rowSums(z$P), z$mk),
           induced = sum(z2$E) - sum(z$E), new_E = sum(z2$E[idx_new]))
  }) |> mutate(share_induced = induced / new_E)
})
print(as.data.frame(mutate(prof_lambda, across(where(is.numeric), \(v) round(v, 3)))))
print(as.data.frame(mutate(prof_a0, across(where(is.numeric), \(v) round(v, 3)))))

# ---- Permutation null -----------------------------------------------------------------
# Is the gain from the attraction variables more than a flexible model would find in
# noise? Bus, retail and food counts are shuffled among stations of the same mode, so
# that the LRT-MRT contrast survives and only their link to ridership is broken.
obs_r2 <- fit_stats(dest$daily_in, sel$H$E + sel$W$E)$r2_log
perm_null <- cache_rds("perm_null", {
  X <- mkX(X_VARS)
  seeds <- 2000 + seq_len(99)
  unlist(in_cluster(as.list(seeds), function(s) {
    set.seed(s)
    Xb <- mkX(X_VARS)
    vv <- c("l_bus", "l_ret", "l_food")
    for (grp in list(which(dest$lrt == 1), which(dest$lrt == 0))) Xb[grp, vv] <- Xb[sample(grp), vv]
    sh <- modifyList(make_spec(dest$T_H, MH), list(X = Xb)); sw <- modifyList(make_spec(dest$T_W, MW), list(X = Xb))
    fit_stats(dest$daily_in, fit_huff(sh)$E + fit_huff(sw)$E)$r2_log
  }))
})
perm_p <- (1 + sum(perm_null >= obs_r2)) / (1 + length(perm_null))
message("permutation: observed ", round(obs_r2, 3), ", null max ", round(max(perm_null), 3), ", p ", perm_p)

# ---- Bootstrap --------------------------------------------------------------------------
# A Bayesian bootstrap: stations are reweighted rather than resampled, so every station
# stays in the choice set, as it does in reality, and only its say in the calibration
# varies. The same weights go to both ends in a draw.
boot <- cache_rds("boot", {
  draws <- in_cluster(as.list(seq_len(N_BOOT)), function(b) {
    set.seed(5000 + b)
    w <- rexp(ncol(D)); w <- w / mean(w)
    sh <- modifyList(make_spec(dest$T_H, MH), list(w = w)); sw <- modifyList(make_spec(dest$T_W, MW), list(w = w))
    so <- modifyList(make_spec(dest$T_H, MH, use_A0 = TRUE), list(w = w))
    list(H = fit_huff(sh, maxit = 150)$par, W = fit_huff(sw, maxit = 150)$par, O = fit_huff(so, maxit = 150)$par)
  })
  list(H = do.call(rbind, map(draws, "H")), W = do.call(rbind, map(draws, "W")), O = do.call(rbind, map(draws, "O")))
})

ci_tbl <- function(par, b, names) {
  q <- apply(b, 2, quantile, c(.025, .975), na.rm = TRUE)
  tibble(term = names, estimate = par, lo = q[1, ], hi = q[2, ])
}
coef_H <- ci_tbl(sel$H$par, boot$H, c("lambda", X_VARS, paste0("log_kappa_", colnames(MH))))
coef_W <- ci_tbl(sel$W$par, boot$W, c("lambda", X_VARS, paste0("log_kappa_", colnames(MW))))
coef_O <- ci_tbl(open$H$par, boot$O, c("lambda", X_VARS, paste0("log_kappa_", colnames(MH)), "log_A0"))
print(as.data.frame(mutate(coef_H, across(where(is.numeric), \(v) round(v, 3)))))
print(as.data.frame(mutate(coef_W, across(where(is.numeric), \(v) round(v, 3)))))
print(as.data.frame(mutate(coef_O, across(where(is.numeric), \(v) round(v, 3)))))

# ---- Exact calibration, for comparison --------------------------------------------------
# What huff.fit() in MCI does: find the attractiveness that reproduces each station's
# ridership exactly, whatever it takes. Its effect on how far a station reaches is what
# persuaded this report to pivot on observed ridership instead.
catch_km <- function(P, DD, mk) as.vector(colSums(mk * P * DD) / pmax(colSums(mk * P), 1e-9))
sh_inf <- ladder$m7$H$spec; fh_inf <- ladder$m7$H
cal_inf <- calibrate_r(fh_inf$par, sh_inf)
cal_tbl <- dest |>
  transmute(destination, station_name, r = cal_inf$r,
            km_structural = catch_km(model_parts(fh_inf$par, sh_inf)$P, D, model_parts(fh_inf$par, sh_inf)$mk),
            km_calibrated = catch_km(model_parts(cal_inf$par, sh_inf, r = cal_inf$r)$P, D,
                                     model_parts(cal_inf$par, sh_inf, r = cal_inf$r)$mk))

# ---- Residual structure ------------------------------------------------------------------
st_pts <- st_as_sf(dest, coords = c("x", "y"), crs = 3414)
lw <- nb2listw(knn2nb(knearneigh(st_pts, k = 6)), style = "W")
E_sel <- sel$H$E + sel$W$E
resid_log <- log(dest$daily_in) - log(E_sel)
moran <- moran.mc(resid_log, lw, nsim = 999)
residuals_tbl <- dest |>
  transmute(destination, station_name, x, y, lrt, daily_in, phi, E = E_sel, E_H = sel$H$E, E_W = sel$W$E,
            resid = resid_log, lolo = sel$lolo, fold)

# ---- Against MCI ---------------------------------------------------------------------------
# huff.shares() and shares.total() on the interaction frame, using the selected home-end
# model's attractiveness and decay. Both functions need a complete rectangular table, so
# the frame cannot simply be cut at 5 km: a station out of range is moved to a distance so
# great that its share vanishes, and hexagons with no station in range, which send no one
# in the matrix form, are left out.
mci_check <- cache_rds("mci_check", {
  si <- readRDS(file.path(DER, "si_existing.rds"))
  q <- unpack(sel$H$par, SP_H)
  att <- tibble(destination = dest$destination, attraction = exp(as.vector(SP_H$X %*% q$beta)))
  mk <- tibble(origin = origins$origin, market = as.vector(SP_H$M %*% q$kappa))
  has <- si |> group_by(origin) |> summarise(any_in = any(distance <= DMAX), .groups = "drop") |> filter(any_in)
  frame <- si |>
    select(origin, destination, distance) |>
    filter(origin %in% has$origin) |>
    mutate(distance = ifelse(distance <= DMAX, distance, 1e9)) |>
    left_join(att, by = "destination") |>
    left_join(mk, by = "origin")
  t0 <- Sys.time()
  sh <- huff.shares(as.data.frame(frame), "origin", "destination", "attraction", "distance",
                    gamma = 1, lambda = q$lambda, atype = "pow", dtype = "pow", check_df = FALSE)
  tot <- shares.total(sh, "origin", "destination", "p_ij", "market", check_df = FALSE)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  e <- tot$sum_E_j[match(dest$destination, tot$suppliers_single)]
  list(max_abs_gap = max(abs(e - sel$H$E)), max_rel_gap = max(abs(e / sel$H$E - 1)), seconds = secs, rows = nrow(frame))
})
print(mci_check)

# ---- Fill the attraction column of the frames ---------------------------------------------
attr_of <- function(d, X) tibble(destination = d$destination,
  A_home = exp(as.vector(X %*% unpack(sel$H$par, SP_H)$beta)),
  A_activity = exp(as.vector(X %*% unpack(sel$W$par, SP_W)$beta)))
si_e <- readRDS(file.path(DER, "si_existing.rds")) |> select(origin, destination, distance) |>
  left_join(attr_of(dest, SP_H$X), by = "destination") |> mutate(in_range = distance <= DMAX)
si_p <- readRDS(file.path(DER, "si_proposed.rds")) |> select(origin, destination, distance) |>
  left_join(attr_of(dP, mkX(X_VARS, dP)), by = "destination") |> mutate(in_range = distance <= DMAX)
saveRDS(si_e, file.path(DER, "si_existing.rds")); saveRDS(si_p, file.path(DER, "si_proposed.rds"))

saveRDS(list(
  DMAX = DMAX, X_VARS = X_VARS,
  ladder_tbl = ladder_tbl, ladder = ladder, drm_coef = coef(drm$model), nearest_E = nearest_E,
  sel = sel, open = open, MH = MH, MW = MW,
  fold_tbl = fold_tbl, auc_low = auc_low, low_cut = low_cut, hold_bias = hold_bias, folds = dest$fold, km_new = km_new,
  prof_dmax = prof_dmax, prof_lambda = prof_lambda, prof_a0 = prof_a0,
  perm_null = perm_null, perm_p = perm_p, obs_r2 = obs_r2,
  boot = boot, coef_H = coef_H, coef_W = coef_W, coef_O = coef_O,
  cal_tbl = cal_tbl, residuals_tbl = residuals_tbl,
  moran = list(I = unname(moran$statistic), p = moran$p.value),
  mci_check = mci_check
), file.path(DER, "model.rds"))
message("ALLDONE")
