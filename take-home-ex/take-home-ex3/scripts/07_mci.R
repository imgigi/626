# =============================================================================
# 07_mci.R -- value-added analysis: where the riders go, and what that loads
# Take-home Exercise 3 | ISSS626
#
# The Huff model says where a trip enters the network and nothing about where it goes.
# A multiplicative competitive interaction model fitted to the observed OD shares
# supplies the second half: given a trip starts at station i, which destination does
# it choose? Combined with the Huff boardings and a rail graph, it shows what happens
# after the gate: how many riders stay on the new line, where the rest change, and
# which existing links gain or lose. Falling tap-ins at an interchange can hide rising
# transfers there, which a station-level model cannot see.
# Run from take-home-ex/take-home-ex3/.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(igraph); library(MCI)
})
source("R/utils.R"); source("R/rail.R")
DER <- "data/derived"
set.seed(1234)

codes <- readRDS(file.path(DER, "station_codes.rds"))
stn   <- readRDS(file.path(DER, "stations.rds"))
dE    <- readRDS(file.path(DER, "dest_existing.rds"))
dP    <- readRDS(file.path(DER, "dest_proposed.rds"))
odst  <- readRDS(file.path(DER, "od_station.rds"))
SC    <- readRDS(file.path(DER, "scenario.rds"))
MD    <- readRDS(file.path(DER, "model.rds"))

# ---- The rail graph ------------------------------------------------------------
# Built by R/rail.R; the running speeds and the transfer penalty are assumptions, and the
# penalty is varied below.
cx <- codes |>
  left_join(select(stn, station_id, x, y), by = "station_id") |>
  select(code = station_code, station_id, station_name, x, y)
stopifnot(all(unlist(map(c(NETWORK, PROPOSED), 1)) %in% cx$code))

g_old <- build_graph(NETWORK, cx)
g_new <- build_graph(c(NETWORK, PROPOSED), cx)
id_old <- dE$destination; id_new <- dP$destination
stopifnot(all(paste0("IN_", id_old) %in% V(g_old)$name),
          all(paste0("IN_", id_new) %in% V(g_new)$name))
RM_old <- rail_matrix(g_old, id_old)
RM_new <- rail_matrix(g_new, id_new)
stopifnot(all(is.finite(RM_old)), all(is.finite(RM_new)))

# ---- The MCI model on observed OD shares -----------------------------------------
# log-centring: subtract, within each origin, the mean log of each variable. OLS on the
# result estimates the MCI exponents directly. mci.fit() needs a complete rectangular
# table and stops on this one, because 1,466 of the 33,672 station pairs recorded no
# trips, so the transformation is done here and verified against the package below on
# a complete block of the table.
vars_of <- function(d) {
  d |> transmute(destination,
                 v_bus = 1 + bus_400, v_ret = 1 + retail_400, v_food = 1 + food_400,
                 v_civ = 1 + civic_400, v_lei = 1 + leisure_400,
                 v_work = 1 + work_m2_400 / 1000, v_res = 1 + res_800 / 1000,
                 v_lines = 1 + ic, lrt = lrt)
}
Z_old <- vars_of(dE)
LOG_V <- c("v_bus", "v_ret", "v_food", "v_civ", "v_lei", "v_work", "v_res", "v_lines")

long_rail <- function(RM) {
  tibble(o = rep(rownames(RM), times = ncol(RM)), d = rep(colnames(RM), each = nrow(RM)),
         minutes = as.vector(RM))
}
pair_tbl <- odst |>
  left_join(long_rail(RM_old), by = c("o", "d")) |>
  left_join(Z_old, by = c("d" = "destination")) |>
  left_join(select(dE, o = destination, ox = x, oy = y), by = "o") |>
  left_join(select(dE, d = destination, dx = x, dy = y), by = "d") |>
  mutate(km = pmax(sqrt((ox - dx)^2 + (oy - dy)^2) / 1000, 0.1)) |>
  group_by(o) |>
  mutate(share = daily / sum(daily)) |>
  ungroup()

log_centre <- function(df, v) log(df[[v]]) - ave(log(df[[v]]), df$o)
design <- function(df, dist = "minutes", attr = LOG_V, use_lrt = TRUE) {
  X <- sapply(c(dist, attr), \(v) log_centre(df, v))
  if (use_lrt) X <- cbind(X, lrt = df$lrt - ave(df$lrt, df$o))
  X
}
fit_mci <- function(df, dist = "minutes", attr = LOG_V, use_lrt = TRUE, wls = FALSE) {
  y <- log_centre(df, "share")
  X <- design(df, dist, attr, use_lrt)
  if (wls) lm(y ~ 0 + X, weights = df$daily) else lm(y ~ 0 + X)
}
m_full   <- fit_mci(pair_tbl)
m_wls    <- fit_mci(pair_tbl, wls = TRUE)
m_euclid <- fit_mci(pair_tbl, dist = "km")
m_noattr <- lm(log_centre(pair_tbl, "share") ~ 0 + log_centre(pair_tbl, "minutes"))
m_nodist <- local({
  y <- log_centre(pair_tbl, "share")
  X <- cbind(sapply(LOG_V, \(v) log_centre(pair_tbl, v)), lrt = pair_tbl$lrt - ave(pair_tbl$lrt, pair_tbl$o))
  lm(y ~ 0 + X)
})
mci_null <- tibble(model = c("Rail time only", "Attraction only", "Euclidean distance + attraction", "Rail time + attraction"),
                   r2 = c(summary(m_noattr)$r.squared, summary(m_nodist)$r.squared,
                          summary(m_euclid)$r.squared, summary(m_full)$r.squared))
r2_wls <- summary(m_wls)$r.squared

# Verification: the package and the hand transformation on a complete block.
big <- dE |> arrange(desc(daily_in)) |> pull(destination)
blk <- pair_tbl |> filter(o %in% big[1:25], d %in% big[26:50]) |>
  group_by(o) |> mutate(share = daily / sum(daily)) |> ungroup()
stopifnot(nrow(blk) == 625)
pk <- mci.fit(as.data.frame(blk |> select(o, d, share, minutes, v_bus, v_ret, v_food)), "o", "d", "share",
              "minutes", "v_bus", "v_ret", "v_food")
hd <- lm(log_centre(blk, "share") ~ 0 + I(log_centre(blk, "minutes")) + I(log_centre(blk, "v_bus")) +
           I(log_centre(blk, "v_ret")) + I(log_centre(blk, "v_food")))
mci_check <- max(abs(unname(coef(pk)) - unname(coef(hd))))

# Log-centred OLS, the standard estimator, treats a pair with 0.02 trips a day like one
# with 5,000, and the many rarely used pairs pull the distance exponent too steep: trips
# come out far too short. Weighting pairs by their trips helps little and damages the fit
# to station arrivals. So the attraction estimates are kept and the distance exponent
# alone is recalibrated to the observed mean trip length, the classical calibration of a
# trip distribution model (Hyman, 1969). This is done below, once the observed flows are
# in place.
clean <- function(m) { b <- coef(m); names(b) <- sub("^X", "", names(b)); b }
beta_ols <- clean(m_full); beta_wls <- clean(m_wls)
coef_mci <- tibble(term = names(beta_ols), ols = unname(beta_ols), ols_se = summary(m_full)$coefficients[, 2],
                   wls = unname(beta_wls), wls_se = summary(m_wls)$coefficients[, 2])
beta <- beta_ols

# Log-centred OLS treats a pair with 0.02 trips a day like one with 5,000. The cluster of
# rarely used pairs is mostly noise, so the fit is also read where it matters: at the
# level of each station's total arrivals, and on the pairs that carry real volume.
share_pred <- function(RM, Z, delta = NULL, gam = rep(beta["minutes"], nrow(RM))) {
  lm_ <- gam * log(RM)
  za <- sapply(LOG_V, \(v) beta[v] * log(Z[[v]])) |> rowSums() + beta["lrt"] * Z$lrt
  if (!is.null(delta)) za <- za + delta
  u <- sweep(lm_, 2, za, "+")
  diag(u) <- -Inf
  mx <- apply(u, 1, max)
  v <- exp(u - mx)
  v / rowSums(v)
}
O_obs <- pair_tbl |> group_by(o) |> summarise(O = sum(daily), .groups = "drop")
O_vec <- setNames(O_obs$O, O_obs$o)[id_old]
D_obs <- pair_tbl |> group_by(d) |> summarise(A = sum(daily), .groups = "drop")
A_vec <- setNames(D_obs$A, D_obs$d)[id_old]

obs_mat <- matrix(0, length(id_old), length(id_old), dimnames = list(id_old, id_old))
obs_mat[cbind(match(pair_tbl$o, id_old), match(pair_tbl$d, id_old))] <- pair_tbl$daily
TRIP_OBS <- sum(obs_mat * RM_old) / sum(obs_mat)

# One exponent fits the network's mean trip but not its parts: it overstates trips from the
# periphery and understates them from the centre, because people far out do not make
# proportionately longer trips. All the proposed stations are on the periphery, so each
# ring of rail time from the CBD gets its own exponent, calibrated to its own observed mean.
RINGS <- c(0, 15, 30, Inf)
ring_of <- function(m) as.integer(cut(m, RINGS))
ring_old <- ring_of(dE$min_cbd)

mean_trip <- function(gam, rows = TRUE) {
  f <- O_vec * share_pred(RM_old, Z_old, gam = gam)
  sum(f[rows, ] * RM_old[rows, ]) / sum(f[rows, ])
}
calibrate_rings <- function() {
  sapply(sort(unique(ring_old)), function(k) {
    rows <- ring_old == k
    target <- sum(obs_mat[rows, ] * RM_old[rows, ]) / sum(obs_mat[rows, ])
    uniroot(\(g) mean_trip(rep(g, nrow(RM_old)), rows) - target, c(-3.5, -0.05), tol = 1e-4)$root
  })
}

beta <- beta_ols
trip_ols <- mean_trip(rep(beta_ols["minutes"], nrow(RM_old)))
beta <- beta_wls
trip_wls <- mean_trip(rep(beta_wls["minutes"], nrow(RM_old)))
beta <- beta_ols
one_gamma <- uniroot(\(g) mean_trip(rep(g, nrow(RM_old))) - TRIP_OBS, c(-3.5, -0.05), tol = 1e-4)$root
trip_one <- function(k) { rows <- ring_old == k; mean_trip(rep(one_gamma, nrow(RM_old)), rows) }
GAM <- calibrate_rings()
gam_old <- GAM[ring_old]

P0 <- share_pred(RM_old, Z_old, gam = gam_old)
F0 <- O_vec * P0
fit_arrivals <- tibble(r2_log = 1 - sum((log(A_vec) - log(colSums(F0)))^2) / sum((log(A_vec) - mean(log(A_vec)))^2),
                       rho = cor(A_vec, colSums(F0), method = "spearman"))
big_pairs <- obs_mat >= 10
fit_pairs <- tibble(n = sum(big_pairs), share_trips = sum(obs_mat[big_pairs]) / sum(obs_mat),
                    rho = cor(obs_mat[big_pairs], F0[big_pairs], method = "spearman"),
                    r2_log = 1 - sum((log(obs_mat[big_pairs]) - log(F0[big_pairs]))^2) /
                      sum((log(obs_mat[big_pairs]) - mean(log(obs_mat[big_pairs])))^2))

trip_check <- tibble(observed = TRIP_OBS, ols = trip_ols, wls = trip_wls, one_exponent = TRIP_OBS,
                     rings = sum(F0 * RM_old) / sum(F0))
trip_rings <- tibble(ring = c("Under 15 min", "15 to 30 min", "Over 30 min"), n = as.integer(table(ring_old)),
                     exponent = GAM,
                     observed = sapply(1:3, \(k) { r <- ring_old == k; sum(obs_mat[r, ] * RM_old[r, ]) / sum(obs_mat[r, ]) }),
                     one_exponent = sapply(1:3, trip_one),
                     rings = sapply(1:3, \(k) { r <- ring_old == k; sum(F0[r, ] * RM_old[r, ]) / sum(F0[r, ]) }))
trip_by_origin <- tibble(destination = id_old, station_name = dE$station_name, min_cbd = dE$min_cbd, O = O_vec,
                         obs = rowSums(obs_mat * RM_old) / rowSums(obs_mat),
                         mod = rowSums(F0 * RM_old) / rowSums(F0)) |>
  mutate(gap = mod - obs)

# Hold out each line again, now as a destination: its stations are left out of the share
# regression, and the arrivals predicted for them are compared with what they recorded.
folds <- setNames(MD$folds, id_old)
lolo_arr <- rep(NA_real_, length(id_old)); names(lolo_arr) <- id_old
for (f in setdiff(unique(folds), NA)) {
  te <- names(folds)[which(folds == f)]
  sub <- pair_tbl |> filter(!o %in% te, !d %in% te) |> group_by(o) |> mutate(share = daily / sum(daily)) |> ungroup()
  beta <- clean(fit_mci(sub))
  Pf <- share_pred(RM_old, Z_old, gam = calibrate_rings()[ring_old])
  lolo_arr[te] <- colSums(O_vec * Pf)[te]
}
beta <- beta_ols
ok <- !is.na(lolo_arr)
fit_lolo <- tibble(r2_log = 1 - sum((log(A_vec[ok]) - log(lolo_arr[ok]))^2) / sum((log(A_vec[ok]) - mean(log(A_vec[ok])))^2),
                   rho = cor(A_vec[ok], lolo_arr[ok], method = "spearman"),
                   mape = median(abs(lolo_arr[ok] / A_vec[ok] - 1)))

# ---- Destination flows on the proposed network -------------------------------------
# An existing destination keeps its own residual: Raffles Place draws more than its
# variables say and is assumed to go on doing so. Existing pairs are pivoted from the
# observed flow, so only the modelled change moves them; pairs that involve a new
# station come from the model alone.
delta_old <- log(pmax(A_vec, 1)) - log(colSums(F0))
Z_new <- vars_of(dP)
delta_new <- setNames(numeric(length(id_new)), id_new)
delta_new[id_old] <- delta_old

# A proposed station takes the exponent of its own ring; an existing one keeps its own,
# so that the pivot compares like with like.
ring_new <- ring_of(dP$min_cbd)
ring_new[match(id_old, id_new)] <- ring_old
gam_new <- GAM[ring_new]
P0d <- share_pred(RM_old, Z_old, delta = delta_old, gam = gam_old)
P1d <- share_pred(RM_new, Z_new, delta = delta_new, gam = gam_new)
# Boardings on the proposed network, from the closed Huff model. OD totals run a little
# under tap-ins, so each station's OD total is scaled by the same proportion.
od_ratio <- sum(O_vec) / sum(dE$daily_in)
new_board <- bind_rows(
  transmute(SC$old_tbl, destination, O_closed = E_closed * (O_vec[destination] / daily_in),
            O_open = E_open * (O_vec[destination] / daily_in)),
  transmute(SC$new_tbl, destination, O_closed = E_closed * od_ratio, O_open = E_open * od_ratio))
O_new <- setNames(new_board$O_closed, new_board$destination)[id_new]
O_new_open <- setNames(new_board$O_open, new_board$destination)[id_new]

flow_new <- function(O) {
  Fm <- O * P1d
  ii <- match(id_old, id_new)
  # Observed flow, moved by how the model's flow changes between the two networks.
  base_model <- O_vec * P0d
  ratio <- Fm[ii, ii] / pmax(base_model, 1e-9)
  ratio[!is.finite(ratio)] <- 1
  Fm[ii, ii] <- ifelse(base_model > 0, obs_mat * ratio, 0)
  # Scale each origin back to its boardings, because observed pairs with no trips are 0.
  Fm * (O / pmax(rowSums(Fm), 1e-9))
}
# Doubly constraining the flows. Left alone, a destination-choice model treats each new
# station as an extra place to go and draws trips to it from everywhere, when the person who
# now lives nearer Gek Poh simply stops arriving at Boon Lay. Stations balance: tap-ins and
# tap-outs are almost equal at each of them over a day. So each station's arrivals are made
# to move in step with its boardings, and the matrix is balanced to both (Furness).
in_out_cor <- cor(log(dE$daily_in), log(dE$daily_out))
furness <- function(Fm, O, A, its = 50) {
  A <- A * sum(O) / sum(A)
  for (k in seq_len(its)) {
    Fm <- Fm * (O / pmax(rowSums(Fm), 1e-9))
    Fm <- sweep(Fm, 2, A / pmax(colSums(Fm), 1e-9), "*")
  }
  Fm
}
arrivals_for <- function(O) {
  A <- setNames(O, id_new)                      # a new station balances
  A[id_old] <- A_vec[id_old] * (O[match(id_old, id_new)] / O_vec[id_old])
  A
}
F_new <- furness(flow_new(O_new), O_new, arrivals_for(O_new))
F_new_open <- furness(flow_new(O_new_open), O_new_open, arrivals_for(O_new_open))
balance_gap <- max(abs(rowSums(F_new) / O_new - 1))

# ---- Where the new lines' riders go --------------------------------------------------
new_ids <- SC$new_tbl$destination
proj <- setNames(SC$new_tbl$project, SC$new_tbl$destination)
on_line <- function(p) {
  # A new line's stations, with the existing interchanges it reaches.
  ids <- names(proj)[proj == p]
  extra <- if (p == "JRL") dP$destination[dP$station_name %in% c("Choa Chu Kang", "Boon Lay", "Jurong East")]
           else dP$destination[dP$station_name %in% c("Pasir Ris", "Hougang", "Ang Mo Kio", "Bright Hill")]
  list(new = ids, all = c(ids, extra), hubs = extra)
}
line_summary <- map_dfr(c("JRL", "CRL Phase 1"), function(p) {
  L <- on_line(p)
  fl <- F_new[L$new, , drop = FALSE]
  tibble(project = p, boardings = sum(fl),
         stays_on_line = sum(fl[, L$all]) / sum(fl),
         to_interchange_itself = sum(fl[, L$hubs]) / sum(fl),
         beyond_line = 1 - sum(fl[, L$all]) / sum(fl))
})
top_dest <- map_dfr(c("JRL", "CRL Phase 1"), function(p) {
  L <- on_line(p)
  fl <- colSums(F_new[L$new, , drop = FALSE])
  tibble(project = p, destination = names(fl), trips = fl, share = fl / sum(fl)) |>
    slice_max(trips, n = 8)
}) |>
  left_join(select(dP, destination, station_name), by = "destination")

# ---- All-or-nothing assignment -------------------------------------------------------
# Each flow takes the fastest path. No crowding, no route choice, and no link counts to
# validate against, so the loads are a comparison between the two networks and not a
# forecast of crowding on any link.
assign_flows <- function(g, Fm, ids) {
  en <- E(g)$type; el <- E(g)$line
  ends_e <- ends(g, E(g))
  load <- numeric(ecount(g))
  w <- E(g)$minutes
  for (o in ids) {
    fl <- Fm[o, ids]
    keep <- which(fl > 0)
    if (!length(keep)) next
    sp <- shortest_paths(g, from = paste0("IN_", o), to = paste0("OUT_", ids[keep]), mode = "out",
                         weights = w, output = "epath")$epath
    for (k in seq_along(keep)) {
      idx <- as.integer(sp[[k]])
      load[idx] <- load[idx] + fl[keep[k]]
    }
  }
  tibble(from = ends_e[, 1], to = ends_e[, 2], line = el, type = en, km = E(g)$km, minutes = w, load = load)
}
as_old <- assign_flows(g_old, obs_mat, id_old)
as_new <- assign_flows(g_new, F_new, id_new)
as_new_open <- assign_flows(g_new, F_new_open, id_new)

station_of <- setNames(cx$station_id, cx$code)
xfer_at <- function(a) {
  a |> filter(type == "transfer") |>
    mutate(station_id = station_of[from]) |>
    group_by(station_id) |>
    summarise(transfers = sum(load), .groups = "drop")
}
xfer <- full_join(rename(xfer_at(as_old), base = transfers),
                  rename(xfer_at(as_new), after = transfers), by = "station_id") |>
  left_join(rename(xfer_at(as_new_open), after_open = transfers), by = "station_id") |>
  mutate(across(c(base, after, after_open), \(v) coalesce(v, 0)), change = after - base) |>
  left_join(select(stn, station_id, station_name), by = "station_id") |>
  arrange(desc(change))

und <- function(a) {
  a |> filter(type == "ride") |>
    mutate(a_ = pmin(from, to), b_ = pmax(from, to)) |>
    group_by(line, a_, b_) |>
    summarise(load = max(sum(load[from == a_]), sum(load[from == b_])),
              both = sum(load), km = first(km), .groups = "drop")
}
links <- full_join(rename(und(as_old), base = load, base_both = both),
                   rename(und(as_new), after = load, after_both = both) |> select(-km),
                   by = c("line", "a_", "b_")) |>
  mutate(across(c(base, after, base_both, after_both), \(v) coalesce(v, 0)),
         km = coalesce(km, und(as_new)$km[match(paste(line, a_, b_), paste(und(as_new)$line, und(as_new)$a_, und(as_new)$b_))]),
         change = after - base, pct = ifelse(base > 0, after / base - 1, NA_real_)) |>
  left_join(select(cx, code, ax = x, ay = y, a_name = station_name), by = c("a_" = "code")) |>
  left_join(select(cx, code, bx = x, by = y, b_name = station_name), by = c("b_" = "code"))

line_km <- links |>
  group_by(line) |>
  summarise(base_pkm = sum(base_both * km), after_pkm = sum(after_both * km), .groups = "drop") |>
  mutate(change = after_pkm / base_pkm - 1)
pkm <- list(base = sum(line_km$base_pkm), after = sum(line_km$after_pkm))

# How much the headline depends on the transfer penalty.
g_new10 <- build_graph(c(NETWORK, PROPOSED), cx, transfer = 10)
g_old10 <- build_graph(NETWORK, cx, transfer = 10)
as_old10 <- assign_flows(g_old10, obs_mat, id_old)
as_new10 <- assign_flows(g_new10, F_new, id_new)
x10 <- full_join(rename(xfer_at(as_old10), b10 = transfers), rename(xfer_at(as_new10), a10 = transfers), by = "station_id") |>
  mutate(across(c(b10, a10), \(v) coalesce(v, 0)), change10 = a10 - b10)
xfer <- left_join(xfer, select(x10, station_id, change10), by = "station_id")

print(as.data.frame(mutate(mci_null, r2 = round(r2, 3)))); print(trip_check)
print(as.data.frame(trip_by_origin |> filter(station_name %in% c("Pioneer", "Lakeside", "Boon Lay", "Chinese Garden", "Joo Koon", "Bukit Batok", "Bukit Gombak", "Choa Chu Kang", "Jurong East", "Yew Tee", "Hougang", "Ang Mo Kio", "Serangoon", "Tampines", "Pasir Ris", "Kovan", "Raffles Place", "Bishan", "Woodlands", "Punggol")) |> transmute(station_name, mins = round(min_cbd), obs = round(obs, 1), mod = round(mod, 1), gap = round(gap, 1))))
print(as.data.frame(mutate(trip_rings, across(where(is.numeric), \(v) round(v, 2)))))
print(trip_by_origin |> mutate(band = cut(min_cbd, c(0, 15, 25, 35, 45, 90))) |> group_by(band) |> summarise(n = n(), obs = weighted.mean(obs, O), mod = weighted.mean(mod, O)))
new_trip <- sum(F_new[new_ids, ] * RM_new[new_ids, ]) / sum(F_new[new_ids, ])
cat("mean rail time from new stations:", new_trip, " in/out cor:", in_out_cor, " balance gap:", balance_gap, "
")
print(as.data.frame(mutate(coef_mci, across(where(is.numeric), \(v) round(v, 3)))))
cat("package check, max coef gap:", mci_check, "\n"); print(fit_arrivals); print(fit_pairs); print(fit_lolo)
print(line_summary); print(as.data.frame(top_dest |> transmute(project, station_name, trips = round(trips), share = round(share, 3))))
print(as.data.frame(xfer |> transmute(station_name, base = round(base), after = round(after), change = round(change), change10 = round(change10)) |> head(14)))
print(as.data.frame(links |> arrange(change) |> transmute(line, a_name, b_name, base = round(base), after = round(after), pct = round(pct, 3)) |> head(12)))
print(as.data.frame(links |> filter(base > 0) |> arrange(desc(change)) |> transmute(line, a_name, b_name, base = round(base), after = round(after), pct = round(pct, 3)) |> head(10)))
print(line_km); print(pkm); cat("transfers total base", sum(xfer$base), "after", sum(xfer$after), "\n")

saveRDS(list(
  coef_mci = coef_mci, mci_null = mci_null, r2_wls = r2_wls, in_out_cor = in_out_cor, balance_gap = balance_gap, gam = GAM, one_gamma = one_gamma, trip_rings = trip_rings, trip_by_origin = trip_by_origin, new_trip = new_trip, mci_check = mci_check, fit_arrivals = fit_arrivals,
  fit_pairs = fit_pairs, fit_lolo = fit_lolo, line_summary = line_summary, top_dest = top_dest,
  xfer = xfer, links = links, line_km = line_km, pkm = pkm,
  n_pairs = nrow(pair_tbl), n_possible = length(id_old) * (length(id_old) - 1),
  assumptions = list(speed = SPEED, transfer = TRANSFER_MIN),
  trip = trip_check
), file.path(DER, "mci.rds"))
message("ALLDONE")
