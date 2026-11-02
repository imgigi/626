# =============================================================================
# 06_scenario.R -- the proposed network: ridership, impact, and market areas
# Take-home Exercise 3 | ISSS626
#
# A new station has no history, so its ridership comes from the model directly. An
# existing station does, so the model is used only for the ratio between the two
# networks, and that ratio is applied to what the station actually carries today: a
# pivot. Both ends of the trip are carried separately. The closed model, in which every
# new rider comes from an existing station, and the open model, with an outside option,
# are both reported, because the data do not separate them and they differ exactly
# where the question is asked.
# Run from take-home-ex/take-home-ex3/.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(sf)
})
source("R/utils.R"); source("R/huff.R")
DER <- "data/derived"
set.seed(1234)

MD      <- readRDS(file.path(DER, "model.rds"))
origins <- readRDS(file.path(DER, "origins.rds"))
dE      <- readRDS(file.path(DER, "dest_existing.rds"))
dP      <- readRDS(file.path(DER, "dest_proposed.rds"))
codes   <- readRDS(file.path(DER, "station_codes.rds"))
D_old   <- readRDS(file.path(DER, "D_existing.rds"))
D_new   <- readRDS(file.path(DER, "D_proposed.rds"))
AF      <- readRDS(file.path(DER, "audit_frames.rds"))
X_VARS  <- MD$X_VARS
DMAX    <- MD$DMAX
km_new  <- MD$km_new

X_of <- function(d, v = X_VARS) sweep(as.matrix(d[, v, drop = FALSE]), 2, AF$centre[v], "-")
X_new <- X_of(dP)
idx_old <- match(dE$destination, dP$destination)
idx_new <- which(dP$role == "new access point")
stopifnot(!anyNA(idx_old), length(idx_old) + length(idx_new) == nrow(dP))

# ---- Line and stage of each proposed station ----------------------------------------
JRL_S1 <- c(paste0("JS", 1:8), "JW1", "JW2")
JRL_S2 <- paste0("JE", 1:7)
JRL_S3 <- c(paste0("JS", 9:12), paste0("JW", 3:5))
CRL_P1 <- paste0("CR", 2:13)
new_meta <- codes |>
  filter(station_code %in% c(JRL_S1, JRL_S2, JRL_S3, CRL_P1)) |>
  mutate(project = ifelse(station_code %in% CRL_P1, "CRL Phase 1", "JRL"),
         stage = case_when(station_code %in% JRL_S1 ~ "JRL Stage 1 (2028)",
                           station_code %in% JRL_S2 ~ "JRL Stage 2 (2028)",
                           station_code %in% JRL_S3 ~ "JRL Stage 3 (2029)",
                           TRUE ~ "CRL Phase 1 (2030)")) |>
  group_by(destination = station_id) |>
  summarise(code = first(station_code), project = first(project), stage = first(stage), .groups = "drop")

# ---- The pivot ------------------------------------------------------------------------
safe_ratio <- function(n, b) ifelse(b > 1e-6, n / b, 1)

# One network comparison at one end of the trip.
compare <- function(par, spec, X_n = X_new, D_n = D_new, M_n = spec$M, r_new = NULL) {
  base <- model_parts(par, spec)
  r <- NULL
  if (!is.null(r_new)) { r <- numeric(ncol(D_n)); r[idx_new] <- r_new }
  new  <- model_parts(par, spec, D = D_n, X = X_n, M = M_n, r = r)
  list(base = base, new = new, ratio = safe_ratio(new$E[idx_old], base$E))
}
pivot <- function(cH, cW) {
  list(old = dE$T_H * cH$ratio + dE$T_W * cW$ratio,
       new = cH$new$E[idx_new] + cW$new$E[idx_new],
       new_H = cH$new$E[idx_new], new_W = cW$new$E[idx_new])
}

fH <- MD$sel$H; fW <- MD$sel$W; fO <- MD$open$H
SP_H <- fH$spec; SP_W <- fW$spec; SP_O <- fO$spec
cH <- compare(fH$par, SP_H); cW <- compare(fW$par, SP_W); cO <- compare(fO$par, SP_O)
main <- pivot(cH, cW)
open <- pivot(cO, cW)

# ---- Simulation -------------------------------------------------------------------------
# Each draw takes parameters from the bootstrap and gives each new station an unexplained
# pull drawn from the hold-out errors on MRT stations: the model's own record on stations
# it did not see, which carries the bias as well as the spread.
hold_err <- with(MD$residuals_tbl, log(daily_in / lolo)[!is.na(fold) & lrt == 0])

simulate <- function(bH, spH, bW, spW, n = nrow(bH)) {
  old <- matrix(NA_real_, n, length(idx_old)); new <- matrix(NA_real_, n, length(idx_new))
  for (b in seq_len(n)) {
    e <- sample(hold_err, length(idx_new), replace = TRUE)
    a <- pivot(compare(bH[b, ], spH, r_new = e), compare(bW[b, ], spW, r_new = e))
    old[b, ] <- a$old; new[b, ] <- a$new
  }
  list(old = old, new = new)
}
sim_c <- cache_rds("sim_closed", simulate(MD$boot$H, SP_H, MD$boot$W, SP_W))
sim_o <- cache_rds("sim_open", simulate(MD$boot$O, SP_O, MD$boot$W, SP_W))

q_of <- function(m, p) apply(m, 2, quantile, p, na.rm = TRUE)

# Thresholds for a station that would struggle: Hume, the newest station held back for
# lack of demand, and the lowest fifth of MRT stations today.
HUME  <- dE$daily_in[dE$station_name == "Hume"]
LOW20 <- unname(quantile(dE$daily_in[dE$lrt == 0], 0.2))

new_tbl <- dP[idx_new, ] |>
  transmute(destination, station_name, x, y, bus_400, retail_400, food_400, res_800, work_m2_400, min_cbd) |>
  mutate(E_closed = main$new, E_H = main$new_H, E_W = main$new_W,
         lo_closed = q_of(sim_c$new, .05), hi_closed = q_of(sim_c$new, .95),
         E_open = open$new, lo_open = q_of(sim_o$new, .05), hi_open = q_of(sim_o$new, .95),
         p_below_hume = colMeans(sim_c$new < HUME),
         p_below_low20 = colMeans(sim_c$new < LOW20)) |>
  left_join(new_meta, by = "destination") |>
  arrange(desc(E_closed))

old_tbl <- dE |>
  transmute(destination, station_name, x, y, daily_in, T_H, T_W, phi, lrt, role, codes_all, lines_now, lines_proposed) |>
  mutate(km_new = km_new, E_closed = main$old,
         lo_closed = q_of(sim_c$old, .05), hi_closed = q_of(sim_c$old, .95),
         E_open = open$old, lo_open = q_of(sim_o$old, .05), hi_open = q_of(sim_o$old, .95),
         chg_closed = E_closed / daily_in - 1, chg_open = E_open / daily_in - 1,
         chg_H = cH$ratio - 1, chg_W = cW$ratio - 1,
         abs_closed = E_closed - daily_in, abs_open = E_open - daily_in)

# In model units the closed model conserves its market exactly and the open model does
# not, which is how the induced share is measured. In observed units, the pivot makes
# each station's loss proportional to what it really carries, so the sums differ.
induced_model <- sum(cO$new$E) - sum(cO$base$E)
totals <- tibble(
  model = c("Closed", "Open"),
  base = sum(dE$daily_in),
  new_stations = c(sum(main$new), sum(open$new)),
  existing_after = c(sum(main$old), sum(open$old))) |>
  mutate(existing_loss = base - existing_after,
         system_after = existing_after + new_stations,
         net = system_after - base,
         share_new_model = c(0, induced_model / sum(cO$new$E[idx_new] + cW$new$E[idx_new])))
sim_tot <- function(s) list(new = quantile(rowSums(s$new), c(.05, .95)),
                            loss = quantile(sum(dE$daily_in) - rowSums(s$old), c(.05, .95)))
totals_ci <- list(closed = sim_tot(sim_c), open = sim_tot(sim_o))

by_project <- new_tbl |>
  group_by(project) |>
  summarise(n = n(), closed = sum(E_closed), open = sum(E_open),
            median_closed = median(E_closed), median_open = median(E_open),
            home_share = sum(E_H) / sum(E_closed), .groups = "drop")
by_stage <- new_tbl |>
  group_by(stage) |>
  summarise(n = n(), closed = sum(E_closed), open = sum(E_open), .groups = "drop")

by_distance <- old_tbl |>
  mutate(band = cut(km_new, c(0, 1, 2, 3, 5, 8, 12, 50), labels = c("0-1", "1-2", "2-3", "3-5", "5-8", "8-12", "12+"))) |>
  group_by(band) |>
  summarise(n = n(), median_chg = median(chg_closed), median_open = median(chg_open), worst = min(chg_closed),
            loss = sum(daily_in - E_closed), .groups = "drop")

# ---- Market areas -------------------------------------------------------------------------
# Two readings of "market area", both at the home end, where the idea has a physical
# meaning. The dominant station is the one a hexagon is most likely to use, which draws a
# map. The probabilistic market is the share of a hexagon's residents a station can
# expect, which is what the Huff model really says and which registers a loss long
# before a boundary moves.
Rs <- origins$residents
P0 <- cH$base$P; P1 <- cH$new$P
dom <- function(P, ids) { o <- apply(P, 1, which.max); o[rowSums(P) == 0] <- NA; ids[o] }
id_old <- dE$destination; id_new <- dP$destination
hex_impact <- tibble(
  origin = origins$origin, residents = Rs, market = cH$base$mk,
  dominant_old = dom(P0, id_old), dominant_new = dom(P1, id_new),
  p_dom_old = apply(P0, 1, max), p_dom_new = apply(P1, 1, max),
  share_new = rowSums(P1[, idx_new, drop = FALSE]),
  share_new_open = rowSums(cO$new$P[, idx_new, drop = FALSE]),
  share_outside_old = 1 - rowSums(cO$base$P),
  share_outside_new = 1 - rowSums(cO$new$P)) |>
  mutate(switched = coalesce(dominant_old != dominant_new, FALSE),
         to_new = dominant_new %in% id_new[idx_new])

catch_km <- function(P, DD, mk) as.vector(colSums(mk * P * DD) / pmax(colSums(mk * P), 1e-9))
sum_by <- function(key, ids) { v <- tapply(Rs, factor(key, levels = ids), sum); v[is.na(v)] <- 0; as.numeric(v) }
area_tbl <- tibble(destination = id_old,
                   res_dom_old = sum_by(hex_impact$dominant_old, id_old),
                   res_dom_new = sum_by(hex_impact$dominant_new, id_new)[idx_old],
                   res_prob_old = as.vector(crossprod(Rs, P0)),
                   res_prob_new = as.vector(crossprod(Rs, P1[, idx_old])),
                   km_old = catch_km(P0, D_old, cH$base$mk),
                   km_new_net = catch_km(P1[, idx_old], D_new[, idx_old], cH$new$mk)) |>
  mutate(dom_change = res_dom_new - res_dom_old,
         prob_change = res_prob_new / pmax(res_prob_old, 1) - 1)
old_tbl <- left_join(old_tbl, area_tbl, by = "destination")
new_tbl <- new_tbl |>
  left_join(tibble(destination = id_new[idx_new],
                   res_dom = sum_by(hex_impact$dominant_new, id_new)[idx_new],
                   res_prob = as.vector(crossprod(Rs, P1[, idx_new])),
                   km = catch_km(P1[, idx_new, drop = FALSE], D_new[, idx_new, drop = FALSE], cH$new$mk)),
            by = "destination")

coverage <- function(D, km = c(.4, .8, 1.2)) {
  n <- apply(D, 1, min)
  setNames(sapply(km, \(k) sum(Rs[n <= k]) / sum(Rs)), paste0("within_", km))
}
cov_old <- coverage(D_old); cov_new <- coverage(D_new)

# ---- Sensitivity --------------------------------------------------------------------------
jrl_rows <- which(dP$destination[idx_new] %in% new_meta$destination[new_meta$project == "JRL"])
sens_row <- function(name, a) {
  chg <- a$old / dE$daily_in - 1
  tibble(case = name, new_total = sum(a$new), jrl = sum(a$new[jrl_rows]), crl = sum(a$new[-jrl_rows]),
         rho_new = cor(a$new, main$new, method = "spearman"),
         worst_existing = min(chg), near_median = median(chg[km_new <= 2]), far_median = median(chg[km_new > 8]),
         n_over_5 = sum(chg < -0.05), n_below_hume = sum(a$new < HUME))
}
sens <- list(sens_row("Main: two purposes, 5 km choice set", main),
             sens_row("Open model at the home end", open))

two <- function(f, s_only = NULL) pivot(compare(f$H$par, f$H$spec), compare(f$W$par, f$W$spec))
L <- MD$ladder
sens <- c(sens, list(
  sens_row("No maximum access distance", two(L$m7)),
  sens_row("Fixed lambda = -2 (Zhu et al.)", two(L$m8f)),
  sens_row("Exponential decay", two(L$m8e))))

# A single model of all tap-ins, as in the model ladder's M5.
sens <- c(sens, list(local({
  f <- L$m5$fit
  sp <- f$spec
  Xn <- X_of(dP, setdiff(X_VARS, "l_cbd"))
  b <- model_parts(f$par, sp); n <- model_parts(f$par, sp, D = D_new, X = Xn)
  sens_row("One purpose, all tap-ins (M5)", list(old = dE$daily_in * safe_ratio(n$E[idx_old], b$E), new = n$E[idx_new]))
})))

# Refits under other settings.
refit <- function(hs, ws) {
  fh <- fit_huff(hs, par0 = fH$par); fw <- fit_huff(ws, par0 = fW$par)
  list(H = fh, W = fw)
}
sens <- c(sens, list(local({
  z <- refit(modifyList(SP_H, list(dmax = 3)), modifyList(SP_W, list(dmax = 3)))
  # At 3 km a station can have no market in range at one end; the fit needs a positive E.
  sens_row("3 km choice set", two(z))
})))
sens <- c(sens, list(local({
  v <- setdiff(X_VARS, "l_cbd")
  sh <- modifyList(SP_H, list(X = X_of(dE, v))); sw <- modifyList(SP_W, list(X = X_of(dE, v)))
  fh <- fit_huff(sh); fw <- fit_huff(sw); Xn <- X_of(dP, v)
  sens_row("No rail time to the CBD", pivot(compare(fh$par, sh, X_n = Xn), compare(fw$par, sw, X_n = Xn)))
})))

# What huff.fit() does: attractiveness calibrated to reproduce each station exactly, with
# no maximum access distance, so that stations with an unexplained pull reach far.
cal_case <- local({
  sh <- L$m7$H$spec; sw <- L$m7$W$spec
  ch <- calibrate_r(L$m7$H$par, sh); cw <- calibrate_r(L$m7$W$par, sw)
  run <- function(cal, sp) { r <- numeric(ncol(D_new)); r[idx_old] <- cal$r; expected(cal$par, sp, D = D_new, X = X_new, r = r) }
  eh <- run(ch, sh); ew <- run(cw, sw)
  list(old = eh[idx_old] + ew[idx_old], new = eh[idx_new] + ew[idx_new])
})
sens <- c(sens, list(sens_row("Exact calibration, as in huff.fit()", cal_case)))

# The Jurong Region Line is medium capacity. Half the LRT effect is applied to it.
jrl_all <- which(dP$destination %in% new_meta$destination[new_meta$project == "JRL"])
X_jrl <- X_new; X_jrl[jrl_all, "lrt"] <- X_jrl[jrl_all, "lrt"] + 0.5
a_jrl <- pivot(compare(fH$par, SP_H, X_n = X_jrl), compare(fW$par, SP_W, X_n = X_jrl))
sens <- c(sens, list(sens_row("JRL at half the LRT effect", a_jrl)))

# An existing station that gains a line gets no extra pull from it. The interchange
# coefficient comes from comparing stations, not from watching a line arrive, and
# interchanges were sited in town centres, so the estimated uplift may overstate.
d_noup <- dP |> mutate(ic = ifelse(role == "existing, gains a line", pmax(lines_now - 1, 0), ic))
a_noup <- pivot(compare(fH$par, SP_H, X_n = X_of(d_noup)), compare(fW$par, SP_W, X_n = X_of(d_noup)))
sens <- c(sens, list(sens_row("No uplift at interchanges", a_noup)))
old_tbl <- old_tbl |> mutate(E_noup = a_noup$old, chg_noup = E_noup / daily_in - 1)

# Distance floor.
sens <- c(sens, list(local({
  Do <- D_old; Dn <- D_new; Do[Do < .25] <- .25; Dn[Dn < .25] <- .25
  sh <- modifyList(SP_H, list(D = Do)); sw <- modifyList(SP_W, list(D = Do))
  fh <- fit_huff(sh, par0 = fH$par); fw <- fit_huff(sw, par0 = fW$par)
  sens_row("Distance floor 250 m", pivot(compare(fh$par, sh, D_n = Dn), compare(fw$par, sw, D_n = Dn)))
})))

# The grid, rebuilt at 250 m and 1,000 m from the same people and activity points.
people <- readRDS(file.path(DER, "people_points.rds"))
acts   <- readRDS(file.path(DER, "activity_points.rds"))
land   <- readRDS(file.path(DER, "land.rds"))
grid_case <- function(size) {
  sf::sf_use_s2(FALSE)
  g <- st_make_grid(land, cellsize = size, square = FALSE) |> st_sf(geometry = _)
  g$id <- seq_len(nrow(g))
  to_cell <- function(d) {
    st_join(st_as_sf(d, coords = c("x", "y"), crs = 3414, remove = FALSE), g, join = st_intersects) |>
      st_drop_geometry() |>
      filter(!is.na(id))
  }
  a <- to_cell(filter(people, res > 0)) |> group_by(id) |>
    summarise(x = sum(x * res) / sum(res), y = sum(y * res) / sum(res), residents = sum(res), .groups = "drop")
  b <- to_cell(acts) |> group_by(id) |>
    summarise(ax = mean(x), ay = mean(y), work = sum(area[kind == "work"]), dorm = sum(area[kind == "dorm"]),
              poi = sum(kind == "poi"), .groups = "drop")
  m <- full_join(a, b, by = "id") |>
    mutate(x = coalesce(x, ax), y = coalesce(y, ay), across(c(residents, work, dorm, poi), \(v) coalesce(v, 0)))
  tot <- sum(m$residents); nz <- function(v) v / sum(v) * tot
  MHg <- cbind(residents = m$residents, dorm = nz(m$dorm)); MWg <- cbind(work = nz(m$work), poi = nz(m$poi))
  dd <- function(s) { x <- sqrt(outer(m$x, s$x, "-")^2 + outer(m$y, s$y, "-")^2) / 1000; x[x < .1] <- .1; x }
  Do <- dd(dE); Dn <- dd(dP)
  sh <- modifyList(SP_H, list(D = Do, M = MHg)); sw <- modifyList(SP_W, list(D = Do, M = MWg))
  fh <- fit_huff(sh, par0 = fH$par); fw <- fit_huff(sw, par0 = fW$par)
  sens_row(sprintf("Grid %s m", format(size, big.mark = ",")),
           pivot(compare(fh$par, sh, D_n = Dn), compare(fw$par, sw, D_n = Dn))) |>
    mutate(lambda_H = fh$par[1], lambda_W = fw$par[1], r2 = fit_stats(dE$daily_in, fh$E + fw$E)$r2_log, cells = nrow(m))
}
sens <- c(sens, list(grid_case(250), grid_case(1000)))

# ---- The market of 2030 ---------------------------------------------------------------------
# Land use today is the model's only evidence, and the stations held back in the past were
# held back for what stood around them then. Two additions test how much that matters: flats
# completed in 2026, which the base market leaves out, and Tengah at its planned size of about
# 30,000 HDB flats, spread evenly over the town's land.
blocks  <- readRDS(file.path(DER, "hdb_blocks.rds"))
sz      <- readRDS(file.path(DER, "subzone.rds"))
hex_all <- readRDS(file.path(DER, "hex_markets.rds"))
TENGAH_TARGET <- 30000
PPF <- 3
tengah_du  <- sum(blocks$total_dwelling_units[blocks$PA == "TENGAH"])
tengah_add <- max(TENGAH_TARGET - tengah_du, 0) * PPF
tengah_poly <- st_union(sz[sz$PA == "TENGAH", ])
hex_cent <- st_centroid(st_geometry(hex_all))
t_hex <- hex_all$hex_id[lengths(st_intersects(hex_cent, tengah_poly)) > 0]

miss <- setdiff(t_hex, origins$origin)
xy_miss <- st_coordinates(hex_cent[match(miss, hex_all$hex_id)])
o2 <- bind_rows(select(origins, origin, x, y), tibble(origin = miss, x = xy_miss[, 1], y = xy_miss[, 2]))
add26 <- blocks |>
  filter(year_completed == 2026) |>
  st_join(select(hex_all, hex_id), join = st_intersects) |>
  st_drop_geometry() |>
  group_by(hex_id) |>
  summarise(add = sum(res), .groups = "drop")
o2 <- o2 |>
  left_join(add26, by = c("origin" = "hex_id")) |>
  mutate(add26 = coalesce(add, 0), addT = ifelse(origin %in% t_hex, tengah_add / length(t_hex), 0))
Dn2 <- { x <- sqrt(outer(o2$x, dP$x, "-")^2 + outer(o2$y, dP$y, "-")^2) / 1000; x[x < .1] <- .1; x }
pad <- function(M) rbind(M, matrix(0, nrow(o2) - nrow(M), ncol(M)))
MH2 <- pad(SP_H$M); MH2[, "residents"] <- MH2[, "residents"] + o2$add26 + o2$addT
MW2 <- pad(SP_W$M)
a_2030 <- pivot(compare(fH$par, SP_H, D_n = Dn2, M_n = MH2), compare(fW$par, SP_W, D_n = Dn2, M_n = MW2))
sens <- c(sens, list(sens_row("2026 flats and Tengah at 30,000 flats", a_2030)))
new_tbl <- new_tbl |>
  left_join(tibble(destination = dP$destination[idx_new], E_2030 = a_2030$new), by = "destination")

# ---- Committed stations -----------------------------------------------------------------------
# Bedok South, Sungei Bedok and Xilin are due in 2026 and were not open in the reference months.
# Leaving them out isolates the two lines in the brief; this checks it does not change the reading.
stn <- readRDS(file.path(DER, "stations.rds"))
comm <- stn |>
  filter(station_name %in% c("Bedok South", "Sungei Bedok", "Xilin")) |>
  mutate(l_bus = log1p(bus_400), l_ret = log1p(retail_400), l_food = log1p(food_400),
         lrt = 0, ic = ifelse(station_name == "Sungei Bedok", 1, 0),
         l_cbd = log(pmax(sqrt((x - dE$x[dE$station_name == "Raffles Place"])^2 +
                                 (y - dE$y[dE$station_name == "Raffles Place"])^2) / 1000 / 35 * 60, 1)))
dC <- bind_rows(select(dP, station_name, x, y, all_of(X_VARS)), select(comm, station_name, x, y, all_of(X_VARS)))
Dc <- { x <- sqrt(outer(origins$x, dC$x, "-")^2 + outer(origins$y, dC$y, "-")^2) / 1000; x[x < .1] <- .1; x }
a_comm <- pivot(compare(fH$par, SP_H, X_n = X_of(dC), D_n = Dc), compare(fW$par, SP_W, X_n = X_of(dC), D_n = Dc))
sens <- c(sens, list(sens_row("With Bedok South, Sungei Bedok, Xilin", a_comm)))
comm_E <- local({
  eh <- expected(fH$par, SP_H, D = Dc, X = X_of(dC)); ew <- expected(fW$par, SP_W, D = Dc, X = X_of(dC))
  k <- nrow(dP) + seq_len(nrow(comm))
  tibble(station_name = comm$station_name, E = eh[k] + ew[k])
})

sens_tbl <- bind_rows(sens)
print(as.data.frame(mutate(sens_tbl, across(where(is.numeric), \(v) round(v, 3)))))
print(as.data.frame(new_tbl |> transmute(station_name, code, E = round(E_closed), H = round(E_H), W = round(E_W), lo = round(lo_closed), hi = round(hi_closed), O = round(E_open), p_h = round(p_below_hume, 2), p_l = round(p_below_low20, 2), E30 = round(E_2030), mins = round(min_cbd), km = round(km, 1))))
print(as.data.frame(old_tbl |> arrange(chg_closed) |> head(16) |> transmute(station_name, T = round(daily_in), E = round(E_closed), chg = round(chg_closed, 3), lo = round(lo_closed / daily_in - 1, 2), hi = round(hi_closed / daily_in - 1, 2), chg_o = round(chg_open, 3), km = round(km_new, 1), dom = round(dom_change))))
print(as.data.frame(old_tbl |> filter(role == "existing, gains a line") |> transmute(station_name, T = round(daily_in), chg = round(chg_closed, 3), chg_o = round(chg_open, 3), noup = round(chg_noup, 3))))
print(as.data.frame(mutate(by_distance, across(where(is.numeric), \(v) round(v, 3)))))
print(totals); print(by_project); print(by_stage); print(c(HUME = HUME, LOW20 = LOW20)); print(cov_old); print(cov_new); print(comm_E)
cat("switched residents:", sum(Rs[hex_impact$switched]), " to new:", sum(Rs[hex_impact$to_new]), " no station in range:", sum(Rs[is.na(hex_impact$dominant_old)]), "\n")

saveRDS(list(
  new_tbl = new_tbl, old_tbl = old_tbl, hex_impact = hex_impact, totals = totals, totals_ci = totals_ci,
  by_project = by_project, by_stage = by_stage, by_distance = by_distance, cov_old = cov_old, cov_new = cov_new,
  sens_tbl = sens_tbl, HUME = HUME, LOW20 = LOW20, hold_err = hold_err, comm_E = comm_E,
  tengah = list(du_now = tengah_du, target = TENGAH_TARGET, residents_added = tengah_add, n_hex = length(t_hex), ppf = PPF),
  add26_total = sum(o2$add26),
  sim_new_q = list(lo = q_of(sim_c$new, .05), hi = q_of(sim_c$new, .95))
), file.path(DER, "scenario.rds"))
message("ALLDONE")
