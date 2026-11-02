# =============================================================================
# 04_frames.R -- the spatial interaction data frames the brief asks for
# Take-home Exercise 3 | ISSS626
#
# Two interaction frames (existing network; existing plus proposed), one table of
# origins with their local market potential, and one table of destinations with
# observed ridership and attraction variables. The origin and destination columns
# carry the same names in all of them, as MCI's functions require.
# Run from take-home-ex/take-home-ex3/.
# =============================================================================

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(tidyr); library(igraph)
})
source("R/rail.R")
DER <- "data/derived"
FLOOR_KM <- 0.1   # a hexagon's point can sit on a station; decay needs a positive distance

hex <- readRDS(file.path(DER, "hex_markets.rds")) |> st_drop_geometry()
st  <- readRDS(file.path(DER, "stations.rds"))
rid <- readRDS(file.path(DER, "ridership.rds"))

# ---- Origins ------------------------------------------------------------------
# Only hexagons holding some market can send anyone to a station.
origins <- hex |>
  filter(residents > 0 | work_m2 > 0 | dorm_m2 > 0) |>
  transmute(origin = hex_id, x, y, SZ, PA, REGION, residents, res_hdb, res_private,
            hdb_du, poi, work_m2, dorm_m2)

# ---- Destinations -------------------------------------------------------------
# The attraction variables are logged, so that each coefficient reads as an
# elasticity, and centred on the existing network, so that a station at the network's
# average has an attractiveness of one.
transform_x <- function(d) {
  d |> mutate(l_bus  = log1p(bus_400),
              l_ret  = log1p(retail_400),
              l_food = log1p(food_400),
              l_civ  = log1p(civic_400),
              l_lei  = log1p(leisure_400),
              l_work = log1p(work_m2_400 / 1000),
              l_res8 = log1p(res_800))
}
st <- transform_x(st)

dest_existing <- st |>
  filter(status == "open") |>
  left_join(rid, by = "station_id") |>
  mutate(lrt = as.numeric(is_lrt_only),
         ic  = pmax(lines_now - 1, 0),
         # Two kinds of tap-in. In the morning peak a station's tap-ins are people
         # leaving home and its tap-outs are people arriving at work or study, so the
         # balance between them says what share of its ridership starts from home.
         # Raffles Place is near 0 and Buangkok near 1. The rest of its ridership is
         # the return leg, tapping in near work, school or shops.
         phi = am_in / (am_in + am_out),
         T_H = daily_in * phi,
         T_W = daily_in * (1 - phi)) |>
  arrange(station_id) |>
  rename(destination = station_id)
stopifnot(!anyNA(dest_existing$daily_in))

# The second destination table describes the network once the proposed lines are in
# place. New access points have no ridership yet. They are classed as MRT: the Cross
# Island Line is heavy rail, and the Jurong Region Line, though medium capacity, runs
# as a separate line with its own codes, which a sensitivity test relaxes. An
# existing station that gains a line stays where it is, with one more line.
dest_proposed <- st |>
  filter(status == "open" | role == "new access point") |>
  left_join(rid, by = "station_id") |>
  mutate(lrt = as.numeric(is_lrt_only & status == "open"),
         ic  = ifelse(status == "open", pmax(lines_proposed - 1, 0), 0)) |>
  arrange(station_id) |>
  rename(destination = station_id)

# ---- Rail time to the CBD ------------------------------------------------------
# A station's pull depends on where it leads as well as on what surrounds it, and a
# station that needs a transfer to reach the city is worth less than one on a through
# line. Rail time to Raffles Place, with a penalty for each change, carries that. The
# proposed lines shorten it for some existing stations as well, so it is recomputed.
cx <- readRDS(file.path(DER, "station_codes.rds")) |>
  left_join(select(st, station_id, x, y), by = "station_id") |>
  select(code = station_code, station_id, station_name, x, y)
g_old <- build_graph(NETWORK, cx)
g_new <- build_graph(c(NETWORK, PROPOSED), cx)
RM_old <- rail_matrix(g_old, dest_existing$destination)
RM_new <- rail_matrix(g_new, dest_proposed$destination)
stopifnot(all(is.finite(RM_old)), all(is.finite(RM_new)))
cbd <- dest_existing$destination[dest_existing$station_name == "Raffles Place"]
dest_existing <- mutate(dest_existing, min_cbd = pmax(RM_old[, cbd], 1), l_cbd = log(min_cbd))
dest_proposed <- mutate(dest_proposed, min_cbd = pmax(RM_new[, cbd], 1), l_cbd = log(min_cbd))

X_VARS <- c("l_bus", "l_ret", "l_food", "lrt", "ic", "l_cbd")
centre <- colMeans(dest_existing[, X_VARS])

# ---- Interaction frames -------------------------------------------------------
pairs_for <- function(dst) {
  d <- sqrt(outer(origins$x, dst$x, "-")^2 + outer(origins$y, dst$y, "-")^2) / 1000
  d[d < FLOOR_KM] <- FLOOR_KM
  dimnames(d) <- list(origins$origin, dst$destination)
  d
}
D_existing <- pairs_for(dest_existing)
D_proposed <- pairs_for(dest_proposed)

as_frame <- function(D) {
  tibble(origin = rep(rownames(D), times = ncol(D)),
         destination = rep(colnames(D), each = nrow(D)),
         distance = as.vector(D))
}
si_existing <- as_frame(D_existing)
si_proposed <- as_frame(D_proposed)

audit <- list(
  n_origins = nrow(origins),
  n_origins_residents = sum(origins$residents > 0),
  n_existing = nrow(dest_existing),
  n_new = sum(dest_proposed$role == "new access point"),
  n_gains = sum(dest_proposed$role == "existing, gains a line"),
  n_pairs_existing = nrow(si_existing),
  n_pairs_proposed = nrow(si_proposed),
  res_by_nearest = with(list(nearest = apply(D_existing, 1, min)),
    round(tapply(origins$residents, cut(nearest, c(0, .4, .8, 1.2, 2, 3, 30)), sum) /
            sum(origins$residents), 4)),
  res_by_nearest_new = with(list(nearest = apply(D_proposed, 1, min)),
    round(tapply(origins$residents, cut(nearest, c(0, .4, .8, 1.2, 2, 3, 30)), sum) /
            sum(origins$residents), 4)),
  X_VARS = X_VARS,
  centre = centre
)

saveRDS(origins, file.path(DER, "origins.rds"))
saveRDS(dest_existing, file.path(DER, "dest_existing.rds"))
saveRDS(dest_proposed, file.path(DER, "dest_proposed.rds"))
saveRDS(si_existing, file.path(DER, "si_existing.rds"))
saveRDS(si_proposed, file.path(DER, "si_proposed.rds"))
saveRDS(D_existing, file.path(DER, "D_existing.rds"))
saveRDS(D_proposed, file.path(DER, "D_proposed.rds"))
saveRDS(RM_old, file.path(DER, "RM_existing.rds"))
saveRDS(RM_new, file.path(DER, "RM_proposed.rds"))
saveRDS(audit, file.path(DER, "audit_frames.rds"))
str(audit, max.level = 1, give.attr = FALSE)
print(audit$res_by_nearest); print(audit$res_by_nearest_new)
