# =============================================================================
# 01_markets.R -- origin markets: where the potential MRT users live
# Take-home Exercise 3 | ISSS626
#
# Builds the local market potential that the spatial interaction model allocates
# to stations. HDB's building layer has footprints but no dwelling counts, and
# its property table has dwelling counts but no geometry and no shared street
# field, so the two are first linked through a crosswalk recovered from the data.
# Subzone residents are then spread onto HDB blocks by flat type and onto private
# housing through OpenStreetMap buildings, and summed to a hexagon grid.
# Run from take-home-ex/take-home-ex3/.
# =============================================================================

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(tidyr); library(readr); library(stringr)
})
sf_use_s2(FALSE)   # planar operations: everything is moved to SVY21 first

RAW <- "data/raw"; DER <- "data/derived"
dir.create(DER, recursive = TRUE, showWarnings = FALSE)
SVY21 <- 3414
HEX_M <- 500       # flat-to-flat width of the origin hexagons, in metres
audit <- list()

# ---- Subzones ---------------------------------------------------------------
# The Master Plan 2019 layer has six invalid polygons and a Z dimension. Both
# would break the spatial joins below, so they are repaired on import.
sz <- st_read(file.path(RAW, "boundary/MasterPlan2019SubzoneBoundaryNoSea.geojson"),
              quiet = TRUE) |>
  st_zm()
audit$sz_invalid <- sum(!st_is_valid(sz))
sz <- sz |>
  st_make_valid() |>
  st_transform(SVY21) |>
  transmute(SZ = SUBZONE_N, PA = PLN_AREA_N, REGION = REGION_N)
stopifnot(!anyDuplicated(sz$SZ))
land <- st_union(sz)

# ---- HDB: the street crosswalk ----------------------------------------------
# The building layer identifies a street by a code (YIR00Z) and the property
# table by a name (YISHUN RING RD). No published table links them, but both list
# the same block numbers, so a street code and a street name that share the same
# set of blocks are the same street. Block numbers repeat across the island, so
# overlap alone can tie (KG KAYU RD, SAGO LANE); the code's letters, which
# abbreviate the street, break those ties.
bld <- st_read(file.path(RAW, "hdb/HDBExistingBuilding.geojson"), quiet = TRUE) |>
  st_zm()
audit$bld_n       <- nrow(bld)
audit$bld_invalid <- sum(!st_is_valid(bld))
bld <- bld |>
  st_make_valid() |>
  st_transform(SVY21) |>
  transmute(blk = toupper(trimws(BLK_NO)), ST_COD, POSTAL_COD)

prop <- read_csv(file.path(RAW, "hdb/HDBPropertyInformation.csv"),
                 show_col_types = FALSE) |>
  mutate(blk = toupper(trimws(blk_no)))
audit$prop_n <- nrow(prop)
audit$prop_du <- sum(prop$total_dwelling_units)

n_cod <- count(st_drop_geometry(bld), ST_COD, name = "n_cod")
n_str <- count(prop, street, name = "n_str")

prefix_score <- function(code, street) {
  pre   <- substr(code, 1, 3)
  words <- strsplit(gsub("[^A-Z0-9 ]", "", street), " ")[[1]]
  best  <- 0
  for (w in words) for (k in 3:1) {
    if (substr(w, 1, k) == substr(pre, 1, k)) { best <- max(best, k); break }
  }
  best / 3
}

pairs <- st_drop_geometry(bld) |>
  select(blk, ST_COD) |>
  inner_join(select(prop, blk, street), by = "blk", relationship = "many-to-many") |>
  count(ST_COD, street, name = "n_common") |>
  left_join(n_cod, by = "ST_COD") |>
  left_join(n_str, by = "street") |>
  mutate(jaccard = n_common / (n_cod + n_str - n_common),
         letters = mapply(prefix_score, ST_COD, street),
         score   = jaccard + 0.5 * letters) |>
  arrange(desc(score), desc(n_common))

# Greedy one-to-one assignment, best pair first.
take <- logical(nrow(pairs)); seen_c <- character(); seen_s <- character()
for (i in seq_len(nrow(pairs))) {
  if (!pairs$ST_COD[i] %in% seen_c && !pairs$street[i] %in% seen_s) {
    take[i] <- TRUE
    seen_c  <- c(seen_c, pairs$ST_COD[i]); seen_s <- c(seen_s, pairs$street[i])
  }
}
crosswalk <- pairs[take, ]
audit$cw_codes    <- nrow(n_cod)
audit$cw_streets  <- nrow(n_str)
audit$cw_assigned <- nrow(crosswalk)
audit$cw_perfect  <- sum(crosswalk$jaccard == 1)
audit$cw_unassigned_codes   <- setdiff(n_cod$ST_COD, crosswalk$ST_COD)
audit$cw_unassigned_streets <- setdiff(n_str$street, crosswalk$street)

bld_pts <- bld |>
  st_point_on_surface() |>
  left_join(select(crosswalk, ST_COD, street), by = "ST_COD")

blocks <- prop |>
  filter(residential == "Y") |>
  left_join(select(bld_pts, blk, street, geometry), by = c("blk", "street")) |>
  mutate(has_footprint = !st_is_empty(st_sfc(geometry, crs = SVY21)))

audit$bld_matched     <- sum(paste(bld_pts$blk, bld_pts$street) %in%
                               paste(prop$blk, prop$street))
audit$res_blocks      <- nrow(blocks)
audit$res_du          <- sum(blocks$total_dwelling_units)
audit$res_no_footprint    <- sum(!blocks$has_footprint)
audit$res_no_footprint_du <- sum(blocks$total_dwelling_units[!blocks$has_footprint])
audit$res_no_footprint_years <- sort(unique(blocks$year_completed[!blocks$has_footprint]))

# Blocks completed after the footprint layer was last revised have no geometry.
# Each is placed at the centre of the other blocks on its street, or left for the
# geocoded list below when the whole street is new.
street_centre <- blocks |>
  filter(has_footprint) |>
  st_as_sf(crs = SVY21) |>
  mutate(x = st_coordinates(geometry)[, 1], y = st_coordinates(geometry)[, 2]) |>
  st_drop_geometry() |>
  group_by(street) |>
  summarise(sx = mean(x), sy = mean(y), .groups = "drop")

xy <- st_coordinates(st_as_sf(filter(blocks, has_footprint), crs = SVY21))
blocks$x <- NA_real_; blocks$y <- NA_real_
blocks$x[blocks$has_footprint] <- xy[, 1]
blocks$y[blocks$has_footprint] <- xy[, 2]
blocks <- blocks |>
  select(-geometry) |>
  left_join(street_centre, by = "street") |>
  mutate(placed = case_when(has_footprint ~ "footprint",
                            !is.na(sx) ~ "street centre",
                            TRUE ~ "none"),
         x = coalesce(x, sx), y = coalesce(y, sy)) |>
  select(-sx, -sy)
audit$res_placed <- count(blocks, placed, wt = total_dwelling_units)
unplaced <- filter(blocks, placed == "none")
audit$res_unplaced <- unplaced |> select(blk, street, year_completed, total_dwelling_units)
blocks <- filter(blocks, placed != "none")

# Flat types, grouped to match the dwelling categories of the population table.
blocks <- blocks |>
  mutate(g12 = `1room_sold` + `2room_sold` + studio_apartment_sold + `1room_rental` +
               `2room_rental` + other_room_rental,
         g3  = `3room_sold` + `3room_rental`,
         g4  = `4room_sold`,
         g5  = `5room_sold` + exec_sold + multigen_sold) |>
  st_as_sf(coords = c("x", "y"), crs = SVY21) |>
  st_join(select(sz, SZ, PA), join = st_nearest_feature)

# ---- Residents by subzone ---------------------------------------------------
pop <- read_csv(file.path(RAW, "popn/respopagesextod2024.csv"), show_col_types = FALSE) |>
  mutate(SZ = toupper(SZ)) |>
  group_by(SZ, TOD) |>
  summarise(Pop = sum(Pop), .groups = "drop") |>
  mutate(group = case_when(
    TOD == "HDB 1- and 2-Room Flats" ~ "g12",
    TOD == "HDB 3-Room Flats" ~ "g3",
    TOD == "HDB 4-Room Flats" ~ "g4",
    TOD %in% c("HDB 5-Room and Executive Flats",
               "HUDC Flats (excluding those privatised)") ~ "g5",
    TOD == "Condominiums and Other Apartments" ~ "condo",
    TOD == "Landed Properties" ~ "landed",
    TRUE ~ "other")) |>
  group_by(SZ, group) |>
  summarise(Pop = sum(Pop), .groups = "drop")

audit$pop_total <- sum(pop$Pop)
audit$pop_by_group <- pop |> group_by(group) |> summarise(Pop = sum(Pop), .groups = "drop")
audit$pop_sz_unmatched <- setdiff(unique(pop$SZ[pop$Pop > 0]), sz$SZ)

# ---- Residents onto HDB blocks ----------------------------------------------
# Persons per flat are computed within each subzone and flat type, from blocks
# standing at the population reference date (June 2024). Blocks completed since
# then are not in that count, so they receive the island-wide persons per flat of
# their type instead: the ridership being modelled is from 2026, when they are
# occupied. Blocks completed during 2026 are left out of the base market because
# their keys are still being collected, and kept for the 2030 scenario.
REF_YEAR <- 2024
GROUPS <- c("g12", "g3", "g4", "g5")

long <- blocks |>
  st_drop_geometry() |>
  select(blk, street, SZ, year_completed, all_of(GROUPS)) |>
  pivot_longer(all_of(GROUPS), names_to = "group", values_to = "flats") |>
  filter(flats > 0)

stand <- long |>
  filter(year_completed <= REF_YEAR) |>
  group_by(SZ, group) |>
  summarise(flats_sz = sum(flats), .groups = "drop")

ppf_sz <- pop |>
  filter(group %in% GROUPS) |>
  full_join(stand, by = c("SZ", "group")) |>
  mutate(Pop = coalesce(Pop, 0), flats_sz = coalesce(flats_sz, 0),
         ppf = ifelse(flats_sz > 0, Pop / flats_sz, NA_real_))

ppf_city <- ppf_sz |>
  group_by(group) |>
  summarise(ppf_city = sum(Pop[flats_sz > 0]) / sum(flats_sz), .groups = "drop")
audit$ppf_city <- ppf_city

# Residents recorded in a flat type a subzone does not have: the two sources
# disagree slightly because HUDC estates and en-bloc sites are classified
# differently. These are spread over the subzone's other blocks.
orphans <- ppf_sz |> filter(flats_sz == 0, Pop > 0)
audit$pop_orphan_hdb <- sum(orphans$Pop)

long <- long |>
  left_join(select(ppf_sz, SZ, group, ppf), by = c("SZ", "group")) |>
  left_join(ppf_city, by = "group") |>
  mutate(is_new = year_completed > REF_YEAR,
         res = flats * ifelse(is_new | is.na(ppf), ppf_city, ppf))

orphan_sz <- orphans |> group_by(SZ) |> summarise(orphan = sum(Pop), .groups = "drop")
block_res <- long |>
  group_by(blk, street, SZ) |>
  summarise(res = sum(res), .groups = "drop") |>
  left_join(orphan_sz, by = "SZ") |>
  group_by(SZ) |>
  mutate(res = res + coalesce(orphan, 0) * res / sum(res)) |>
  ungroup() |>
  select(-orphan)

blocks <- blocks |>
  left_join(select(block_res, blk, street, res), by = c("blk", "street")) |>
  mutate(res = coalesce(res, 0))
audit$hdb_res_base <- sum(blocks$res[blocks$year_completed <= 2025])
audit$hdb_res_2026 <- sum(blocks$res[blocks$year_completed == 2026])
audit$hdb_du_2026  <- sum(blocks$total_dwelling_units[blocks$year_completed == 2026])
audit$hdb_du_new_since_ref <- sum(blocks$total_dwelling_units[blocks$year_completed > REF_YEAR &
                                                                blocks$year_completed <= 2025])

# ---- Residents in private housing -------------------------------------------
# Nearly a quarter of residents live outside HDB estates, and no official layer
# locates them. Condominium residents are placed on OpenStreetMap apartment and residential
# buildings by footprint area, and landed residents on house-type buildings by
# count, within each subzone. Buildings on an HDB footprint are removed first,
# because OSM tags HDB blocks as residential too.
Z <- paste0("/vsizip/", RAW, "/osm/malaysia-singapore-brunei-latest-free.shp.zip/")
bb <- st_as_text(st_as_sfc(st_bbox(st_transform(sz, 4326))))
osm_b <- st_read(paste0(Z, "gis_osm_buildings_a_free_1.shp"), wkt_filter = bb, quiet = TRUE) |>
  st_make_valid() |>
  st_transform(SVY21) |>
  select(type)
osm_b$area <- as.numeric(st_area(osm_b))
osm_pts <- st_point_on_surface(osm_b)
osm_pts <- st_filter(osm_pts, land)
audit$osm_b_n <- nrow(osm_pts)

on_hdb <- lengths(st_intersects(osm_pts, bld)) > 0
audit$osm_b_on_hdb <- sum(on_hdb)
priv <- osm_pts[!on_hdb, ] |> st_join(select(sz, SZ), join = st_intersects)

# Dormitories are left out: their occupants are mostly foreign workers and students
# who are not in the resident count. They are carried as a separate market in
# 02_stations.R.
condo_types  <- c("apartments", "residential")
landed_types <- c("house", "terrace", "semidetached_house", "detached", "bungalow")
cand <- bind_rows(
  priv |> filter(type %in% condo_types)  |> mutate(group = "condo",  w = area),
  priv |> filter(type %in% landed_types) |> mutate(group = "landed", w = 1)
)

# "Others" (shophouses, attap and so on) has no building type of its own, and is
# treated with condominiums. Where a subzone has residents but no candidate of the
# right type, the other private type is tried, then OSM residential land, and
# finally the subzone's interior, each of which is counted in the audit.
pop_priv <- pop |>
  filter(group %in% c("condo", "landed", "other")) |>
  mutate(group = ifelse(group == "other", "condo", group)) |>
  group_by(SZ, group) |>
  summarise(Pop = sum(Pop), .groups = "drop") |>
  filter(Pop > 0)

landuse <- st_read(paste0(Z, "gis_osm_landuse_a_free_1.shp"), wkt_filter = bb, quiet = TRUE) |>
  filter(fclass == "residential") |>
  st_make_valid() |>
  st_transform(SVY21)
landuse$w <- as.numeric(st_area(landuse))
lu_pts <- st_point_on_surface(landuse) |> st_join(select(sz, SZ), join = st_intersects)

allocate_private <- function(s, g, p) {
  this <- cand |> filter(SZ == s, group == g)
  tier <- "same type"
  if (nrow(this) == 0) { this <- cand |> filter(SZ == s); tier <- "other private type" }
  if (nrow(this) == 0) { this <- lu_pts |> filter(SZ == s); tier <- "OSM residential land" }
  if (nrow(this) == 0) {
    this <- st_sf(w = 1, geometry = st_point_on_surface(st_geometry(sz[sz$SZ == s, ])))
    tier <- "subzone interior"
  }
  xy <- st_coordinates(this)
  tibble(x = xy[, 1], y = xy[, 2], res = p * this$w / sum(this$w), tier = tier, group = g, SZ = s)
}
priv_res <- bind_rows(Map(allocate_private, pop_priv$SZ, pop_priv$group, pop_priv$Pop))
audit$priv_total <- sum(priv_res$res)
audit$priv_by_tier <- priv_res |> group_by(tier) |> summarise(res = sum(res), .groups = "drop")

# ---- Hexagon grid -----------------------------------------------------------
# Hexagons are equidistant from their six neighbours and small enough for a walk
# to a station to differ across them. Each hexagon's distance point is its
# population-weighted centre rather than its geometric centre, which removes most
# of the error of treating a cell as a point.
hex <- st_make_grid(land, cellsize = HEX_M, square = FALSE) |>
  st_sf(geometry = _)
hex <- hex[lengths(st_intersects(hex, land)) > 0, ]
hex$hex_id <- sprintf("H%04d", seq_len(nrow(hex)))
audit$hex_n <- nrow(hex)

hdb_pts <- blocks |>
  filter(year_completed <= 2025) |>
  mutate(x = st_coordinates(geometry)[, 1], y = st_coordinates(geometry)[, 2]) |>
  st_join(select(hex, hex_id), join = st_intersects) |>
  st_drop_geometry() |>
  select(hex_id, x, y, res, du = total_dwelling_units)

priv_pts <- priv_res |>
  st_as_sf(coords = c("x", "y"), crs = SVY21, remove = FALSE) |>
  st_join(select(hex, hex_id), join = st_nearest_feature) |>
  st_drop_geometry() |>
  select(hex_id, x, y, res)

people <- bind_rows(mutate(hdb_pts, src = "hdb"), mutate(priv_pts, src = "private", du = 0))

hex_tot <- people |>
  group_by(hex_id) |>
  summarise(wx = sum(x * res) / sum(res), wy = sum(y * res) / sum(res),
            residents = sum(res),
            res_hdb = sum(res[src == "hdb"]), res_private = sum(res[src == "private"]),
            hdb_du = sum(du), .groups = "drop")

hex <- hex |>
  left_join(hex_tot, by = "hex_id") |>
  mutate(across(c(residents, res_hdb, res_private, hdb_du), \(v) coalesce(v, 0)))
hex_c <- st_coordinates(st_centroid(st_geometry(hex)))
hex <- hex |>
  mutate(cx = hex_c[, 1], cy = hex_c[, 2],
         x = coalesce(wx, cx), y = coalesce(wy, cy)) |>
  select(-wx, -wy, -cx, -cy)
# A hexagon can straddle subzones, so it is labelled by its population-weighted
# centre rather than by whichever subzone its polygon meets first.
hex_lab <- st_as_sf(st_drop_geometry(hex), coords = c("x", "y"), crs = SVY21) |>
  st_join(select(sz, SZ, PA, REGION), join = st_nearest_feature) |>
  st_drop_geometry() |>
  select(hex_id, SZ, PA, REGION)
hex <- left_join(hex, hex_lab, by = "hex_id")

audit$res_in_hex <- sum(hex$residents)
audit$hex_populated <- sum(hex$residents > 0)
audit$share_private <- sum(hex$res_private) / sum(hex$residents)

# Blocks completed in 2026, kept for the 2030 scenario.
new_2026 <- blocks |>
  filter(year_completed == 2026) |>
  select(blk, street, res, du = total_dwelling_units)

saveRDS(sz, file.path(DER, "subzone.rds"))
saveRDS(st_union(sz), file.path(DER, "land.rds"))
saveRDS(blocks, file.path(DER, "hdb_blocks.rds"))
saveRDS(crosswalk, file.path(DER, "street_crosswalk.rds"))
saveRDS(new_2026, file.path(DER, "hdb_new_2026.rds"))
saveRDS(hex, file.path(DER, "hex.rds"))
saveRDS(select(people, x, y, res, src, du), file.path(DER, "people_points.rds"))
saveRDS(audit, file.path(DER, "audit_markets.rds"))

str(audit, max.level = 1, give.attr = FALSE)
