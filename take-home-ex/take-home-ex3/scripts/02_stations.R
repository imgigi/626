# =============================================================================
# 02_stations.R -- destinations: existing and proposed stations, and their pull
# Take-home Exercise 3 | ISSS626
#
# A station is a place, not a code. An interchange has several codes, and a future
# line can bring a new code to a station that already exists, so codes are first
# collapsed into physical stations. The attraction variables are then measured
# around each station, on the same definitions for existing and proposed ones,
# because a proposed station has nothing else to be predicted from.
# Run from take-home-ex/take-home-ex3/.
# =============================================================================

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(tidyr); library(readr); library(stringr); library(purrr)
})
sf_use_s2(FALSE)

RAW <- "data/raw"; DER <- "data/derived"; SVY21 <- 3414
audit <- list()
land <- readRDS(file.path(DER, "land.rds"))
hex  <- readRDS(file.path(DER, "hex.rds"))

# ---- LTA station layer ------------------------------------------------------
# Station identifiers sit in STN_NAM as codes, and the layer's SVY21 definition has
# no EPSG code. It also holds depots and ancillary buildings with no passengers
# (their codes are not station codes), several polygons for some single stations,
# and shared polygons named "EW24 / NS1".
ly <- st_read(paste0("/vsizip/", RAW, "/lta/TrainStation_Sep2026.zip/TrainStation_Sep2026/",
                     "RapidTransitSystemStation.shp"), quiet = TRUE) |>
  st_set_crs(NA) |> st_set_crs(SVY21) |>
  st_make_valid()
audit$layer_n <- nrow(ly)
audit$layer_dup <- st_drop_geometry(ly) |> count(STN_NAM) |> filter(n > 1)
audit$layer_shared <- ly$STN_NAM[grepl("/", ly$STN_NAM)]

is_station_code <- function(x) str_detect(x, "^(NS|EW|CG|NE|CC|CE|DT|TE|BP|SE|SW|PE|PW|STC|PTC)[0-9]*[A-Z]?$")
ly_codes <- st_drop_geometry(ly) |>
  mutate(code = str_split(STN_NAM, "\\s*/\\s*")) |>
  unnest(code) |>
  mutate(code = str_trim(code))
audit$layer_nonstation <- sort(unique(ly_codes$code[!is_station_code(ly_codes$code)]))
layer_codes <- unique(ly_codes$code[is_station_code(ly_codes$code)])

# ---- Station lists ----------------------------------------------------------
# The coordinate lists carry a code, a name and a point for every station, built
# from URA's Master Plan. They date from March 2026, so three stations they call
# future (Circle Line Stage 6) had opened by the September layer and the OD months.
op <- read_csv(file.path(RAW, "stations/operational_stations.csv"), show_col_types = FALSE)
fu <- read_csv(file.path(RAW, "stations/future_stations.csv"), show_col_types = FALSE)
all_codes <- bind_rows(mutate(op, listed = "operational"), mutate(fu, listed = "future")) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) |>
  st_transform(SVY21) |>
  mutate(x = st_coordinates(geometry)[, 1], y = st_coordinates(geometry)[, 2],
         line = str_extract(station_code, "^[A-Z]+")) |>
  st_drop_geometry()

audit$codes_operational <- nrow(op)
audit$codes_future <- nrow(fu)
audit$op_not_in_layer <- setdiff(op$station_code, layer_codes)
audit$future_in_layer <- intersect(fu$station_code, layer_codes)

# What each code is, as of the September 2026 layer and LTA's project pages.
# The handout's proposed stations are the Jurong Region Line stations opening by
# 2029 (24, leaving out the mid-2030s infill at JS2A) and Cross Island Line Phase 1
# (12 stations, CR2 to CR13, 2030).
CCL6    <- c("CC30", "CC31", "CC32")
SHELLS  <- c("CC18", "TE10", "TE21", "TE22A")
JRL     <- c(paste0("JS", 1:12), paste0("JE", 1:7), paste0("JW", 1:5))
CRL1    <- paste0("CR", 2:13)
all_codes <- all_codes |>
  mutate(status = case_when(
    listed == "operational" | station_code %in% CCL6 ~ "open",
    station_code %in% SHELLS ~ "shell",
    station_code %in% c(JRL, CRL1) ~ "proposed",
    TRUE ~ "later"))

# ---- Physical stations ------------------------------------------------------
# Codes sharing a name are one station. The check below confirms that the name
# never joins two places: the largest gap inside a name is reported in the audit.
gap <- all_codes |>
  group_by(station_name) |>
  summarise(n = n(), gap_m = max(as.numeric(dist(cbind(x, y))), 0), .groups = "drop")
audit$name_gap <- gap |> filter(n > 1) |> arrange(desc(gap_m)) |> head(12)

# An existing station's position stays as it is when a future line arrives: a
# passenger reaches the same entrances, so the new code adds a line and not a place.
stations <- all_codes |>
  mutate(rank_status = match(status, c("open", "shell", "proposed", "later"))) |>
  group_by(station_name) |>
  summarise(
    codes_all = paste(station_code, collapse = "/"),
    lines_now      = n_distinct(line[status == "open"]),
    lines_proposed = n_distinct(line[status %in% c("open", "proposed")]),
    is_lrt_only    = any(status == "open") &&
      all(line[status == "open"] %in% c("BP", "SE", "SW", "PE", "PW", "STC", "PTC")),
    has_new_code   = any(status == "proposed"),
    x = mean(x[rank_status == min(rank_status)]),
    y = mean(y[rank_status == min(rank_status)]),
    # Assigned last: from here on, `status` would mean the station's, not the code's.
    status    = status[which.min(rank_status)],
    .groups = "drop") |>
  mutate(station_id = sprintf("S%03d", row_number()),
         # LRT stations and MRT lines count as one mode change each
         is_interchange_now = lines_now >= 2,
         role = case_when(
           status == "open" & has_new_code ~ "existing, gains a line",
           status == "open" ~ "existing",
           status == "proposed" ~ "new access point",
           status == "shell" ~ "built, not open",
           TRUE ~ "after 2030"))
audit$role_counts <- count(stations, role)
audit$proposed_codes_n <- sum(all_codes$station_code %in% c(JRL, CRL1))

code_map <- all_codes |>
  select(station_code, station_name, line, status) |>
  left_join(select(stations, station_name, station_id), by = "station_name")

st_pts <- stations |>
  st_as_sf(coords = c("x", "y"), crs = SVY21, remove = FALSE)

# ---- Attraction variables ---------------------------------------------------
Z  <- paste0("/vsizip/", RAW, "/osm/malaysia-singapore-brunei-latest-free.shp.zip/")
bb <- st_as_text(st_as_sfc(st_bbox(st_transform(st_as_sf(land), 4326))))
read_osm <- function(l) {
  st_read(paste0(Z, l, ".shp"), wkt_filter = bb, quiet = TRUE) |>
    st_make_valid() |>
    st_transform(SVY21)
}

# Bus stops. The file has 21 codes entered twice and five stops across the Causeway
# in Johor Bahru, which are removed; a stop's own code is what makes it one stop.
bus <- st_read(paste0("/vsizip/", RAW, "/lta/BusStopLocation_Mar2026.zip/",
                      "BusStopLocation_Mar2026/BusStop.shp"), quiet = TRUE) |>
  st_set_crs(NA) |> st_set_crs(SVY21)
audit$bus_n <- nrow(bus)
audit$bus_dup <- sum(duplicated(bus$BUS_STOP_N))
out <- lengths(st_intersects(bus, st_buffer(land, 100))) == 0
audit$bus_outside <- st_drop_geometry(bus)[out, c("BUS_STOP_N", "LOC_DESC")]
bus <- bus[!out, ] |> distinct(BUS_STOP_N, .keep_all = TRUE)
audit$bus_clean <- nrow(bus)

# Points of interest, grouped into the kinds of trip they generate. Street furniture
# that OSM also records (benches, waste baskets, vending machines) is left out.
POI_GROUPS <- list(
  food    = c("restaurant", "cafe", "fast_food", "food_court", "bar", "pub", "biergarten"),
  retail  = c("mall", "supermarket", "convenience", "department_store", "clothes", "shoe_shop",
              "jeweller", "optician", "beauty_shop", "furniture_shop", "bakery", "doityourself",
              "general", "gift_shop", "bookshop", "computer_shop", "mobile_phone_shop",
              "sports_shop", "toy_shop", "stationery", "florist", "marketplace", "beverages",
              "butcher", "greengrocer", "kiosk", "chemist", "garden_centre", "outdoor_shop",
              "bicycle_shop", "car_dealership", "hairdresser", "laundry", "bank", "travel_agent"),
  civic   = c("hospital", "clinic", "doctors", "dentist", "pharmacy", "nursing_home", "veterinary",
              "school", "college", "university", "kindergarten", "library", "community_centre",
              "courthouse", "town_hall", "post_office", "police", "fire_station", "embassy"),
  leisure = c("attraction", "museum", "theatre", "cinema", "arts_centre", "zoo", "theme_park",
              "stadium", "sports_centre", "fitness_centre", "hotel", "hostel", "guesthouse",
              "motel", "nightclub", "golf_course", "swimming_pool"))
poi_lookup <- tibble(group = rep(names(POI_GROUPS), lengths(POI_GROUPS)),
                     fclass = unlist(POI_GROUPS))
poi <- bind_rows(
  read_osm("gis_osm_pois_free_1") |> select(osm_id, fclass),
  read_osm("gis_osm_pois_a_free_1") |> select(osm_id, fclass) |> st_point_on_surface()
) |>
  st_filter(land) |>
  inner_join(poi_lookup, by = "fclass")
audit$poi_n <- nrow(poi)
audit$poi_by_group <- st_drop_geometry(poi) |> count(group)

# Floor space for work and retail. OSM has no heights, so this is footprint, which
# understates the tall offices of the CBD; that bias is examined in the report.
WORK_TYPES <- c("commercial", "office", "retail", "industrial", "warehouse", "hotel", "public",
                "government", "hospital", "school", "university", "college", "supermarket")
osm_b <- read_osm("gis_osm_buildings_a_free_1") |> select(type)
osm_b$area <- as.numeric(st_area(osm_b))
work <- osm_b |>
  filter(type %in% WORK_TYPES) |>
  st_point_on_surface() |>
  st_filter(land)
audit$work_n <- nrow(work)

# Dormitories house people the resident count leaves out: foreign workers above all,
# and students. OSM is the only available record of where they are.
dorm <- osm_b |>
  filter(type == "dormitory") |>
  st_point_on_surface() |>
  st_filter(land)
audit$dorm_n <- nrow(dorm)

lu <- read_osm("gis_osm_landuse_a_free_1") |>
  filter(fclass %in% c("commercial", "retail", "industrial")) |>
  st_intersection(st_as_sf(land)) |>
  select(fclass)

count_within <- function(pts, centres, r) lengths(st_intersects(st_buffer(centres, r), pts))
sum_within <- function(pts, val, centres, r) {
  hits <- st_intersects(st_buffer(centres, r), pts)
  vapply(hits, function(i) sum(val[i]), numeric(1))
}
lu_area_within <- function(centres, r) {
  buf <- st_buffer(centres, r) |> select(station_id)
  st_intersection(buf, st_geometry(st_as_sf(st_union(lu)))) |>
    mutate(a = as.numeric(st_area(geometry))) |>
    st_drop_geometry() |>
    group_by(station_id) |>
    summarise(a = sum(a), .groups = "drop") |>
    right_join(st_drop_geometry(buf), by = "station_id") |>
    mutate(a = coalesce(a, 0)) |>
    arrange(match(station_id, buf$station_id)) |>
    pull(a)
}

hex_pts <- hex |>
  st_drop_geometry() |>
  filter(residents > 0) |>
  st_as_sf(coords = c("x", "y"), crs = SVY21)

RAFFLES <- all_codes |> filter(station_name == "Raffles Place") |> summarise(x = mean(x), y = mean(y))

attr_tbl <- st_pts |>
  transmute(
    station_id,
    bus_400     = count_within(bus, geometry, 400),
    food_400    = count_within(filter(poi, group == "food"), geometry, 400),
    retail_400  = count_within(filter(poi, group == "retail"), geometry, 400),
    civic_400   = count_within(filter(poi, group == "civic"), geometry, 400),
    leisure_400 = count_within(filter(poi, group == "leisure"), geometry, 400),
    poi_800     = count_within(poi, geometry, 800),
    work_m2_400 = sum_within(work, work$area, geometry, 400),
    work_m2_800 = sum_within(work, work$area, geometry, 800),
    res_400     = sum_within(hex_pts, hex_pts$residents, geometry, 400),
    res_800     = sum_within(hex_pts, hex_pts$residents, geometry, 800),
    km_cbd      = sqrt((x - RAFFLES$x)^2 + (y - RAFFLES$y)^2) / 1000
  ) |>
  st_drop_geometry()
attr_tbl$lu_m2_400 <- lu_area_within(select(st_pts, station_id), 400)

stations_full <- stations |>
  left_join(attr_tbl, by = "station_id")

# ---- Origin-side activity ---------------------------------------------------
# The second market. Tap-ins away from home come from people returning from work,
# school and shops, so each hexagon also carries the activity located in it.
hex_poly <- select(hex, hex_id)
add_hex <- function(pts, val = NULL) {
  j <- st_join(pts, hex_poly, join = st_intersects) |> st_drop_geometry()
  j$v <- if (is.null(val)) 1 else val
  j |> filter(!is.na(hex_id)) |> group_by(hex_id) |> summarise(v = sum(v), .groups = "drop")
}
activity <- select(st_drop_geometry(hex), hex_id) |>
  left_join(rename(add_hex(poi), poi = v), by = "hex_id") |>
  left_join(rename(add_hex(work, work$area), work_m2 = v), by = "hex_id") |>
  left_join(rename(add_hex(dorm, dorm$area), dorm_m2 = v), by = "hex_id") |>
  mutate(across(c(poi, work_m2, dorm_m2), \(v) coalesce(v, 0)))
hex2 <- left_join(hex, activity, by = "hex_id")
audit$hex_with_market <- sum(hex2$residents > 0 | hex2$poi > 0 | hex2$work_m2 > 0 | hex2$dorm_m2 > 0)

saveRDS(stations_full, file.path(DER, "stations.rds"))
saveRDS(code_map, file.path(DER, "station_codes.rds"))
saveRDS(hex2, file.path(DER, "hex_markets.rds"))
activity_pts <- bind_rows(
  mutate(st_drop_geometry(work), x = st_coordinates(work)[, 1], y = st_coordinates(work)[, 2], kind = "work"),
  mutate(st_drop_geometry(dorm), x = st_coordinates(dorm)[, 1], y = st_coordinates(dorm)[, 2], kind = "dorm"),
  mutate(st_drop_geometry(poi), x = st_coordinates(poi)[, 1], y = st_coordinates(poi)[, 2], area = 1, kind = "poi")
) |> select(x, y, area, kind)
saveRDS(activity_pts, file.path(DER, "activity_points.rds"))
saveRDS(audit, file.path(DER, "audit_stations.rds"))

str(audit, max.level = 1, give.attr = FALSE)
print(audit$name_gap); print(audit$role_counts); print(audit$poi_by_group)
print(as.data.frame(stations_full |> filter(role %in% c("existing, gains a line", "built, not open")) |>
                      select(station_name, codes_all, role, lines_now, lines_proposed)))
print(summary(select(stations_full, bus_400, food_400, retail_400, civic_400, leisure_400, work_m2_400, res_800, lu_m2_400)))
