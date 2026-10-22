# =============================================================================
# 02_acled.R -- import, audit and geocode the ACLED event records
# Take-home Exercise 2 | ISSS626
#
# Turns the ACLED extract into an audited event-level sf object, the township
# totals the local measures of spatial association are computed on, and the
# township-by-quarter and township-by-month panels the space-time cube is
# built from. Every record that is dropped is counted and the reason is
# recorded, so the report can state what the analysis is run on rather than
# implying it used everything.
# Run from take-home-ex/take-home-ex2/.
# =============================================================================

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(tidyr); library(stringr)
  library(lubridate); library(readr)
})

RAW <- "data/raw"; DER <- "data/derived"
dir.create(DER, recursive = TRUE, showWarnings = FALSE)

STUDY_START <- as.Date("2021-01-01")
STUDY_END   <- as.Date("2025-09-30")

# ---- Locate the extract ------------------------------------------------------
# The file is distributed through the course rather than fetched from an API,
# so its name is not fixed. Anything in data/raw/ that looks like a tabular
# ACLED export is accepted, and the chosen file is recorded in the audit.
cand <- list.files(RAW, pattern = "\\.(csv|xlsx|xls)$", full.names = TRUE,
                   ignore.case = TRUE, recursive = TRUE)
stopifnot("No ACLED extract found in data/raw/" = length(cand) > 0)
src <- cand[which.max(file.size(cand))]

read_any <- function(p) {
  if (grepl("\\.csv$", p, ignore.case = TRUE)) {
    readr::read_csv(p, show_col_types = FALSE, guess_max = 200000,
                    name_repair = "minimal")
  } else {
    readxl::read_excel(p, guess_max = 200000)
  }
}
raw <- read_any(src)

# ACLED has shipped the same fields under different capitalisation across
# releases; names are normalised once so nothing downstream depends on which
# release produced the file.
names(raw) <- tolower(gsub("[^A-Za-z0-9]+", "_", names(raw)))

need <- c("event_date", "event_type", "latitude", "longitude", "fatalities",
          "admin3", "geo_precision", "sub_event_type")
missing_cols <- setdiff(need, names(raw))
stopifnot("ACLED extract is missing required columns" = length(missing_cols) == 0)

n_raw <- nrow(raw)

# ---- Typing ------------------------------------------------------------------
# event_date arrives either as a date or as a string in one of two orders;
# parsing is explicit so that a silent NA never turns into a dropped month.
parse_event_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXct")) return(as.Date(x))
  d <- suppressWarnings(lubridate::dmy(x, quiet = TRUE))
  if (sum(is.na(d)) > 0.5 * length(d)) d <- suppressWarnings(lubridate::ymd(x, quiet = TRUE))
  d
}

ev <- raw |>
  mutate(
    event_date = parse_event_date(event_date),
    latitude   = suppressWarnings(as.numeric(latitude)),
    longitude  = suppressWarnings(as.numeric(longitude)),
    fatalities = suppressWarnings(as.numeric(fatalities)),
    across(any_of(c("geo_precision", "time_precision", "year")),
           ~ suppressWarnings(as.integer(.x)))
  )

# ---- Audit: the reasons a record cannot be used ------------------------------
# Counted before any filtering so the proportions in the report refer to the
# extract as delivered.
audit <- list(
  source_file   = basename(src),
  n_raw         = n_raw,
  n_cols        = ncol(raw),
  n_dup_id      = if ("event_id_cnty" %in% names(ev)) sum(duplicated(ev$event_id_cnty)) else NA_integer_,
  n_bad_date    = sum(is.na(ev$event_date)),
  n_out_window  = sum(!is.na(ev$event_date) &
                      (ev$event_date < STUDY_START | ev$event_date > STUDY_END)),
  n_no_coord    = sum(is.na(ev$latitude) | is.na(ev$longitude)),
  n_zero_coord  = sum(!is.na(ev$latitude) & !is.na(ev$longitude) &
                      ev$latitude == 0 & ev$longitude == 0),
  n_na_fatal    = sum(is.na(ev$fatalities)),
  geo_precision = table(ev$geo_precision, useNA = "ifany"),
  time_precision= if ("time_precision" %in% names(ev)) table(ev$time_precision, useNA = "ifany") else NULL,
  event_type    = sort(table(ev$event_type), decreasing = TRUE),
  disorder_type = if ("disorder_type" %in% names(ev)) sort(table(ev$disorder_type), decreasing = TRUE) else NULL,
  country       = if ("country" %in% names(ev)) sort(table(ev$country), decreasing = TRUE) else NULL,
  date_range    = range(ev$event_date, na.rm = TRUE),
  # ACLED's timestamp is the last time a record was modified. Its maximum dates
  # the extract, which matters because ACLED revises and back-fills records.
  extracted     = if ("timestamp" %in% names(ev))
                    as.Date(as.POSIXct(max(as.numeric(ev$timestamp), na.rm = TRUE),
                                       origin = "1970-01-01", tz = "UTC")) else NA
)

# ---- Filtering ---------------------------------------------------------------
# Duplicates are removed on ACLED's own event identifier rather than on the
# whole row: two genuinely distinct events can share a date, a location and an
# event type, and dropping one of them would understate intensity.
ev2 <- ev |>
  { \(d) if ("event_id_cnty" %in% names(d)) distinct(d, event_id_cnty, .keep_all = TRUE) else d }() |>
  filter(!is.na(event_date), event_date >= STUDY_START, event_date <= STUDY_END) |>
  filter(!is.na(latitude), !is.na(longitude), !(latitude == 0 & longitude == 0))

# The Myanmar extract also carries a handful of maritime records coded to the
# "Indian Ocean" pseudo-country. They have no township and are removed here.
if ("country" %in% names(ev2)) {
  audit$n_not_myanmar <- sum(ev2$country != "Myanmar")
  audit$not_myanmar   <- ev2 |> filter(country != "Myanmar") |>
    count(country, event_type, sub_event_type)
  ev2 <- filter(ev2, country == "Myanmar")
}

audit$n_after_filter <- nrow(ev2)

# ---- Geocoding to townships --------------------------------------------------
# The join is spatial rather than by name. ACLED carries its own admin3 string,
# but it is maintained against ACLED's gazetteer and not against MIMU v9.4, so
# matching on text would silently lose every township whose romanisation
# differs. The name field is kept and compared afterwards as a check on the
# spatial join rather than as the join key. The join is made against all 330
# MIMU townships, so that an event falling in the township later excluded from
# the analytical frame would be seen rather than silently lost.
ts    <- readRDS(file.path(DER, "township.rds"))
ts_an <- readRDS(file.path(DER, "township_an.rds"))
WORK_CRS <- st_crs(ts)

ev_sf <- st_as_sf(ev2, coords = c("longitude", "latitude"), crs = 4326,
                  remove = FALSE) |>
  st_transform(WORK_CRS)

sf_use_s2(FALSE)
idx <- st_within(ev_sf, ts, sparse = TRUE)
hit <- lengths(idx) > 0L

# Events that fall just outside a polygon -- coastal rounding, or a coordinate
# on a boundary line -- are recovered by nearest township, but only within a
# tolerance. Beyond it the record is dropped rather than attached to a
# township it may not belong to.
NEAR_TOL <- 5000  # metres

ev_sf$ts_idx <- NA_integer_
ev_sf$ts_idx[hit] <- vapply(idx[hit], function(i) i[1], integer(1))

if (any(!hit)) {
  out <- which(!hit)
  nn  <- st_nearest_feature(ev_sf[out, ], ts)
  dd  <- as.numeric(st_distance(ev_sf[out, ], ts[nn, ], by_element = TRUE))
  ok  <- dd <= NEAR_TOL
  ev_sf$ts_idx[out[ok]] <- nn[ok]
  audit$n_outside_poly  <- length(out)
  audit$n_snapped       <- sum(ok)
  audit$n_dropped_far   <- sum(!ok)
  audit$max_snap_dist_m <- if (any(ok)) max(dd[ok]) else 0
} else {
  audit$n_outside_poly <- 0L; audit$n_snapped <- 0L
  audit$n_dropped_far <- 0L; audit$max_snap_dist_m <- 0
}

ev_sf <- ev_sf |>
  filter(!is.na(ts_idx)) |>
  mutate(TS_PCODE = ts$TS_PCODE[ts_idx],
         TS       = ts$TS[ts_idx],
         DT       = ts$DT[ts_idx],
         ST       = ts$ST[ts_idx],
         ym       = floor_date(event_date, "month"),
         yq       = floor_date(event_date, "quarter"))

audit$n_geocoded    <- nrow(ev_sf)
audit$n_in_excluded <- sum(!ev_sf$TS_PCODE %in% ts_an$TS_PCODE)

# Agreement between the spatial join and ACLED's own admin3 label. This is a
# check, not a correction: disagreement is expected where the two gazetteers
# romanise differently, and the figure tells the reader how much.
norm   <- function(x) tolower(gsub("[^a-z]", "", tolower(x)))
has_a3 <- !is.na(ev_sf$admin3)
audit$admin3_agreement_pct <-
  100 * mean(norm(ev_sf$admin3[has_a3]) == norm(ev_sf$TS[has_a3]))
audit$admin3_disagree <- st_drop_geometry(ev_sf)[has_a3, ] |>
  filter(norm(admin3) != norm(TS)) |>
  count(acled_admin3 = admin3, mimu_township = TS, sort = TRUE)

# ---- What counts as an armed-conflict event -----------------------------------
# ACLED's six event types are not all armed conflict. Protests and riots are
# demonstrations; strategic developments record non-violent activity such as
# arrests, troop movements and territorial transfers. The three event types
# ACLED itself groups as political violence carry the analysis, and each is kept
# as its own variable because they are produced by different actors and should
# not be assumed to share a geography.
CONFLICT_TYPES <- c("Battles", "Explosions/Remote violence",
                    "Violence against civilians")

ev_sf <- ev_sf |>
  mutate(
    is_conflict = event_type %in% CONFLICT_TYPES,
    # Two location problems make a record unusable at township level even
    # though the spatial join succeeds. Records whose admin3 is missing were
    # placed by ACLED at a city reference point ("Mandalay - Unspecified",
    # "Yangon - Unspecified", "Nay Pyi Taw"): the city is known, the township
    # is not, and the join puts every one of them into whichever township holds
    # the reference point. Records at geo_precision 3 are known only to a
    # district or state and are coded at its capital. Both would manufacture a
    # concentration that the underlying events do not have.
    loc_city_only   = is.na(admin3),
    loc_region_only = geo_precision == 3,
    usable_loc      = !loc_city_only & !loc_region_only
  )

audit$city_only <- st_drop_geometry(ev_sf) |>
  filter(loc_city_only) |>
  count(location, mimu_township = TS, conflict = is_conflict, sort = TRUE)

audit$loc_flow <- st_drop_geometry(ev_sf) |>
  filter(is_conflict) |>
  summarise(
    conflict_geocoded = n(),
    city_only         = sum(loc_city_only),
    region_only       = sum(loc_region_only & !loc_city_only),
    analysed          = sum(usable_loc),
    fatal_geocoded    = sum(fatalities),
    fatal_analysed    = sum(fatalities[usable_loc])
  )

audit$exclusion_by_type <- st_drop_geometry(ev_sf) |>
  filter(is_conflict) |>
  group_by(event_type) |>
  summarise(geocoded = n(), city_only = sum(loc_city_only),
            region_only = sum(loc_region_only & !loc_city_only),
            analysed = sum(usable_loc),
            prec1 = sum(usable_loc & geo_precision == 1),
            prec2 = sum(usable_loc & geo_precision == 2),
            .groups = "drop")

# "Excessive force against protesters" is lethal state violence, but ACLED
# files it under the Protests event type. It is counted here so that the report
# can say how much of the 2021 crackdown falls outside the conflict variables.
audit$excessive_force <- st_drop_geometry(ev_sf) |>
  filter(sub_event_type == "Excessive force against protesters") |>
  summarise(n = n(), fatalities = sum(fatalities),
            first = min(event_date), last = max(event_date))

conf <- ev_sf |> filter(is_conflict, usable_loc)

audit$time_precision_conf <- table(conf$time_precision)

# ---- The space-time panels ---------------------------------------------------
# A space-time cube has to be complete: every township present in every period,
# including the townships and periods with no recorded event. Those zeroes are
# data -- a township with no reported violence in a quarter is exactly the
# observation a cold spot is built from -- so each panel is completed against
# the full grid rather than left as the sparse set of rows ACLED supplies.
# Quarters are the primary time unit and months are kept for the sensitivity
# analysis; both run over the 329 townships of the analytical frame.
make_panel <- function(d, time_col, periods) {
  d |>
    st_drop_geometry() |>
    group_by(TS_PCODE, period = .data[[time_col]]) |>
    summarise(
      conflict   = n(),
      battles    = sum(event_type == "Battles"),
      erv        = sum(event_type == "Explosions/Remote violence"),
      vac        = sum(event_type == "Violence against civilians"),
      fatalities = sum(fatalities),
      airstrike  = sum(sub_event_type == "Air/drone strike"),
      .groups = "drop"
    ) |>
    complete(TS_PCODE = ts_an$TS_PCODE, period = periods,
             fill = list(conflict = 0L, battles = 0L, erv = 0L, vac = 0L,
                         fatalities = 0, airstrike = 0L)) |>
    arrange(period, TS_PCODE)
}

months   <- seq(floor_date(STUDY_START, "month"),   floor_date(STUDY_END, "month"),   by = "month")
quarters <- seq(floor_date(STUDY_START, "quarter"), floor_date(STUDY_END, "quarter"), by = "quarter")

panel_m <- make_panel(conf, "ym", months)
panel_q <- make_panel(conf, "yq", quarters)

stopifnot(nrow(panel_m) == nrow(ts_an) * length(months),
          nrow(panel_q) == nrow(ts_an) * length(quarters),
          sum(panel_q$conflict) == nrow(conf),
          sum(panel_m$conflict) == nrow(conf))

# Whole-period township totals: the input to the local measures of spatial
# association.
tot <- conf |>
  st_drop_geometry() |>
  group_by(TS_PCODE) |>
  summarise(conflict   = n(),
            battles    = sum(event_type == "Battles"),
            erv        = sum(event_type == "Explosions/Remote violence"),
            vac        = sum(event_type == "Violence against civilians"),
            fatalities = sum(fatalities),
            airstrike  = sum(sub_event_type == "Air/drone strike"),
            .groups = "drop")

ts_tot <- ts_an |>
  left_join(tot, by = "TS_PCODE") |>
  mutate(across(c(conflict, battles, erv, vac, fatalities, airstrike),
                ~ tidyr::replace_na(.x, 0)),
         area_km2 = as.numeric(st_area(geometry)) / 1e6)

# Calendar years are kept as well, for the year-by-year LISA. 2025 covers nine
# months only, so any comparison of levels across years is made per month.
panel_y <- conf |>
  st_drop_geometry() |>
  mutate(year = year(event_date)) |>
  count(TS_PCODE, year, name = "conflict") |>
  complete(TS_PCODE = ts_an$TS_PCODE, year = 2021:2025,
           fill = list(conflict = 0L)) |>
  mutate(months = ifelse(year == 2025, 9L, 12L))

# ---- Aggregates for the report --------------------------------------------------
# ACLED's terms of use do not allow the event records themselves to be
# republished, so everything the report draws from is aggregated here: no
# derived object read by the report holds an individual event.
# Location precision by township: the share of each township's events that
# ACLED could only place "near" a town rather than at a named village.
ts_prec <- conf |>
  st_drop_geometry() |>
  group_by(TS_PCODE) |>
  summarise(n = n(), share_prec2 = mean(geo_precision == 2), .groups = "drop")
ts_tot <- ts_tot |> left_join(ts_prec |> select(TS_PCODE, share_prec2), by = "TS_PCODE")

# Who initiates each kind of event, by year. ACLED's inter1 code classifies the
# first actor; it is the closest the data come to attributing an event.
audit$actor_year <- conf |>
  st_drop_geometry() |>
  mutate(year = year(event_date)) |>
  count(year, event_type, actor_type = inter1)

audit$sub_event_year <- conf |>
  st_drop_geometry() |>
  mutate(year = year(event_date)) |>
  count(year, event_type, sub_event_type)

# State-level totals by quarter and event type, for the trend figures.
audit$state_quarter <- conf |>
  st_drop_geometry() |>
  count(ST, yq, event_type)

audit$n_months   <- length(months)
audit$n_quarters <- length(quarters)
audit$pct_zero_q <- 100 * mean(panel_q$conflict == 0)
audit$pct_zero_m <- 100 * mean(panel_m$conflict == 0)
audit$n_ts_zero  <- sum(ts_tot$conflict == 0)

# ---- Save --------------------------------------------------------------------
saveRDS(ev_sf,   file.path(DER, "events_sf.rds"))
saveRDS(conf,    file.path(DER, "conflict_sf.rds"))
saveRDS(panel_q, file.path(DER, "panel_q.rds"))
saveRDS(panel_m, file.path(DER, "panel_m.rds"))
saveRDS(panel_y, file.path(DER, "panel_y.rds"))
saveRDS(ts_tot,  file.path(DER, "ts_tot.rds"))
saveRDS(audit,   file.path(DER, "acled_audit.rds"))
unlink(file.path(DER, "panel.rds"))

cat("\n== Source ==\n"); cat(audit$source_file, " extracted", format(audit$extracted), "\n")
cat("\n== Record flow ==\n")
cat(sprintf("%-28s %8d\n", "raw rows",             audit$n_raw))
cat(sprintf("%-28s %8d\n", "not Myanmar",          audit$n_not_myanmar))
cat(sprintf("%-28s %8d\n", "after date/coord/dup", audit$n_after_filter))
cat(sprintf("%-28s %8d\n", "geocoded to township", audit$n_geocoded))
cat(sprintf("%-28s %8d\n", "outside any polygon",  audit$n_outside_poly))
cat(sprintf("%-28s %8d\n", "snapped within 5 km",  audit$n_snapped))
cat(sprintf("%-28s %8d\n", "in excluded township", audit$n_in_excluded))
cat("\n== Conflict events, location usability ==\n"); print(audit$loc_flow)
print(audit$exclusion_by_type)
cat("\n== admin3 agreement ==", round(audit$admin3_agreement_pct, 3), "%\n")
print(audit$admin3_disagree)
cat("\n== Excessive force against protesters ==\n"); print(audit$excessive_force)
cat("\n== Panels ==\n")
cat(sprintf("quarterly: %d x %d, %.1f%% zero | monthly: %d x %d, %.1f%% zero | townships with no conflict event: %d\n",
            nrow(ts_an), audit$n_quarters, audit$pct_zero_q,
            nrow(ts_an), audit$n_months, audit$pct_zero_m, audit$n_ts_zero))
print(summary(st_drop_geometry(ts_tot)[, c("conflict", "battles", "erv", "vac", "fatalities")]))
