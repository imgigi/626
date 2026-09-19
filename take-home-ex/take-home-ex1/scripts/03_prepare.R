# =============================================================================
# 03_prepare.R
# Take-home Exercise 1 | ISSS626
#
# PURPOSE
#   Turn raw FIRMS detection records + official administrative boundaries into
#   the analysis-ready objects used by the technical report, while recording a
#   complete, auditable account of every record that is dropped and why.
#
# OUTPUTS (data/derived/)
#   qa_audit.rds        record-level audit trail (step, n_in, n_out, reason)
#   fire_sf.rds         retained detections as sf POINT, EPSG:32749
#   fire_all_sf.rds     retained detections, all four sensors (sensitivity)
#   win.rds             study window as spatstat owin, EPSG:32749
#   ppp_det.rds         ppp: every retained detection
#   ppp_pxday.rds       ppp: unique 375 m pixel x local day  <-- PRIMARY UNIT
#   ppp_event.rds       ppp: spatio-temporally clustered fire events
# =============================================================================

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(readr); library(tidyr); library(purrr)
  library(lubridate); library(stringr); library(spatstat.geom)
})
sf_use_s2(FALSE)
set.seed(1234)

DER <- "data/derived"; dir.create(DER, recursive = TRUE, showWarnings = FALSE)

# ---- Analytical parameters (single place, so the report can cite them) ------
P <- list(
  crs_geo         = 4326,          # FIRMS delivers WGS84 geographic coordinates
  crs_proj        = 32749,         # WGS 84 / UTM 49S: metre units, conformal
  tz_local        = "Asia/Jakarta",# Central Kalimantan observes WIB = UTC+07
  primary_sensor  = "VIIRS_SNPP",  # single, consistent observation process
  drop_confidence = "low",         # VIIRS low class: elevated false-alarm rate
  pixel_grid_m    = 375,           # VIIRS I-band nominal resolution
  event_space_m   = 750,           # = 2 pixels: linkage radius for events
  event_time_day  = 1,             # detections <=1 day apart may be one event
  persist_days    = 6              # >= this many distinct days => review
)
saveRDS(P, file.path(DER, "params.rds"))

audit <- tibble(step = character(), n_in = integer(), n_out = integer(),
                n_removed = integer(), reason = character())
log_step <- function(step, n_in, n_out, reason) {
  audit <<- add_row(audit, step = step, n_in = n_in, n_out = n_out,
                    n_removed = n_in - n_out, reason = reason)
  invisible(NULL)
}

# ---- 1. Locate and import the raw detection file ---------------------------
# Whichever acquisition script ran (01 = MAP_KEY archive, 01b = open 7-day
# files), the raw product is one CSV named firms_pulangpisau_<from>_<to>.csv
raw_files <- list.files("data/raw", "^firms_pulangpisau_.*[.]csv$", full.names = TRUE)
stopifnot(length(raw_files) >= 1)
RAW_FILE <- raw_files[which.max(file.info(raw_files)$size)]
message("Using raw file: ", RAW_FILE)

# Read every column as character first: this prevents readr from silently
# coercing a malformed coordinate to NA before we can count it.
fires_raw <- read_csv(RAW_FILE, show_col_types = FALSE,
                      col_types = cols(.default = col_character()))
n0 <- nrow(fires_raw)
log_step("00 import", n0, n0, paste("raw records read from", basename(RAW_FILE)))

# ---- 2. Harmonise the VIIRS / MODIS schemas --------------------------------
# VIIRS reports bright_ti4 / bright_ti5 (4 um and 11 um channel brightness
# temperatures); MODIS reports brightness / bright_t31. They are renamed to a
# common pair so the sensitivity analysis can pool sensors, but the originating
# product is always retained in `firms_source`.
has_col <- function(nm) nm %in% names(fires_raw)
getcol  <- function(nm) if (has_col(nm)) fires_raw[[nm]] else NA_character_

fires <- fires_raw |>
  mutate(
    sensor_family  = if_else(str_detect(firms_source, "^VIIRS"), "VIIRS", "MODIS"),
    bt_fire        = suppressWarnings(as.numeric(coalesce(getcol("bright_ti4"),
                                                          getcol("brightness")))),
    bt_bg          = suppressWarnings(as.numeric(coalesce(getcol("bright_ti5"),
                                                          getcol("bright_t31")))),
    frp            = suppressWarnings(as.numeric(frp)),
    scan           = suppressWarnings(as.numeric(scan)),
    track          = suppressWarnings(as.numeric(track)),
    lon            = suppressWarnings(as.numeric(longitude)),
    lat            = suppressWarnings(as.numeric(latitude)),
    type           = suppressWarnings(as.integer(getcol("type"))),
    confidence_raw = confidence
  )

# MODIS confidence is a 0-100 integer; VIIRS confidence is l/n/h. Both are
# mapped onto one ordered class so a single filtering rule can be stated.
fires <- fires |>
  mutate(
    conf_num = suppressWarnings(as.numeric(confidence_raw)),
    confidence_cls = case_when(
      sensor_family == "VIIRS" & str_starts(str_to_lower(confidence_raw), "l") ~ "low",
      sensor_family == "VIIRS" & str_starts(str_to_lower(confidence_raw), "n") ~ "nominal",
      sensor_family == "VIIRS" & str_starts(str_to_lower(confidence_raw), "h") ~ "high",
      sensor_family == "MODIS" & conf_num <  30 ~ "low",
      sensor_family == "MODIS" & conf_num <  80 ~ "nominal",
      sensor_family == "MODIS" & conf_num >= 80 ~ "high",
      TRUE ~ NA_character_),
    confidence_cls = factor(confidence_cls,
                            levels = c("low", "nominal", "high"), ordered = TRUE))

# ---- 3. Missing and invalid coordinates ------------------------------------
n_miss_coord <- sum(is.na(fires$lon) | is.na(fires$lat))
fires <- fires |> filter(!is.na(lon), !is.na(lat))
log_step("01 missing coordinates", n0, nrow(fires),
         sprintf("%d record(s) with absent or non-numeric latitude/longitude",
                 n_miss_coord))

n <- nrow(fires)
fires <- fires |>
  filter(between(lat, -90, 90), between(lon, -180, 180), !(lon == 0 & lat == 0))
log_step("02 impossible coordinates", n, nrow(fires),
         "outside the valid geographic range, or the (0,0) null island")

# ---- 4. Timestamp construction, in LOCAL time ------------------------------
# FIRMS `acq_date` + `acq_time` are UTC. Central Kalimantan is UTC+07 and the
# VIIRS night overpass is ~01:30 local = ~18:30 UTC on the PREVIOUS calendar
# day. Taking the "day" of a fire from the UTC date would therefore misassign
# every night detection by one day and corrupt any daily time series. All
# temporal variables below are derived from the LOCAL timestamp.
n <- nrow(fires)
fires <- fires |>
  mutate(acq_time4 = str_pad(acq_time, 4, "left", "0"),
         dt_utc    = ymd_hm(paste(acq_date, acq_time4), tz = "UTC", quiet = TRUE),
         dt_local  = with_tz(dt_utc, P$tz_local)) |>
  filter(!is.na(dt_utc))
log_step("03 unparseable timestamp", n, nrow(fires),
         "acq_date/acq_time could not be parsed into a valid UTC instant")

fires <- fires |>
  mutate(date_local = as_date(dt_local),
         date_utc   = as_date(acq_date),
         day_shifted = date_local != date_utc,
         doy        = yday(date_local),
         month      = month(date_local, label = TRUE),
         iso_week   = isoweek(date_local),
         hour_local = hour(dt_local) + minute(dt_local) / 60,
         overpass   = factor(if_else(daynight == "D", "Day", "Night"),
                             levels = c("Day", "Night")))

# ---- 5. Static / non-vegetation sources ------------------------------------
# FIRMS Standard-Processing products carry a `type` flag:
#   0 presumed vegetation fire | 1 active volcano
#   2 other static land source | 3 offshore
# Only type 0 is a candidate wildland fire. The open 7-day NRT files do NOT
# carry this flag; when it is absent the filter cannot be applied, and the
# limitation is logged rather than silently ignored.
n <- nrow(fires)
if (all(is.na(fires$type))) {
  log_step("04 static-source flag", n, n,
           "NOT APPLIED: `type` flag absent from NRT product; persistent-pixel screening used instead")
} else {
  fires <- fires |> filter(is.na(type) | type == 0)
  log_step("04 static-source flag", n, nrow(fires),
           "removed volcano (1), other static land source (2) and offshore (3) detections")
}

# ---- 6. Confidence screening -----------------------------------------------
n <- nrow(fires)
fires <- fires |> filter(is.na(confidence_cls) | confidence_cls != P$drop_confidence)
log_step("05 low confidence", n, nrow(fires),
         "removed 'low' confidence detections (elevated false-alarm rate)")

# ---- 7. Exact duplicate records --------------------------------------------
n <- nrow(fires)
fires <- fires |>
  distinct(firms_source, lon, lat, acq_date, acq_time, .keep_all = TRUE)
log_step("06 exact duplicates", n, nrow(fires),
         "identical product + coordinate + acquisition timestamp (repeat rows across download chunks)")

# ---- 8. Project, then clip to the administrative study window --------------
bnd <- readRDS(file.path(DER, "pp_adm2_cod.rds"))
stopifnot(all(st_is_valid(bnd)))
bnd_proj <- st_transform(bnd, P$crs_proj)

fire_all_sf <- st_as_sf(fires, coords = c("lon", "lat"),
                        crs = P$crs_geo, remove = FALSE) |>
  st_transform(P$crs_proj)

n <- nrow(fire_all_sf)
inside <- lengths(st_intersects(fire_all_sf, bnd_proj)) > 0
fire_all_sf <- fire_all_sf[inside, ]
log_step("07 outside study window", n, nrow(fire_all_sf),
         "detections inside the buffered download box but outside Pulang Pisau Regency")

# ---- 9. Snap to the VIIRS pixel grid; screen persistent pixels -------------
# Repeat detections of the same ground pixel are the dominant source of
# *apparent* small-scale clustering. Snapping to a 375 m grid makes the
# repetition explicit and countable.
xy <- st_coordinates(fire_all_sf)
fire_all_sf <- fire_all_sf |>
  mutate(gx = floor(xy[, 1] / P$pixel_grid_m),
         gy = floor(xy[, 2] / P$pixel_grid_m),
         pixel_id = paste(gx, gy, sep = "_"))

persist <- fire_all_sf |>
  st_drop_geometry() |>
  filter(str_detect(firms_source, P$primary_sensor)) |>
  group_by(pixel_id) |>
  summarise(n_days = n_distinct(date_local), n_det = n(), .groups = "drop") |>
  arrange(desc(n_days), desc(n_det))
saveRDS(persist, file.path(DER, "persistent_pixels.rds"))

# ---- 10. Primary sensor subset ---------------------------------------------
fire_sf <- fire_all_sf |> filter(str_detect(firms_source, P$primary_sensor))
log_step("08 primary sensor", nrow(fire_all_sf), nrow(fire_sf),
         sprintf("restricted to %s for the main analysis; other products retained for sensitivity testing",
                 P$primary_sensor))

# ---- 11. Build the study window (owin) -------------------------------------
win <- as.owin(st_geometry(bnd_proj))
stopifnot(inherits(win, "owin"))

# ---- 12. Three analysis units ----------------------------------------------
mk_ppp <- function(d, win) {
  ppp(x = st_coordinates(d)[, 1], y = st_coordinates(d)[, 2],
      window = win, check = FALSE)
}

# (a) every retained detection
ppp_det <- mk_ppp(fire_sf, win)

# (b) PRIMARY UNIT: one point per 375 m pixel per local day. This removes the
#     purely observational duplication created when two overpasses on the same
#     day see the same pixel, while keeping genuine day-to-day persistence.
fire_pxday <- fire_sf |>
  group_by(pixel_id, date_local) |>
  slice_max(frp, n = 1, with_ties = FALSE) |>   # keep most energetic record
  ungroup()
ppp_pxday <- mk_ppp(fire_pxday, win)

# (c) spatio-temporal event clustering: single-linkage over a space-time
#     neighbourhood (<= event_space_m AND <= event_time_day), giving an
#     approximate count of distinct *fire events* rather than detections.
cl_xy <- st_coordinates(fire_pxday)
nb  <- dbscan::frNN(cl_xy[, 1:2], eps = P$event_space_m)
i   <- rep(seq_along(nb$id), lengths(nb$id))
j   <- unlist(nb$id)
if (length(i) == 0) {
  comp <- seq_len(nrow(fire_pxday))
} else {
  dt   <- abs(as.numeric(fire_pxday$date_local[i] - fire_pxday$date_local[j]))
  keep <- dt <= P$event_time_day
  g <- igraph::make_empty_graph(n = nrow(fire_pxday), directed = FALSE)
  if (any(keep)) g <- igraph::add_edges(g, as.vector(rbind(i[keep], j[keep])))
  comp <- igraph::components(g)$membership
}
fire_pxday$event_id <- comp

fire_event <- fire_pxday |>
  st_drop_geometry() |>
  mutate(Xm = cl_xy[, 1], Ym = cl_xy[, 2]) |>
  group_by(event_id) |>
  summarise(x = mean(Xm), y = mean(Ym),
            start = min(date_local), end = max(date_local),
            n_px = n(), frp_sum = sum(frp, na.rm = TRUE),
            frp_max = suppressWarnings(max(frp, na.rm = TRUE)),
            .groups = "drop") |>
  mutate(dur_days = as.numeric(end - start) + 1)
ppp_event <- ppp(fire_event$x, fire_event$y, window = win, check = FALSE)

# ---- 13. Duplicated coordinates in the ppp objects -------------------------
# Several spatstat summary functions assume a simple point pattern. Remaining
# coincident points are reported, so the report states the fact rather than
# leaving it buried in a warning.
dup_report <- tibble(
  unit = c("detections", "pixel-day", "events"),
  n_points = c(npoints(ppp_det), npoints(ppp_pxday), npoints(ppp_event)),
  n_duplicated = c(sum(duplicated(ppp_det)), sum(duplicated(ppp_pxday)),
                   sum(duplicated(ppp_event))),
  n_outside_window = c(sum(!inside.owin(ppp_det$x, ppp_det$y, win)),
                       sum(!inside.owin(ppp_pxday$x, ppp_pxday$y, win)),
                       sum(!inside.owin(ppp_event$x, ppp_event$y, win)))
)

# ---- 14. Persist -----------------------------------------------------------
saveRDS(fire_all_sf, file.path(DER, "fire_all_sf.rds"))
saveRDS(fire_sf,     file.path(DER, "fire_sf.rds"))
saveRDS(fire_pxday,  file.path(DER, "fire_pxday_sf.rds"))
saveRDS(fire_event,  file.path(DER, "fire_event.rds"))
saveRDS(win,         file.path(DER, "win.rds"))
saveRDS(bnd_proj,    file.path(DER, "bnd_proj.rds"))
saveRDS(ppp_det,     file.path(DER, "ppp_det.rds"))
saveRDS(ppp_pxday,   file.path(DER, "ppp_pxday.rds"))
saveRDS(ppp_event,   file.path(DER, "ppp_event.rds"))
saveRDS(audit,       file.path(DER, "qa_audit.rds"))
saveRDS(dup_report,  file.path(DER, "dup_report.rds"))

message("\n================ QA AUDIT TRAIL ================")
print(as.data.frame(audit), row.names = FALSE)
message("\n================ ANALYSIS UNITS ================")
print(as.data.frame(dup_report), row.names = FALSE)
message("\n========= PERSISTENT PIXELS (top 10) ===========")
print(as.data.frame(head(persist, 10)), row.names = FALSE)
message("\nLocal date range: ", paste(range(fire_sf$date_local), collapse = " .. "))
message("Records whose LOCAL date differs from the UTC date: ",
        sum(fire_sf$day_shifted), " of ", nrow(fire_sf))
message("\nDay/night composition:"); print(table(fire_sf$overpass))
message("\nConfidence composition:"); print(table(fire_sf$confidence_cls))
message("\nEvents: ", nrow(fire_event),
        " | median pixels per event: ", median(fire_event$n_px),
        " | max: ", max(fire_event$n_px))
