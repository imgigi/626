# =============================================================================
# 01b_assemble_public7d.R
# Take-home Exercise 1 | ISSS626
#
# PURPOSE
#   Fallback / cross-check acquisition path. FIRMS publishes the most recent
#   7 days of detections as open regional CSVs that need no MAP_KEY:
#     https://firms.modaps.eosdis.nasa.gov/data/active_fire/<product>/csv/
#   This script harvests the "SouthEast_Asia" 7-day files for four sensors and
#   rewrites them into exactly the schema produced by 01_download_firms.R, so
#   that the downstream preparation and analysis code is identical whichever
#   acquisition route supplied the data.
#
#   Note the products differ: the 7-day open files are NRT only and therefore
#   carry no `type` flag (static-source / volcano / offshore classification).
#   That limitation is carried forward and reported in the quality assessment.
# =============================================================================

suppressPackageStartupMessages({library(readr); library(dplyr); library(glue)})

RAW <- "data/raw"; dir.create(RAW, recursive = TRUE, showWarnings = FALSE)

PRODUCTS <- tibble::tribble(
  ~firms_source,        ~path,
  "VIIRS_SNPP_NRT",     "suomi-npp-viirs-c2/csv/SUOMI_VIIRS_C2_SouthEast_Asia_7d.csv",
  "VIIRS_NOAA20_NRT",   "noaa-20-viirs-c2/csv/J1_VIIRS_C2_SouthEast_Asia_7d.csv",
  "VIIRS_NOAA21_NRT",   "noaa-21-viirs-c2/csv/J2_VIIRS_C2_SouthEast_Asia_7d.csv",
  "MODIS_NRT",          "modis_c6_1/csv/MODIS_C6_1_SouthEast_Asia_7d.csv"
)
BASE <- "https://firms.modaps.eosdis.nasa.gov/data/active_fire"
BBOX <- c(west = 113.43, south = -3.62, east = 114.50, north = -1.39)

out <- purrr::pmap_dfr(PRODUCTS, function(firms_source, path) {
  dest <- file.path(RAW, basename(path))
  if (!file.exists(dest)) {
    message("downloading ", basename(path))
    utils::download.file(file.path(BASE, path), dest, mode = "wb", quiet = TRUE)
  }
  d <- read_csv(dest, show_col_types = FALSE, col_types = cols(.default = col_character()))
  d$firms_source <- firms_source
  # Clip to the buffered study bbox at read time: the regional file covers all
  # of South-East Asia and is ~5 MB per sensor.
  d |> mutate(.lon = as.numeric(longitude), .lat = as.numeric(latitude)) |>
    filter(.lon >= BBOX["west"], .lon <= BBOX["east"],
           .lat >= BBOX["south"], .lat <= BBOX["north"]) |>
    select(-.lon, -.lat)
})

rng <- range(as.Date(out$acq_date))
dest <- glue("{RAW}/firms_pulangpisau_{rng[1]}_{rng[2]}.csv")
write_csv(out, dest)
message(glue("\n{nrow(out)} records, {rng[1]} to {rng[2]} -> {dest}"))
print(out |> count(firms_source, name = "records"))
print(out |> count(firms_source, acq_date) |> tidyr::pivot_wider(names_from = firms_source, values_from = n))
