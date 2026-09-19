# =============================================================================
# 01_download_firms.R
# Take-home Exercise 1 | ISSS626 Geospatial Analytics and Applications
#
# PURPOSE
#   Acquire NASA FIRMS active-fire detections for the Pulang Pisau study area
#   for the 2026 observation period, via the FIRMS Area API.
#
# WHY A SCRIPT (AND NOT AN INLINE CHUNK IN THE REPORT)
#   The FIRMS Area API requires a personal MAP_KEY and is rate limited.
#   Acquisition is therefore separated from analysis so that the technical
#   report can be re-rendered offline from the cached raw CSVs, while the
#   acquisition step remains fully documented and reproducible.
#
# AUTHENTICATION
#   Obtain a free MAP_KEY at https://firms.modaps.eosdis.nasa.gov/api/map_key/
#   Supply it in ONE of the following ways (checked in this order):
#     1. environment variable  FIRMS_MAP_KEY
#     2. a one-line file       .firms_key   (git-ignored)
#
# USAGE
#   Rscript scripts/01_download_firms.R [start_date] [end_date]
#   e.g. Rscript scripts/01_download_firms.R 2026-01-01 2026-09-11
# =============================================================================

suppressPackageStartupMessages({
  library(httr); library(readr); library(dplyr); library(purrr); library(glue)
})

# ---- 0. Parameters ----------------------------------------------------------

args       <- commandArgs(trailingOnly = TRUE)
START_DATE <- if (length(args) >= 1) as.Date(args[1]) else as.Date("2026-01-01")
END_DATE   <- if (length(args) >= 2) as.Date(args[2]) else Sys.Date() - 1

# Study-area bounding box (WGS84), = Pulang Pisau Regency bbox buffered by
# ~0.15 deg (~17 km) so that (a) points just outside the regency are available
# for boundary-sensitivity checks, and (b) the exact administrative clip is
# performed later, in the analysis, rather than being baked into the download.
BBOX <- c(west = 113.43, south = -3.62, east = 114.50, north = -1.39)

# Sensors. VIIRS 375 m is the primary instrument family for this study:
#   VIIRS_SNPP_*  Suomi-NPP   (375 m, operational since 2012)  -> PRIMARY
#   VIIRS_NOAA20_NRT / VIIRS_NOAA21_NRT                        -> sensitivity
#   MODIS_*       Terra/Aqua  (1 km, longest record)            -> sensitivity
# "_SP" = Standard Processing (science quality, includes the `type` flag);
# "_NRT" = Near Real-Time (last ~2-3 months, no `type` flag).
SOURCES <- c("VIIRS_SNPP_SP",   "VIIRS_SNPP_NRT",
             "VIIRS_NOAA20_NRT","VIIRS_NOAA21_NRT",
             "MODIS_SP",        "MODIS_NRT")

CHUNK_DAYS <- 10L   # hard maximum accepted by the FIRMS Area API
OUT_DIR    <- "data/raw/firms_chunks"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Resolve the MAP_KEY -------------------------------------------------

get_map_key <- function() {
  k <- Sys.getenv("FIRMS_MAP_KEY", unset = "")
  if (nzchar(k)) return(trimws(k))
  if (file.exists(".firms_key")) {
    k <- trimws(readLines(".firms_key", warn = FALSE)[1])
    if (nzchar(k)) return(k)
  }
  stop("No FIRMS MAP_KEY found. Set env var FIRMS_MAP_KEY or create .firms_key\n",
       "Request one (free) at https://firms.modaps.eosdis.nasa.gov/api/map_key/",
       call. = FALSE)
}
MAP_KEY <- get_map_key()

# ---- 2. Interrogate data availability --------------------------------------
# Each FIRMS source exposes a different valid date range (SP lags real time by
# ~2-3 months; NRT covers only the recent window). Querying availability first
# avoids wasted requests and lets us document exactly which sensor-product
# supplied each part of the time series.

availability <- function() {
  u <- glue("https://firms.modaps.eosdis.nasa.gov/api/data_availability/csv/{MAP_KEY}/all")
  r <- GET(u, timeout(60))
  stop_for_status(r)
  txt <- content(r, "text", encoding = "UTF-8")
  if (grepl("Invalid MAP_KEY", txt, fixed = TRUE)) stop("FIRMS rejected the MAP_KEY.", call. = FALSE)
  read_csv(txt, show_col_types = FALSE)
}

avail <- availability()
message("\n== FIRMS data availability ==")
print(as.data.frame(avail))
write_csv(avail, "data/raw/firms_data_availability.csv")

# ---- 3. Chunked download ----------------------------------------------------

fetch_chunk <- function(source, from, n_days) {
  dest <- file.path(OUT_DIR, glue("{source}_{from}_{n_days}d.csv"))
  if (file.exists(dest) && file.info(dest)$size > 0) return(dest)  # cache hit

  u <- glue("https://firms.modaps.eosdis.nasa.gov/api/area/csv/{MAP_KEY}/{source}/",
            "{BBOX['west']},{BBOX['south']},{BBOX['east']},{BBOX['north']}/",
            "{n_days}/{from}")

  for (attempt in 1:4) {                     # simple backoff: API is rate limited
    r <- try(GET(u, timeout(120)), silent = TRUE)
    if (!inherits(r, "try-error") && status_code(r) == 200) {
      txt <- content(r, "text", encoding = "UTF-8")
      if (grepl("Invalid", substr(txt, 1, 40))) {
        warning(glue("{source} {from}: {substr(txt,1,80)}")); return(NA_character_)
      }
      writeLines(txt, dest); return(dest)
    }
    Sys.sleep(5 * attempt)
  }
  warning(glue("Failed after retries: {source} {from}")); NA_character_
}

# Build the (source, chunk-start, length) request grid, restricted to each
# source's advertised availability window.
starts <- seq(START_DATE, END_DATE, by = glue("{CHUNK_DAYS} days"))
grid <- map_dfr(SOURCES, function(src) {
  row <- avail[avail[[1]] == src, ]
  if (nrow(row) == 0) { message(glue("  - {src}: not advertised; skipped")); return(NULL) }
  lo <- as.Date(row$min_date[1]); hi <- as.Date(row$max_date[1])
  tibble(source = src, from = starts) |>
    mutate(n_days = pmin(CHUNK_DAYS, as.integer(END_DATE - from) + 1L)) |>
    filter(from >= lo, from <= hi, n_days > 0)
})

message(glue("\n== Requesting {nrow(grid)} chunks ({length(unique(grid$source))} sources) =="))
grid$file <- pmap_chr(list(grid$source, grid$from, grid$n_days), fetch_chunk)

# ---- 4. Combine -------------------------------------------------------------
# Chunks are concatenated per sensor family. Columns differ between VIIRS
# (bright_ti4/bright_ti5) and MODIS (brightness/bright_t31) and between SP
# (has `type`) and NRT (no `type`), so binding is done with type-tolerant
# coercion and the provenance of every row is retained.

read_chunk <- function(f, src) {
  if (is.na(f)) return(NULL)
  d <- suppressWarnings(read_csv(f, show_col_types = FALSE,
                                 col_types = cols(.default = col_character())))
  if (nrow(d) == 0) return(NULL)
  d$firms_source <- src
  d
}

all_rows <- map2(grid$file, grid$source, read_chunk) |> compact()
if (!length(all_rows)) stop("No FIRMS rows downloaded.", call. = FALSE)
fires <- bind_rows(all_rows)

out <- glue("data/raw/firms_pulangpisau_{START_DATE}_{END_DATE}.csv")
write_csv(fires, out)

message(glue("\n== Done: {nrow(fires)} raw detection records -> {out} =="))
print(fires |> count(firms_source, name = "records"))
