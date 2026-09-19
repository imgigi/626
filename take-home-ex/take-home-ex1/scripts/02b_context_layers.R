# =============================================================================
# 02b_context_layers.R -- small committed inputs for the context map
# Take-home Exercise 1 | ISSS626
#
# The national COD-AB layers are ~230 MB unpacked and are therefore not
# committed. This script distils from them the handful of small objects the
# report actually needs as INPUTS, so that a fresh clone can render without
# re-downloading a national dataset:
#
#   pp_adm2_cod.rds       the study window (Pulang Pisau, full detail)
#   pp_adm3_cod.rds       its eight kecamatan (full detail)
#   kal_prov.rds          the five Kalimantan provinces, SIMPLIFIED to 1 km --
#                         it is a background layer on an island-scale map, so
#                         full coastline detail costs 8 MB for no visible gain
#   kal_fires_context.rds island-wide S-NPP detections for the same eight days
#
# Run after 01b (or 01) and before rendering. Requires the COD-AB extract in
# data/raw/idn_cod_ab/ and the regional FIRMS CSV in data/raw/.
# =============================================================================

suppressPackageStartupMessages({library(sf); library(dplyr); library(readr)})
sf_use_s2(FALSE)
DER <- "data/derived"; dir.create(DER, recursive = TRUE, showWarnings = FALSE)
CRS_PROJ <- 32749

a1 <- st_read("data/raw/idn_cod_ab/idn_admin1.shp", quiet = TRUE) |>
  filter(grepl("Kalimantan", adm1_name)) |> st_make_valid() |>
  select(adm1_name, adm1_pcode)

# 1 km simplification: invisible at the ~1,500 km scale of the inset, and it
# reduces the object from ~8 MB to well under 1 MB.
kal_prov <- a1 |> st_transform(CRS_PROJ) |>
  st_simplify(dTolerance = 1000, preserveTopology = TRUE) |> st_make_valid()
saveRDS(kal_prov, file.path(DER, "kal_prov.rds"))
message("kal_prov: ", nrow(kal_prov), " provinces, ",
        round(file.size(file.path(DER, "kal_prov.rds")) / 1024), " KB")

a2 <- st_read("data/raw/idn_cod_ab/idn_admin2.shp", quiet = TRUE) |> st_make_valid()
pp <- a2 |> filter(grepl("^Pulang Pisau$", adm2_name))
stopifnot(nrow(pp) == 1, st_is_valid(pp))
saveRDS(pp, file.path(DER, "pp_adm2_cod.rds"))

a3 <- st_read("data/raw/idn_cod_ab/idn_admin3.shp", quiet = TRUE) |> st_make_valid()
kec <- a3 |> filter(grepl("^Pulang Pisau$", adm2_name))
stopifnot(nrow(kec) == 8)
saveRDS(kec, file.path(DER, "pp_adm3_cod.rds"))
message("pp_adm2 + pp_adm3 (", nrow(kec), " kecamatan) written")

# Island-wide detections for the context panel, from the regional 7-day file.
f <- "data/raw/SUOMI_VIIRS_C2_SouthEast_Asia_7d.csv"
if (file.exists(f)) {
  fr <- read_csv(f, show_col_types = FALSE) |>
    filter(longitude > 108.5, longitude < 119.5, latitude > -4.5, latitude < 3,
           tolower(substr(confidence, 1, 1)) != "l")
  p <- st_as_sf(fr, coords = c("longitude", "latitude"), crs = 4326)
  p <- p[lengths(st_intersects(p, st_union(st_transform(kal_prov, 4326)))) > 0, ]
  saveRDS(st_transform(p, CRS_PROJ), file.path(DER, "kal_fires_context.rds"))
  message("kal_fires_context: ", nrow(p), " detections")
} else {
  message("regional CSV absent; kal_fires_context.rds left as is")
}
