suppressPackageStartupMessages({library(sf); library(dplyr); library(readr)})
sf_use_s2(FALSE)
raw <- "data/raw"; dir.create("data/derived", showWarnings=FALSE, recursive=TRUE)

adm1 <- st_read(file.path(raw,"geoBoundaries-IDN-ADM1.geojson"), quiet=TRUE) |>
  filter(grepl("Kalimantan", shapeName)) |> st_make_valid()
cat("Kalimantan provinces:", nrow(adm1), "\n")

if (!file.exists("data/derived/kal_adm2.rds")) {
  adm2 <- st_read(file.path(raw,"geoBoundaries-IDN-ADM2.geojson"), quiet=TRUE)
  cat("ADM2 total:", nrow(adm2), "| cols:", paste(names(adm2), collapse=", "), "\n")
  adm2 <- st_make_valid(adm2)
  kalbox <- st_as_sfc(st_bbox(adm1)); st_crs(kalbox) <- st_crs(adm2)
  cand <- adm2[lengths(st_intersects(adm2, kalbox))>0, ]
  ctr  <- st_point_on_surface(cand)
  j    <- st_join(ctr, adm1[,"shapeName"], join=st_within, suffix=c("",".prov"))
  cand$province <- j$shapeName.prov
  kal <- cand |> filter(!is.na(province))
  saveRDS(kal, "data/derived/kal_adm2.rds")
} else kal <- readRDS("data/derived/kal_adm2.rds")
cat("Kalimantan regencies:", nrow(kal), "\n")

fires <- read_csv(file.path(raw,"viirs_snpp_7d_seasia.csv"), show_col_types=FALSE) |>
  filter(longitude>107, longitude<120, latitude>-5, latitude<3)
fp <- st_as_sf(fires, coords=c("longitude","latitude"), crs=4326)
idx <- st_within(fp, kal, sparse=TRUE)
fp$reg <- sapply(idx, function(i) if(length(i)) i[1] else NA_integer_)

kal$area_km2 <- as.numeric(st_area(st_transform(kal, 23845)))/1e6
tab <- fp |> st_drop_geometry() |> filter(!is.na(reg)) |>
  count(reg, name="n_fires") |>
  mutate(regency=kal$shapeName[reg], province=kal$province[reg],
         area_km2=round(kal$area_km2[reg],0),
         per_1000km2=round(n_fires/area_km2*1000,1)) |>
  arrange(desc(n_fires)) |> select(regency, province, n_fires, area_km2, per_1000km2)
cat("\n==== TOP 20 REGENCIES BY VIIRS S-NPP DETECTIONS (2026-09-04..11) ====\n")
print(as.data.frame(head(tab,20)), row.names=FALSE)
write_csv(tab, "data/derived/regency_fire_counts.csv")
