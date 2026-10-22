# =============================================================================
# 01_boundaries.R -- Myanmar township frame for the analysis
# Take-home Exercise 2 | ISSS626
#
# Reads the MIMU v9.4 township polygons, establishes a working projection by
# measuring distortion rather than assuming one, and records the two
# properties that decide whether a contiguity-based local statistic can be
# computed at all: duplicate identifiers and townships with no neighbours.
# Run from take-home-ex/take-home-ex2/.
# =============================================================================

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(spdep)
})

RAW <- "data/raw"; DER <- "data/derived"
dir.create(DER, recursive = TRUE, showWarnings = FALSE)

# ---- Import -----------------------------------------------------------------
# Read straight from the archive as downloaded from the MIMU GeoNode, so that
# the zip is the only boundary file the repository needs to carry. MIMU replaces
# its boundary releases, and v9.4 may not stay downloadable.
ts <- st_read(paste0("/vsizip/", RAW, "/mmr_polbnda_adm3_250k_mimu_v9_4.zip/",
                     "mmr_polbnda_adm3_250k_mimu_1.shp"), quiet = TRUE)

stopifnot(nrow(ts) == 330L)

# ---- Geometry and identifier integrity --------------------------------------
# A local statistic is indexed by polygon, so a duplicated or missing township
# pcode silently corrupts every join downstream. Both are checked before the
# layer is used for anything.
geom_report <- list(
  n          = nrow(ts),
  invalid    = sum(!st_is_valid(ts)),
  empty      = sum(st_is_empty(ts)),
  dup_pcode  = sum(duplicated(ts$TS_PCODE)),
  na_pcode   = sum(is.na(ts$TS_PCODE) | ts$TS_PCODE == ""),
  multipart  = sum(vapply(st_geometry(ts), function(g) length(g) > 1L, logical(1)))
)

if (geom_report$invalid > 0) ts <- st_make_valid(ts)

# ---- Choosing the working projection ----------------------------------------
# Myanmar spans 92.17-101.17 E, which straddles UTM zones 46N and 47N; neither
# zone covers the country without stretching. Rather than adopt one by
# convention, each candidate is scored on how far its planar area departs from
# the geodesic area of the same polygons, which is the quantity an areal
# analysis actually depends on.
# The reference must be the ELLIPSOIDAL area. With s2 enabled sf measures on a
# sphere, which is itself 0.27% away from the ellipsoid and would mask the very
# distortion being measured; with s2 disabled sf routes geographic coordinates
# to lwgeom's geodesic area, which is computed on WGS 84 itself.
sf_use_s2(FALSE)
area_geod <- as.numeric(st_area(st_transform(ts, 4326)))   # m^2, ellipsoidal

candidates <- c(
  "EPSG:32646" = "EPSG:32646",                 # WGS 84 / UTM 46N
  "EPSG:32647" = "EPSG:32647",                 # WGS 84 / UTM 47N
  # Lambert azimuthal equal area centred on the country's own centroid
  "LAEA_MMR"   = paste0("+proj=laea +lat_0=19.1 +lon_0=96.7 ",
                        "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"),
  # Albers equal area with standard parallels set inside Myanmar's latitude span
  "AEA_MMR"    = paste0("+proj=aea +lat_1=12 +lat_2=26 +lat_0=19.1 +lon_0=96.7 ",
                        "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
)

crs_scores <- lapply(names(candidates), function(nm) {
  a <- as.numeric(st_area(st_transform(ts, candidates[[nm]])))
  data.frame(crs = nm,
             total_err_pct = 100 * abs(sum(a) - sum(area_geod)) / sum(area_geod),
             max_abs_err_pct = max(100 * abs(a - area_geod) / area_geod),
             med_abs_err_pct = median(100 * abs(a - area_geod) / area_geod))
}) |> bind_rows() |> arrange(max_abs_err_pct)

WORK_CRS <- candidates[[crs_scores$crs[1]]]
ts_proj  <- st_transform(ts, WORK_CRS)

# ---- Contiguity ---------------------------------------------------------------
# Queen contiguity is the default neighbour definition for township-level LMSA.
# Myanmar's coast includes island townships, and an island has no queen
# neighbour at all: spdep returns a zero-length neighbour set, which propagates
# as NaN through any row-standardised weight and quietly removes those
# townships from the analysis. They are identified here so the decision about
# them is made explicitly in the report rather than by accident.
nb_queen <- poly2nb(ts_proj, queen = TRUE)
no_nb    <- which(card(nb_queen) == 0L)

islands <- if (length(no_nb)) {
  st_drop_geometry(ts_proj)[no_nb, c("ST", "DT", "TS", "TS_PCODE")]
} else {
  st_drop_geometry(ts_proj)[0, c("ST", "DT", "TS", "TS_PCODE")]
}

nb_report <- list(
  n_links      = sum(card(nb_queen)),
  avg_links    = mean(card(nb_queen)),
  min_links    = min(card(nb_queen)),
  max_links    = max(card(nb_queen)),
  n_no_nb      = length(no_nb),
  n_comp       = n.comp.nb(nb_queen)$nc,
  pct_nonzero  = 100 * sum(card(nb_queen) > 0) / nrow(ts_proj)
)

# ---- The analytical frame and its neighbour structure ------------------------
# Each neighbourless township is resolved on its own merits rather than by one
# blanket rule. The edge-to-edge distance to the nearest township decides it:
# Munaung and Chaungzon lie within a few kilometres of the mainland and are
# linked to their single nearest township, which is how they are reached in
# practice. Cocokyun lies more than 100 km offshore; any neighbour assigned to
# it would be invented, and it carries no ACLED record of any kind (checked in
# 02_acled.R), so it is removed from the analytical frame instead and kept only
# as a cartographic feature.
isl_dist <- lapply(no_nb, function(i) {
  d <- as.numeric(st_distance(ts_proj[i, ], ts_proj)); d[i] <- Inf
  j <- which.min(d)
  data.frame(TS_PCODE = ts_proj$TS_PCODE[i], TS = ts_proj$TS[i],
             nearest = ts_proj$TS[j], nearest_pcode = ts_proj$TS_PCODE[j],
             dist_km = d[j] / 1000)
}) |> bind_rows()

OFFSHORE_KM <- 50
excluded <- isl_dist$TS_PCODE[isl_dist$dist_km > OFFSHORE_KM]
linked   <- isl_dist[isl_dist$dist_km <= OFFSHORE_KM, ]

# The frame is sorted by TS_PCODE. sfdep's emerging_hotspot_analysis() labels
# its output rows in geometry order but computes them in sorted-location
# order, so the two must coincide for the labels to be right (see Section 6
# of the report, where this is demonstrated).
ts_an <- ts_proj |>
  filter(!TS_PCODE %in% excluded) |>
  arrange(TS_PCODE)

nb_an <- poly2nb(ts_an, queen = TRUE)
for (r in seq_len(nrow(linked))) {
  i <- match(linked$TS_PCODE[r], ts_an$TS_PCODE)
  j <- match(linked$nearest_pcode[r], ts_an$TS_PCODE)
  nb_an[[i]] <- sort(unique(c(setdiff(nb_an[[i]], 0L), j)))
  nb_an[[j]] <- sort(unique(c(setdiff(nb_an[[j]], 0L), i)))
}
attr(nb_an, "region.id") <- ts_an$TS_PCODE
stopifnot(all(card(nb_an) > 0), n.comp.nb(nb_an)$nc == 1L,
          is.symmetric.nb(nb_an))

# A second, structurally different neighbour definition is kept for the
# sensitivity analysis. k = 6 nearest neighbours by point-on-surface is close
# to the queen average of about 5.2 links, but every township gets exactly six,
# so large Shan and Kachin townships and small Yangon townships are treated
# alike. It is asymmetric, which row-standardised local statistics tolerate.
pos <- st_coordinates(st_point_on_surface(st_geometry(ts_an)))
nb_knn6 <- knn2nb(knearneigh(pos, k = 6), row.names = ts_an$TS_PCODE)

nb_report$islands_dist <- isl_dist
nb_report$excluded     <- excluded
nb_report$an_n         <- nrow(ts_an)
nb_report$an_avg_links <- mean(card(nb_an))
nb_report$an_max_links <- max(card(nb_an))
nb_report$knn6_sym     <- is.symmetric.nb(nb_knn6)

# ---- Dissolved frames for cartography ---------------------------------------
state  <- ts_proj |> group_by(ST, ST_PCODE) |> summarise(.groups = "drop")
nation <- st_union(ts_proj) |> st_sf(geometry = _)

# ---- Save --------------------------------------------------------------------
saveRDS(list(work_crs = WORK_CRS, crs_scores = crs_scores,
             geom_report = geom_report, nb_report = nb_report,
             islands = islands, area_geod_km2 = sum(area_geod) / 1e6),
        file.path(DER, "bnd_report.rds"))
saveRDS(ts_proj,  file.path(DER, "township.rds"))
saveRDS(ts_an,    file.path(DER, "township_an.rds"))
saveRDS(state,    file.path(DER, "state.rds"))
saveRDS(nation,   file.path(DER, "nation.rds"))
saveRDS(nb_queen, file.path(DER, "nb_queen.rds"))
saveRDS(nb_an,    file.path(DER, "nb_an.rds"))
saveRDS(nb_knn6,  file.path(DER, "nb_knn6.rds"))

cat("\n== CRS candidates (ranked by worst-township area error) ==\n")
print(crs_scores, row.names = FALSE, digits = 4)
cat("\nselected:", crs_scores$crs[1], "\n")
cat("\n== Geometry ==\n"); str(geom_report)
cat("\n== Queen contiguity ==\n"); str(nb_report)
if (nrow(islands)) { cat("\n== Townships with no queen neighbour ==\n"); print(isl_dist, row.names = FALSE) }
cat("\nexcluded from analysis:", excluded, "\n")
cat("analytical frame:", nrow(ts_an), "townships, mean", round(mean(card(nb_an)), 2), "links\n")
cat("\ngeodesic area of Myanmar:", format(round(sum(area_geod)/1e6), big.mark = ","), "km2\n")
