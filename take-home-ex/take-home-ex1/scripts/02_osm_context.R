# =============================================================================
# 02_osm_context.R
# Take-home Exercise 1 | ISSS626
#
# PURPOSE
#   Download OpenStreetMap contextual layers for the Pulang Pisau study area.
#   These are NOT used to define the point pattern. They serve two purposes:
#     (a) cartographic reference, so that density surfaces can be read against
#         the canal/road/settlement geography rather than against a blank map;
#     (b) spatial covariates for the inhomogeneous point-process model in the
#         value-added analysis (distance to canal, road, settlement).
#
#   Rationale for these particular layers: Pulang Pisau contains a large part
#   of the abandoned Mega Rice Project drainage network. Drainage lowers the
#   peat water table and is the mechanism most often cited for the recurrence
#   of peat fires, while roads and settlements proxy human access and therefore
#   ignition opportunity.
#
# SOURCE  OpenStreetMap contributors, via the Overpass API (osmdata).
# LICENCE ODbL 1.0.
# =============================================================================

suppressPackageStartupMessages({library(sf); library(osmdata); library(dplyr)})

BBOX <- c(113.43, -3.62, 114.50, -1.39)   # xmin, ymin, xmax, ymax (WGS84)
OUT  <- "data/raw/osm"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

grab <- function(key, value = NULL, geom = "lines", tag = key) {
  dest <- file.path(OUT, paste0("osm_", tag, "_", geom, ".rds"))
  if (file.exists(dest)) { message("cache hit: ", dest); return(invisible(dest)) }
  q <- opq(bbox = BBOX, timeout = 600)
  q <- if (is.null(value)) add_osm_feature(q, key = key)
       else add_osm_feature(q, key = key, value = value)
  for (a in 1:3) {
    res <- try(osmdata_sf(q), silent = TRUE)
    if (!inherits(res, "try-error")) {
      obj <- switch(geom,
                    lines    = res$osm_lines,
                    points   = res$osm_points,
                    polygons = res$osm_polygons,
                    multipolygons = res$osm_multipolygons)
      saveRDS(obj, dest)
      message(sprintf("%-28s %s features -> %s", tag,
                      ifelse(is.null(obj), 0, nrow(obj)), basename(dest)))
      return(invisible(dest))
    }
    message("  retry ", a, " for ", tag); Sys.sleep(20)
  }
  message("FAILED: ", tag); invisible(NA)
}

# Drainage network: canals/drains/ditches are the Mega Rice Project legacy;
# rivers/streams are the natural network.
grab("waterway", c("canal", "drain", "ditch"), "lines", "canal")
grab("waterway", c("river", "stream"),         "lines", "river")
# Road network as an access/ignition-opportunity proxy.
grab("highway",  NULL,                          "lines", "highway")
# Settlements.
grab("place", c("city","town","village","hamlet","suburb"), "points", "place")
# Built-up / plantation land use, for interpretation only.
grab("landuse", c("farmland","plantation","orchard","forest","residential"),
     "multipolygons", "landuse")

message("\nDone. Files in ", OUT)
