# =============================================================================
# 04_covariates.R
# Take-home Exercise 1 | ISSS626
#
# PURPOSE
#   Build the spatial covariate surfaces used to (a) explain first-order
#   intensity variation (rhohat) and (b) fit the inhomogeneous Poisson and
#   cluster point-process models in the value-added analysis.
#
# COVARIATES AND WHY EACH ONE
#   d_canal  distance to the nearest canal / drain / ditch. Pulang Pisau
#            contains a large share of the abandoned Mega Rice Project
#            drainage network. Drainage lowers the peat water table, which is
#            the mechanism most consistently linked in the literature to
#            recurrent peat fire. This is the covariate of primary interest.
#   d_river  distance to the nearest natural river / stream. Separates the
#            *engineered* drainage signal from proximity to water generally.
#   d_road   distance to the nearest road. Proxy for human access, hence
#            ignition opportunity, and also for suppression access.
#   d_settle distance to the nearest mapped settlement. Proxy for agricultural
#            land preparation, the dominant proximate ignition cause.
#
#   All four are *proxies*. None measures peat depth, water-table depth, land
#   cover or land tenure, which the interpretation must acknowledge.
#
# RESOLUTION
#   250 m. Chosen to sit below the 375 m sensor pixel (so the covariate is not
#   the coarser of the two layers) while keeping the raster small enough for
#   repeated ppm/kppm fitting.
# =============================================================================

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(terra); library(spatstat.geom)
})
sf_use_s2(FALSE)

DER <- "data/derived"
P   <- readRDS(file.path(DER, "params.rds"))
RES <- 250   # metres

bnd <- readRDS(file.path(DER, "bnd_proj.rds"))
win <- readRDS(file.path(DER, "win.rds"))

# ---- Template raster covering the study window -----------------------------
bv   <- vect(bnd)
tmpl <- rast(ext(bv), resolution = RES, crs = crs(bv))

# ---- Helper: distance-to-nearest-feature surface ---------------------------
# Lines/points are first rasterised onto the template, then terra::distance()
# returns the Euclidean distance from every cell to the nearest occupied cell.
# This is far cheaper than an exact point-to-geometry distance against ~30,000
# OSM linestrings, and the discretisation error is bounded by the cell size.
dist_surface <- function(g, tmpl) {
  if (is.null(g) || nrow(g) == 0) return(NULL)
  g <- st_transform(st_geometry(g), P$crs_proj)
  r <- rasterize(vect(g), tmpl, field = 1, touches = TRUE)
  terra::distance(r)
}

read_osm <- function(f) {
  p <- file.path("data/raw/osm", f)
  if (!file.exists(p)) return(NULL)
  x <- readRDS(p)
  if (is.null(x) || nrow(x) == 0) return(NULL)
  st_make_valid(x)
}

canal  <- read_osm("osm_canal_lines.rds")
river  <- read_osm("osm_river_lines.rds")
road   <- read_osm("osm_highway_lines.rds")
settle <- read_osm("osm_place_points.rds")

message("OSM features: canal=", ifelse(is.null(canal), 0, nrow(canal)),
        " river=",  ifelse(is.null(river), 0, nrow(river)),
        " road=",   ifelse(is.null(road), 0, nrow(road)),
        " settle=", ifelse(is.null(settle), 0, nrow(settle)))

rs <- list(
  d_canal  = dist_surface(canal,  tmpl),
  d_river  = dist_surface(river,  tmpl),
  d_road   = dist_surface(road,   tmpl),
  d_settle = dist_surface(settle, tmpl)
)
rs <- rs[!vapply(rs, is.null, logical(1))]

# ---- Mask to the study window, convert to spatstat images ------------------
# ppm() requires covariates as `im` objects defined on (at least) the window.
to_im <- function(r, nm) {
  r <- mask(crop(r, bv), bv)
  names(r) <- nm
  m  <- t(as.matrix(r, wide = TRUE))          # terra rows run north->south
  xx <- xFromCol(r, seq_len(ncol(r)))
  yy <- yFromRow(r, seq_len(nrow(r)))
  im(t(m)[rev(seq_len(nrow(r))), ], xcol = xx, yrow = rev(yy), unitname = c("metre", "metres"))
}

covars <- lapply(names(rs), function(nm) to_im(rs[[nm]], nm))
names(covars) <- names(rs)

# Sanity check: every covariate must cover the window and be finite inside it.
for (nm in names(covars)) {
  v <- covars[[nm]][win, drop = TRUE]
  message(sprintf("%-9s n=%d  min=%7.0f  median=%7.0f  max=%8.0f  NA=%d",
                  nm, length(v), min(v, na.rm = TRUE), median(v, na.rm = TRUE),
                  max(v, na.rm = TRUE), sum(is.na(v))))
}

saveRDS(covars, file.path(DER, "covars.rds"))
terra::writeRaster(rast(rs), file.path(DER, "covars.tif"), overwrite = TRUE)

# Keep the vector layers too, clipped to the regency, for cartography.
clip <- function(x) if (is.null(x)) NULL else
  suppressWarnings(st_intersection(st_transform(x, P$crs_proj), bnd))
saveRDS(list(canal = clip(canal), river = clip(river),
             road = clip(road),  settle = clip(settle)),
        file.path(DER, "osm_clipped.rds"))

message("\nSaved covars.rds (", length(covars), " surfaces at ", RES, " m)")
