# =============================================================================
# 03_lmsa.R -- local measures of spatial association on township totals
# Take-home Exercise 2 | ISSS626
#
# Computes, for each of the five armed-conflict variables summed over the whole
# study period, a global Moran's I as the gate, then local Moran's I and
# Getis-Ord Gi* by conditional permutation. A year-by-year local Moran's I of
# all conflict events follows, and then the sensitivity runs that decide how
# much of the cluster map belongs to the data and how much to the analytical
# choices: the neighbour definition, the transformation of the counts, and the
# multiple-testing correction.
# Run from take-home-ex/take-home-ex2/.
# =============================================================================

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(tidyr); library(spdep); library(sfdep)
})
source("R/utils.R")

DER  <- "data/derived"
NSIM <- 999
SEED <- 1234

ts      <- readRDS(file.path(DER, "ts_tot.rds"))
nb_an   <- readRDS(file.path(DER, "nb_an.rds"))
nb_knn6 <- readRDS(file.path(DER, "nb_knn6.rds"))
panel_y <- readRDS(file.path(DER, "panel_y.rds"))

VARS <- c(conflict = "All armed-conflict events", battles = "Battles",
          erv = "Explosions / remote violence",
          vac = "Violence against civilians", fatalities = "Reported fatalities")

# ---- Spatial weights -----------------------------------------------------------
# Row standardisation ("W") makes the spatial lag the average of the neighbours,
# which is what the LISA quadrants are read against. Gi* needs the township
# itself inside its own neighbourhood, so a self-inclusive copy is kept for it.
lw_q    <- nb2listw(nb_an,   style = "W")
lw_k    <- nb2listw(nb_knn6, style = "W")
nb_self <- include_self(nb_an)
wt_self <- st_weights(nb_self)

# ---- Global Moran's I: is there spatial structure to localise? ----------------
glob <- lapply(names(VARS), function(v) {
  set.seed(SEED)
  m <- moran.mc(ts[[v]], lw_q, nsim = NSIM, alternative = "greater")
  data.frame(var = v, I = unname(m$statistic), p = m$p.value)
}) |> bind_rows()

# ---- The local statistics -------------------------------------------------------
# Both statistics are tested by conditional permutation: the township's own value
# is held fixed and its neighbours are redrawn from the other 328 townships.
# That avoids the normal approximation, which is a poor reference for counts as
# skewed as these. The p-value carried forward is spdep's two-sided simulated
# p-value ("Pr(z != E(Ii)) Sim"), not the folded one, which is one-sided and
# would halve every p-value.
local_moran_tbl <- function(x, lw, nsim = NSIM, seed = SEED) {
  set.seed(seed)
  lm <- localmoran_perm(x, lw, nsim = nsim, alternative = "two.sided")
  lag_x <- lag.listw(lw, x)
  tibble(ii = lm[, "Ii"], z_ii = lm[, "Z.Ii"],
         p_ii_sim = lm[, "Pr(z != E(Ii)) Sim"],
         lag = lag_x,
         quad = lisa_quadrant(x, lag_x, lm[, "Pr(z != E(Ii)) Sim"]),
         quad_fdr = lisa_quadrant(x, lag_x, lm[, "Pr(z != E(Ii)) Sim"],
                                  method = "fdr"))
}

local_gi_tbl <- function(x, nb, wt, nsim = NSIM, seed = SEED) {
  set.seed(seed)
  g <- local_gstar_perm(x, nb, wt, nsim = nsim)
  tibble(gi = g$gi_star, p_gi_sim = g$p_sim,
         gi_class = gi_class(g$gi_star, g$p_sim),
         gi_class_fdr = gi_class(g$gi_star, g$p_sim, method = "fdr"))
}

lisa <- lapply(names(VARS), function(v) {
  bind_cols(
    tibble(TS_PCODE = ts$TS_PCODE, var = v, x = ts[[v]]),
    local_moran_tbl(ts[[v]], lw_q),
    local_gi_tbl(ts[[v]], nb_self, wt_self)
  )
}) |> bind_rows()

# ---- Year by year: does the cluster map move? -----------------------------------
# Moran's I is unchanged by rescaling, so the nine-month 2025 needs no per-month
# adjustment here; it would matter only for comparing levels across years.
lisa_year <- lapply(2021:2025, function(y) {
  d <- panel_y |> filter(year == y) |> arrange(TS_PCODE)
  stopifnot(identical(d$TS_PCODE, ts$TS_PCODE))
  bind_cols(tibble(TS_PCODE = d$TS_PCODE, year = y, x = d$conflict),
            local_moran_tbl(d$conflict, lw_q))
}) |> bind_rows()

glob_year <- lapply(2021:2025, function(y) {
  d <- panel_y |> filter(year == y) |> arrange(TS_PCODE)
  set.seed(SEED)
  m <- moran.mc(d$conflict, lw_q, nsim = NSIM, alternative = "greater")
  data.frame(year = y, I = unname(m$statistic), p = m$p.value)
}) |> bind_rows()

# ---- Sensitivity ----------------------------------------------------------------
# Each alternative changes exactly one choice relative to the primary run
# (queen contiguity with the two island links, raw counts, unadjusted p < 0.05)
# and reports how many townships keep the same class. The density variant
# divides by township area, which mainly tests whether the large Shan and
# Kachin townships are hot only because they are large.
sens_spec <- list(
  "Primary (queen, raw counts)"   = list(x = function(d) d$conflict,            lw = lw_q),
  "k = 6 nearest neighbours"      = list(x = function(d) d$conflict,            lw = lw_k),
  "log(1 + count)"                = list(x = function(d) log1p(d$conflict),     lw = lw_q),
  "Events per 1,000 km²"          = list(x = function(d) 1000 * d$conflict / d$area_km2, lw = lw_q)
)

sens <- lapply(names(sens_spec), function(nm) {
  s <- sens_spec[[nm]]
  out <- local_moran_tbl(s$x(ts), s$lw)
  tibble(TS_PCODE = ts$TS_PCODE, spec = nm, quad = out$quad, quad_fdr = out$quad_fdr)
}) |> bind_rows()

prim <- sens |> filter(spec == names(sens_spec)[1]) |> select(TS_PCODE, quad_p = quad)
sens_summary <- sens |>
  left_join(prim, by = "TS_PCODE") |>
  group_by(spec) |>
  summarise(`High-High` = sum(quad == "High-High"), `Low-Low` = sum(quad == "Low-Low"),
            `High-Low` = sum(quad == "High-Low"), `Low-High` = sum(quad == "Low-High"),
            `Survive FDR` = sum(quad_fdr != "Insignificant"),
            `Same class as primary` = sum(quad == quad_p),
            `HH kept` = sum(quad == "High-High" & quad_p == "High-High"),
            .groups = "drop") |>
  mutate(spec = factor(spec, levels = names(sens_spec))) |>
  arrange(spec)

# ---- Save -----------------------------------------------------------------------
saveRDS(list(vars = VARS, nsim = NSIM, seed = SEED, glob = glob, lisa = lisa,
             lisa_year = lisa_year, glob_year = glob_year, sens = sens,
             sens_summary = sens_summary),
        file.path(DER, "lmsa.rds"))

cat("\n== Global Moran's I ==\n"); print(glob, digits = 3)
cat("\n== LISA quadrants (p < 0.05) ==\n")
print(lisa |> count(var, quad) |> pivot_wider(names_from = quad, values_from = n, values_fill = 0))
cat("\n== Gi* classes (p < 0.05) ==\n")
print(lisa |> count(var, gi_class) |> pivot_wider(names_from = gi_class, values_from = n, values_fill = 0))
cat("\n== Surviving FDR ==\n")
print(lisa |> group_by(var) |> summarise(lisa_fdr = sum(quad_fdr != "Insignificant"),
                                          gi_fdr = sum(gi_class_fdr != "Insignificant")))
cat("\n== By year ==\n"); print(glob_year, digits = 3)
print(lisa_year |> count(year, quad) |> pivot_wider(names_from = quad, values_from = n, values_fill = 0))
cat("\n== Sensitivity ==\n"); print(sens_summary)
