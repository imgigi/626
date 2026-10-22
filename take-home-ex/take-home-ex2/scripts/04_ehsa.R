# =============================================================================
# 04_ehsa.R -- space-time cube, Mann-Kendall and emerging hot spot analysis
# Take-home Exercise 2 | ISSS626
#
# Three stages.
#   1. Validation. sfdep::emerging_hotspot_analysis() is run on the quarterly
#      cube of all conflict events and its output is checked against an
#      independent computation of the same statistic (R/ehsa.R). Two defects are
#      measured rather than asserted: the row labels depend on the order of the
#      geometry, and the permutation null pools every period of the cube.
#   2. Primary analysis. For each of the five conflict variables, space-time
#      Gi* with one quarter of temporal lag, a period-conditional permutation
#      null, the Mann-Kendall test on each township's Gi* series, and sfdep's
#      own classification rules at p = 0.05.
#   3. Sensitivity. The time lag, the neighbour definition, the time unit, and
#      the p-value convention are each changed in turn.
# Run from take-home-ex/take-home-ex2/.
# =============================================================================

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(tidyr); library(spdep); library(sfdep)
})
source("R/utils.R")
source("R/ehsa.R")

DER  <- "data/derived"
SEED <- 1234
NSIM <- 999          # period-conditional permutations, primary analysis
NSIM_SFDEP <- 199    # sfdep's own default, used for the validation run
ALPHA <- 0.05

ts_an   <- readRDS(file.path(DER, "township_an.rds"))
ts_raw  <- readRDS(file.path(DER, "township.rds"))
nb_an   <- readRDS(file.path(DER, "nb_an.rds"))
nb_knn6 <- readRDS(file.path(DER, "nb_knn6.rds"))
panel_q <- readRDS(file.path(DER, "panel_q.rds"))
panel_m <- readRDS(file.path(DER, "panel_m.rds"))

VARS <- c(conflict = "All armed-conflict events", battles = "Battles",
          erv = "Explosions / remote violence",
          vac = "Violence against civilians", fatalities = "Reported fatalities")

locs <- ts_an$TS_PCODE
stopifnot(!is.unsorted(locs))     # the order sfdep's labels silently assume
nb_self <- include_self(nb_an)

# ---- The space-time cube --------------------------------------------------------
# sfdep's spacetime object binds the panel to the geometry and is checked for
# completeness: every township in every quarter, exactly once.
geo <- ts_an |> select(TS_PCODE)
geo$nb <- nb_self                    # assigned directly: dplyr >= 1.2 refuses
geo$wt <- st_weights(nb_self)        # an `nb` object inside mutate()

cube_q <- spacetime(as.data.frame(panel_q), geo,
                    .loc_col = "TS_PCODE", .time_col = "period")
stopifnot(is_spacetime_cube(cube_q))

cube_info <- list(n_locs = length(locs), n_times = length(unique(panel_q$period)),
                  n_bins = nrow(panel_q),
                  periods = sort(unique(panel_q$period)))

# ---- 1. Validation against sfdep --------------------------------------------------
X_conf <- panel_matrix(panel_q, "conflict", locs)
G_conf <- gi_spt(X_conf, nb_self, k = 1L)
mk_conf <- apply(G_conf, 1, function(g) Kendall::MannKendall(g)$tau)

ehsa_sf <- cache_rds("ehsa_sfdep_conflict", {
  set.seed(SEED)
  emerging_hotspot_analysis(cube_q, "conflict", k = 1, include_gi = TRUE,
                            nb_col = "nb", wt_col = "wt",
                            nsim = NSIM_SFDEP, threshold = ALPHA)
})
gi_sf <- attr(ehsa_sf, "gi_star")
gi_sf$TS_PCODE <- rep(locs, times = cube_info$n_times)
gi_sf$period   <- rep(cube_info$periods, each = length(locs))

val <- list()
# (a) The deterministic part, Gi* itself, is the same statistic.
val$max_gi_diff <- max(abs(as.vector(G_conf) - gi_sf$gi_star))
# (b) With the geometry sorted by location id, sfdep's tau is the township's own.
tau_sorted <- data.frame(location = ehsa_sf$location, tau_sf = ehsa_sf$tau) |>
  left_join(data.frame(location = locs, tau_me = mk_conf), by = "location")
val$max_tau_diff_sorted <- max(abs(tau_sorted$tau_sf - tau_sorted$tau_me))

# (c) The same cube with the geometry in MIMU's delivered order. Only tau is
# compared, and tau does not depend on the permutations, so a handful of them
# is enough here.
ord <- match(ts_raw$TS_PCODE[ts_raw$TS_PCODE %in% locs], locs)
inv <- integer(length(locs)); inv[ord] <- seq_along(ord)
nb_mimu <- lapply(seq_along(ord), function(p) sort(inv[nb_an[[ord[p]]]]))
class(nb_mimu) <- "nb"; attr(nb_mimu, "region.id") <- locs[ord]
nb_mimu <- include_self(nb_mimu)
geo_mimu <- ts_an[ord, ] |> select(TS_PCODE)
geo_mimu$nb <- nb_mimu
geo_mimu$wt <- st_weights(nb_mimu)
cube_mimu <- spacetime(as.data.frame(panel_q), geo_mimu,
                       .loc_col = "TS_PCODE", .time_col = "period")

ehsa_mimu <- cache_rds("ehsa_sfdep_conflict_mimu_order", {
  set.seed(SEED)
  emerging_hotspot_analysis(cube_mimu, "conflict", k = 1,
                            nb_col = "nb", wt_col = "wt",
                            nsim = 19, threshold = ALPHA)
})
tau_mimu <- data.frame(location = ehsa_mimu$location, tau_sf = ehsa_mimu$tau) |>
  left_join(data.frame(location = locs, tau_me = mk_conf), by = "location")
val$share_wrong_mimu  <- mean(abs(tau_mimu$tau_sf - tau_mimu$tau_me) > 1e-9)
val$cor_tau_mimu      <- cor(tau_mimu$tau_sf, tau_mimu$tau_me)
# Relabelling sfdep's rows in alphabetical order of location recovers the
# right answer exactly, which pins the cause on the labelling step.
tau_relab <- data.frame(location = sort(ehsa_mimu$location), tau_sf = ehsa_mimu$tau) |>
  left_join(data.frame(location = locs, tau_me = mk_conf), by = "location")
val$max_tau_diff_relabelled <- max(abs(tau_relab$tau_sf - tau_relab$tau_me))
val$class_changed_mimu <- sum(
  ehsa_mimu$classification[match(locs, ehsa_mimu$location)] !=
  ehsa_sf$classification[match(locs, ehsa_sf$location)])
val$tau_mimu <- tau_mimu

# (d) Feeding sfdep's own per-bin p-values through ehsa_classify() with sfdep's
# rules unchanged reproduces sfdep's classification exactly, so the classifier
# is not a source of difference either.
P_sf <- matrix(gi_sf$p_sim, nrow = length(locs))
cls_sf_check <- ehsa_classify(G_conf, P_sf, ALPHA, locs, corrected = FALSE)
val$classifier_identical <- identical(
  cls_sf_check$classification,
  sub("oscilating", "oscillating", ehsa_sf$classification[match(locs, ehsa_sf$location)]))

# ---- 2. Primary analysis: period-conditional null --------------------------------
# Only the expensive part, Gi* and its permutation p-values, is cached. The
# classification is cheap and is always recomputed, so a change to the class
# rules can never be hidden behind a stale cache.
run_ehsa <- function(panel, var, nb_s, k = 1L, nsim = NSIM, seed = SEED) {
  X <- panel_matrix(panel, var, locs)
  G <- gi_spt(X, nb_s, k = k)
  P <- gi_spt_perm(X, nb_s, k = k, nsim = nsim, seed = seed)
  list(X = X, G = G, p_fold = P$p_fold, p_two = P$p_two)
}
add_classes <- function(r) {
  r$cls      <- ehsa_classify(r$G, r$p_fold, ALPHA, locs)
  r$cls_two  <- ehsa_classify(r$G, r$p_two,  ALPHA, locs)
  r$cls_sfdep_rules <- ehsa_classify(r$G, r$p_fold, ALPHA, locs, corrected = FALSE)
  r
}

ehsa <- lapply(names(VARS), function(v) {
  add_classes(cache_rds(paste0("ehsa_pc_q_k1_", v), run_ehsa(panel_q, v, nb_self)))
}) |> setNames(names(VARS))

# Townships whose class changes when the one sfdep rule is corrected.
rule_fix <- lapply(names(VARS), function(v) {
  r <- ehsa[[v]]
  ch <- which(r$cls$classification != r$cls_sfdep_rules$classification)
  tibble(var = v, TS_PCODE = locs[ch],
         sfdep_rule = r$cls_sfdep_rules$classification[ch],
         corrected  = r$cls$classification[ch],
         n_hot  = rowSums(r$p_fold[ch, , drop = FALSE] <= ALPHA & r$G[ch, , drop = FALSE] > 0),
         n_cold = rowSums(r$p_fold[ch, , drop = FALSE] <= ALPHA & r$G[ch, , drop = FALSE] < 0))
}) |> bind_rows()

# Per-quarter counts of significant hot and cold bins under the two nulls: the
# figure that shows where sfdep's classes come from.
bins_by_null <- bind_rows(
  gi_sf |> group_by(period) |>
    summarise(hot = sum(p_sim <= ALPHA & gi_star > 0),
              cold = sum(p_sim <= ALPHA & gi_star < 0), .groups = "drop") |>
    mutate(null = "sfdep: pooled over all quarters"),
  tibble(period = rep(cube_info$periods, each = length(locs)),
         g = as.vector(ehsa$conflict$G), p = as.vector(ehsa$conflict$p_fold)) |>
    group_by(period) |>
    summarise(hot = sum(p <= ALPHA & g > 0), cold = sum(p <= ALPHA & g < 0),
              .groups = "drop") |>
    mutate(null = "Period-conditional")
) |>
  left_join(panel_q |> group_by(period) |> summarise(national = sum(conflict)),
            by = "period")

# The two classifications side by side, for all conflict events.
cls_compare <- tibble(
  TS_PCODE = locs,
  sfdep = ehsa_label(sub("oscilating", "oscillating",
                         ehsa_sf$classification[match(locs, ehsa_sf$location)])),
  pc    = ehsa_label(ehsa$conflict$cls$classification)
)

# sfdep's EHSA for the other four variables, so that the comparison of class
# counts covers every variable and not only the one examined in detail.
ehsa_sf_all <- lapply(names(VARS), function(v) {
  e <- if (v == "conflict") ehsa_sf else cache_rds(paste0("ehsa_sfdep_", v), {
    set.seed(SEED)
    emerging_hotspot_analysis(cube_q, v, k = 1, nb_col = "nb", wt_col = "wt",
                              nsim = NSIM_SFDEP, threshold = ALPHA)
  })
  tibble(var = v, TS_PCODE = e$location,
         classification = ehsa_label(sub("oscilating", "oscillating", e$classification)))
}) |> bind_rows()

# ---- How many townships would each class hold by chance? ---------------------------
# Every township is tested once per quarter, 19 times, and a single significant
# hot quarter with no significant cold one is enough for "sporadic hot spot".
# The class counts therefore need a reference. Here each township's complete
# 19-quarter history is moved to a randomly chosen location, which keeps every
# series and its serial dependence intact and destroys only the spatial
# arrangement, and the whole EHSA is re-run. The p-values inside each null cube
# use 99 permutations, which resolves the 0.05 threshold.
# Each null cube also records how many quarters every township spends as a
# significant hot spot, so that the count of significant quarters, which says
# more than the class label does, gets a reference distribution of its own.
NNULL <- 99
class_levels <- names(PAL_EHSA)
null_sim <- cache_rds("ehsa_null_conflict_v3", {
  set.seed(SEED)
  seeds <- sample.int(1e6, NNULL)
  lapply(seq_len(NNULL), function(r) {
    set.seed(seeds[r])
    Xn <- X_conf[sample.int(nrow(X_conf)), , drop = FALSE]
    Gn <- gi_spt(Xn, nb_self, k = 1L)
    Pn <- gi_spt_perm(Xn, nb_self, k = 1L, nsim = 99, seed = seeds[r])
    list(classes = as.vector(table(factor(
           ehsa_label(ehsa_classify(Gn, Pn$p_fold, ALPHA, locs)$classification),
           levels = class_levels))),
         n_hot  = tabulate(rowSums(Pn$p_fold <= ALPHA & Gn > 0) + 1L, nbins = ncol(Xn) + 1L),
         n_cold = tabulate(rowSums(Pn$p_fold <= ALPHA & Gn < 0) + 1L, nbins = ncol(Xn) + 1L))
  })
})
null_counts <- sapply(null_sim, `[[`, "classes")
rownames(null_counts) <- class_levels
# Pooled over all null cubes: the distribution of significant hot (and cold)
# quarters for a township whose neighbours are a random set.
null_nq <- tibble(
  k      = 0:cube_info$n_times,
  p_hot  = rowSums(sapply(null_sim, `[[`, "n_hot"))  / (NNULL * length(locs)),
  p_cold = rowSums(sapply(null_sim, `[[`, "n_cold")) / (NNULL * length(locs))
) |>
  mutate(surv_hot = rev(cumsum(rev(p_hot))), surv_cold = rev(cumsum(rev(p_cold))))

obs_counts <- as.vector(table(factor(ehsa_label(ehsa$conflict$cls$classification),
                                     levels = class_levels)))
null_tbl <- tibble(
  class    = factor(class_levels, levels = class_levels),
  observed = obs_counts,
  null_mean = rowMeans(null_counts),
  null_lo  = apply(null_counts, 1, quantile, 0.025),
  null_hi  = apply(null_counts, 1, quantile, 0.975),
  p_upper  = (rowSums(null_counts >= obs_counts) + 1) / (NNULL + 1)
) |>
  filter(observed > 0 | null_mean > 0.5)

# ---- 3. Sensitivity (all conflict events) ------------------------------------------
# One choice changes at a time; the time-lag variants and the neighbour variant
# use the quarterly cube, and the monthly variant keeps k = 1, i.e. one month.
nb_knn6_self <- include_self(nb_knn6)

sens_runs <- list(
  "k = 0 (no time lag)"  = cache_rds("ehsa_pc_q_k0_conflict", run_ehsa(panel_q, "conflict", nb_self, k = 0L)),
  "k = 2 (two quarters)" = cache_rds("ehsa_pc_q_k2_conflict", run_ehsa(panel_q, "conflict", nb_self, k = 2L)),
  "k = 6 nearest neighbours" = cache_rds("ehsa_pc_q_knn6_conflict", run_ehsa(panel_q, "conflict", nb_knn6_self)),
  "Monthly cube (57 months)" = cache_rds("ehsa_pc_m_k1_conflict", run_ehsa(panel_m, "conflict", nb_self))
) |> lapply(add_classes)

prim_cls <- ehsa_label(ehsa$conflict$cls$classification)
sens_tbl <- bind_rows(
  tibble(spec = "Primary (quarterly, k = 1, queen, folded p)", cls = prim_cls),
  tibble(spec = "Two-sided p ≤ 0.05", cls = ehsa_label(ehsa$conflict$cls_two$classification)),
  lapply(names(sens_runs), function(nm)
    tibble(spec = nm, cls = ehsa_label(sens_runs[[nm]]$cls$classification))) |> bind_rows()
) |>
  group_by(spec) |>
  mutate(same = cls == prim_cls) |>
  summarise(
    `Hot classes`  = sum(grepl("hot spot", cls)),
    `Cold classes` = sum(grepl("cold spot", cls)),
    `No pattern`   = sum(cls == "No pattern detected"),
    `Same class as primary` = sum(same),
    `Hot/cold side kept` = sum(sub(".* (hot|cold) spot$", "\\1", cls) ==
                               sub(".* (hot|cold) spot$", "\\1", prim_cls)),
    .groups = "drop")

# Serial dependence in the Gi* series. The Mann-Kendall test assumes independent
# observations; with k = 1 consecutive bins share a quarter of data, which
# induces positive autocorrelation and makes the trend test too liberal.
acf1 <- function(G) apply(G, 1, function(g) {
  if (sd(g) == 0) return(NA_real_)
  stats::acf(g, lag.max = 1, plot = FALSE)$acf[2]
})
serial <- tibble(
  spec = c("k = 1", "k = 0"),
  median_acf1 = c(median(acf1(ehsa$conflict$G), na.rm = TRUE),
                  median(acf1(sens_runs[["k = 0 (no time lag)"]]$G), na.rm = TRUE)),
  share_sig_mk = c(mean(ehsa$conflict$cls$p_value < ALPHA),
                   mean(sens_runs[["k = 0 (no time lag)"]]$cls$p_value < ALPHA))
)

# ---- Save ----------------------------------------------------------------------------
slim <- function(r) r[c("G", "p_fold", "p_two", "cls", "cls_two")]
saveRDS(list(
  vars = VARS, alpha = ALPHA, nsim = NSIM, nsim_sfdep = NSIM_SFDEP, seed = SEED,
  locs = locs, cube_info = cube_info,
  validation = val, bins_by_null = bins_by_null, cls_compare = cls_compare,
  null_classes = null_tbl, null_counts = null_counts, nnull = NNULL,
  rule_fix = rule_fix,
  null_nq = null_nq,
  sfdep_all = ehsa_sf_all,
  results = lapply(ehsa, slim),
  sens = sens_tbl, sens_monthly = slim(sens_runs[["Monthly cube (57 months)"]]),
  sens_k0 = slim(sens_runs[["k = 0 (no time lag)"]]),
  serial = serial
), file.path(DER, "ehsa.rds"))

cat("\n== Validation ==\n"); str(val[setdiff(names(val), "tau_mimu")])
cat("\n== Significant bins per quarter, by null ==\n")
print(bins_by_null |> pivot_wider(id_cols = c(period, national), names_from = null,
                                  values_from = c(hot, cold)), n = 30)
cat("\n== Classes, primary ==\n")
for (v in names(VARS)) { cat("--", v, "\n"); print(table(ehsa_label(ehsa[[v]]$cls$classification))[table(ehsa_label(ehsa[[v]]$cls$classification)) > 0]) }
cat("\n== sfdep vs period-conditional (conflict) ==\n")
print(table(sfdep = cls_compare$sfdep, pc = cls_compare$pc) |> as.data.frame() |> filter(Freq > 0) |> arrange(desc(Freq)))
cat("\n== Corrected intensifying-hot-spot rule: classes changed ==\n"); print(rule_fix)
cat("\n== Class counts against spatial randomness ==\n"); print(null_tbl, width = 120)
cat("\n== Significant quarters under the null (survival) ==\n"); print(null_nq, n = 25)
cat("\n== Sensitivity ==\n"); print(sens_tbl, width = 160)
cat("\n== Serial dependence ==\n"); print(serial)
