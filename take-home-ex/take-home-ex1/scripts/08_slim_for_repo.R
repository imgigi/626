# =============================================================================
# 08_slim_for_repo.R
# Take-home Exercise 1 | ISSS626
#
# PURPOSE
#   The analysis scripts save their full working objects, which is right for
#   analysis but wrong for a repository: `va_results.rds` comes out at ~121 MB,
#   above GitHub's 100 MB per-file hard limit, almost all of it Monte-Carlo
#   curves stored by `savefuns = TRUE` and quadrature schemes carried inside
#   fitted `ppm` objects.
#
#   This script rewrites the saved results to contain exactly what the report
#   reads and nothing else. It is lossless for every number and every figure in
#   the report; what it discards is only recomputable scaffolding.
#
#   Run AFTER 06 and 07, and BEFORE committing. Idempotent: running it twice is
#   harmless, because it only ever removes components that the report does not
#   read.
# =============================================================================

suppressPackageStartupMessages({library(sf); library(dplyr)})
DER <- "data/derived"

mb <- function(p) round(file.size(p) / 1048576, 1)

# ---- 1. Drop the simulated function curves from every envelope -------------
# An envelope built with savefuns = TRUE keeps all 199 (or 398) simulated
# curves in attr(, "simfuns"). Those were needed once, to compute the DCLF and
# MAD statistics; the test results are stored separately, so the curves are
# dead weight. Plotting needs only r / obs / mmean / lo / hi.
strip_env <- function(x) {
  if (inherits(x, "envelope") || inherits(x, "fv")) {
    attr(x, "simfuns") <- NULL
    attr(x, "simpatterns") <- NULL
  }
  x
}
# Recurse only into BARE lists, and rebuild with `x[] <-` so that attributes
# survive. An earlier version used `lapply()` on anything list-like, which
# silently stripped the class from every htest and every data frame inside one
# (an htest IS a list, and so is a data frame) -- the test objects came back as
# unclassed lists and the report failed on formatC() of a list.
strip_recursive <- function(x) {
  if (inherits(x, c("envelope", "fv"))) return(strip_env(x))
  if (is.list(x) && is.null(attr(x, "class"))) {
    x[] <- lapply(x, strip_recursive)
    return(x)
  }
  x
}

# ---- 2. Value-added results -------------------------------------------------
f_va <- file.path(DER, "va_results.rds")
before <- mb(f_va)
VA <- readRDS(f_va)

# The report read only two numbers out of the fitted model objects --
# length(coef(M1)) and psib(M3) -- yet those two objects accounted for ~74 MB
# of the 121 MB, because a ppm carries its quadrature scheme and a kppm carries
# the fitted covariance structure. Both numbers are therefore computed here,
# once, and stored as scalars; the report reads VA$n_par_M1 and VA$psib.
VA$n_par_M1 <- length(coef(VA$M1))
VA$psib <- tryCatch(unname(spatstat.model::psib(VA$M3)),
                    error = function(e) NA_real_)
# mu is a per-quadrature-point vector but the report only takes its mean.
VA$mu <- mean(VA$mu)

# Components the report actually reads (verified by grepping VA$ in the qmd
# sources).
keep_va <- c("n_par_M1", "psib", "aic_tab", "clust_tab", "clustpar", "eff",
             "gof_M1", "gof_M1_test", "gof_M2", "gof_M2_test",
             "gof_M3", "gof_M3_test", "mu", "parres", "res_sm",
             "trend_M1", "trend_M2")
VA <- VA[intersect(keep_va, names(VA))]   # drops fits, M1, M3, coefs, nsim
VA <- lapply(VA, strip_recursive)
saveRDS(VA, f_va, compress = "xz")
message(sprintf("va_results.rds : %6.1f MB -> %6.1f MB", before, mb(f_va)))

# ---- 3. Second-order results ------------------------------------------------
f_so <- file.path(DER, "so_results.rds")
before <- mb(f_so)
SO <- readRDS(f_so)
SO <- lapply(SO, strip_recursive)
saveRDS(SO, f_so, compress = "xz")
message(sprintf("so_results.rds : %6.1f MB -> %6.1f MB", before, mb(f_so)))

# ---- 4. First-order KDE / STKDE objects ------------------------------------
for (f in c("fo_kde.rds", "fo_stkde.rds", "fo_rhohat.rds")) {
  p <- file.path(DER, f)
  if (!file.exists(p)) next
  before <- mb(p)
  saveRDS(readRDS(p), p, compress = "xz")
  message(sprintf("%-15s: %6.1f MB -> %6.1f MB", f, before, mb(p)))
}

# ---- 5. The geoBoundaries cross-check needs one polygon, not 519 -----------
# Section 2 compares the study-window area against the independently produced
# geoBoundaries ADM2 polygon. Only Pulang Pisau is needed, so the national
# 519-feature object is reduced to that single row.
f_kal <- file.path(DER, "kal_adm2.rds")
if (file.exists(f_kal)) {
  before <- mb(f_kal)
  pp_gb <- readRDS(f_kal) |> filter(shapeName == "Pulang Pisau")
  stopifnot(nrow(pp_gb) == 1)
  saveRDS(pp_gb, f_kal, compress = "xz")
  message(sprintf("kal_adm2.rds   : %6.1f MB -> %6.1f MB (1 of 519 features kept)",
                  before, mb(f_kal)))
}

# ---- 6. Objects no part of the report reads --------------------------------
for (f in c("kal_reg.rds", "covars.tif", "fire_pxday_sf.rds.bak")) {
  p <- file.path(DER, f)
  if (file.exists(p)) { unlink(p); message("removed unused: ", f) }
}

message("\nderived total: ",
        round(sum(file.size(list.files(DER, full.names = TRUE))) / 1048576, 1),
        " MB")
