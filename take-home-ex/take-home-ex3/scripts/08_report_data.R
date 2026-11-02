# =============================================================================
# 08_report_data.R -- what the report reads from the calibration
# Take-home Exercise 3 | ISSS626
#
# model.rds keeps every fitted specification with its matrices and runs to tens of
# megabytes. The report needs only the parameters, scores and tables, which go in a
# small file that can be committed, so the report renders without the slow scripts.
# Run from take-home-ex/take-home-ex3/.
# =============================================================================

suppressPackageStartupMessages(library(dplyr))
DER <- "data/derived"
MD <- readRDS(file.path(DER, "model.rds"))
q_of <- function(f) f$par

saveRDS(list(
  DMAX = MD$DMAX, X_VARS = MD$X_VARS,
  ladder_tbl = MD$ladder_tbl,
  fold_tbl = MD$fold_tbl, auc_low = MD$auc_low, low_cut = MD$low_cut, hold_bias = MD$hold_bias,
  prof_dmax = MD$prof_dmax, prof_lambda = MD$prof_lambda, prof_a0 = MD$prof_a0,
  perm_null = MD$perm_null, perm_p = MD$perm_p, obs_r2 = MD$obs_r2,
  coef_H = MD$coef_H, coef_W = MD$coef_W, coef_O = MD$coef_O,
  cal_tbl = MD$cal_tbl, residuals_tbl = MD$residuals_tbl, moran = MD$moran, mci_check = MD$mci_check,
  par_H = MD$sel$H$par, par_W = MD$sel$W$par, par_O = MD$open$H$par,
  lambda_m7 = c(H = MD$ladder$m7$H$par[1], W = MD$ladder$m7$W$par[1]),
  lolo_H = MD$sel$lolo_H, lolo_W = MD$sel$lolo_W, drm_coef = MD$drm_coef
), file.path(DER, "report.rds"))

R <- readRDS(file.path(DER, "report.rds"))
SC <- readRDS(file.path(DER, "scenario.rds")); MC <- readRDS(file.path(DER, "mci.rds"))
AM <- readRDS(file.path(DER, "audit_markets.rds")); AS <- readRDS(file.path(DER, "audit_stations.rds")); AR <- readRDS(file.path(DER, "audit_ridership.rds"))
print(SC$totals, width = 200); print(SC$totals_ci); print(SC$tengah); print(SC$add26_total)
print(R$cal_tbl |> arrange(desc(km_calibrated)) |> head(6)); print(R$cal_tbl |> summarise(s = median(km_structural), c = median(km_calibrated)))
print(R$residuals_tbl |> arrange(resid) |> select(station_name, daily_in, E, resid) |> head(8)); print(R$residuals_tbl |> arrange(desc(resid)) |> select(station_name, daily_in, E, resid) |> head(10))
print(R$moran); print(R$auc_low); print(R$low_cut); print(R$hold_bias); print(R$perm_p)
print(MC$links |> filter(a_name %in% c("Jurong East", "Clementi") & b_name %in% c("Jurong East", "Clementi")) |> select(line, a_name, b_name, base, after, pct))
print(MC$links |> summarise(max_base = max(base), who = paste(a_name[which.max(base)], b_name[which.max(base)]), line = line[which.max(base)]))
print(MC$links |> filter(line %in% c("JS", "JE", "JW", "CR")) |> arrange(desc(after)) |> select(line, a_name, b_name, after) |> head(6))
print(MC$line_summary); print(MC$r2_wls); print(MC$pkm); print(sum(MC$xfer$base)); print(sum(MC$xfer$after)); print(MC$in_out_cor); print(MC$balance_gap); print(MC$new_trip); print(MC$fit_arrivals); print(MC$fit_lolo)
print(AM$priv_by_tier); print(AM$osm_b_n); print(AM$osm_b_on_hdb); print(AS$layer_nonstation); print(AS$poi_by_group); print(AS$dorm_n); print(AS$work_n)
print(AR$daily_total); print(AR$wd_total); print(AR$we_total); print(AR$stab_aug_sep_cor); print(AR$stab_aug_sep_mape); print(AR$od_pairs); print(AR$od_pairs_possible)
d <- readRDS(file.path(DER, "dest_existing.rds")); print(sum(d$lrt)); print(d |> filter(lrt == 0) |> summarise(med = median(daily_in))); print(d |> filter(lrt == 1) |> summarise(med = median(daily_in)))
print(cor(log(d$daily_in), log(d$res_800 + 1))); print(d |> arrange(desc(daily_in)) |> select(station_name, daily_in, phi) |> head(5)); print(d |> filter(station_name %in% c("Hume", "Buangkok", "Woodleigh", "Keppel", "Teck Lee", "Cantonment", "Prince Edward Road")) |> select(station_name, daily_in, phi))
print(rank(-d$daily_in)[d$station_name %in% c("Buangkok", "Woodleigh", "Hume")])
message("ALLDONE")
