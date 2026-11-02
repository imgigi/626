# =============================================================================
# 03_ridership.R -- observed ridership at each physical station
# Take-home Exercise 3 | ISSS626
#
# Reads the DataMall passenger volumes for July to September 2026, merges the
# several codes a station can appear under, converts monthly totals to average
# daily volumes, and reconciles the OD file against the station file.
# Run from take-home-ex/take-home-ex3/.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr); library(purrr)
})

RAW <- "data/raw/lta"; DER <- "data/derived"
audit <- list()
MONTHS <- c("202607", "202608", "202609")

read_zip <- function(prefix, inner, m, types) {
  read_csv(unz(file.path(RAW, sprintf("%s_%s.zip", prefix, m)), sprintf("%s_%s.csv", inner, m)),
           show_col_types = FALSE, col_types = types)
}
node <- map_dfr(MONTHS, \(m) read_zip("PVTrain", "transport_node_train", m, "cciccdd"))
od   <- map_dfr(MONTHS, \(m) read_zip("PVODTrain", "origin_destination_train", m, "cciccci") |>
                  mutate(TOTAL_TRIPS = as.numeric(TOTAL_TRIPS)))
audit$node_rows <- nrow(node)
audit$od_rows   <- nrow(od)

# ---- Codes to physical stations ---------------------------------------------
# A node code joins the station codes of an interchange ("NS24/NE6/CC1"). Three
# stations appear under two node codes. Circle Line Stage 6 renumbered Marina Bay
# and Bayfront (CE2 to CC33, CE1 to CC34) when it opened on 12 July 2026, so both
# spellings are in the July file; and Stevens is reported as two separate nodes for
# its two lines, DT10 and TE11. The renumbered codes are not in the March
# coordinate list and are added here.
code_map <- readRDS(file.path(DER, "station_codes.rds")) |>
  select(code = station_code, station_id, station_name)
extra <- code_map |>
  filter(code %in% c("CE2", "CE1")) |>
  mutate(code = recode(code, CE2 = "CC33", CE1 = "CC34"))
code_map <- bind_rows(code_map, extra)

raw_codes <- tibble(raw = unique(c(node$PT_CODE, od$ORIGIN_PT_CODE, od$DESTINATION_PT_CODE))) |>
  mutate(code = str_split(raw, "[^A-Z0-9]+")) |>
  unnest(code)
lookup <- raw_codes |>
  left_join(code_map, by = "code") |>
  group_by(raw) |>
  summarise(n_st = n_distinct(station_id), unmatched = any(is.na(station_id)),
            station_id = first(station_id), station_name = first(station_name), .groups = "drop")
stopifnot(all(lookup$n_st == 1), !any(lookup$unmatched))
audit$raw_codes <- nrow(lookup)
audit$stations  <- n_distinct(lookup$station_id)
audit$multi_node <- lookup |>
  group_by(station_name) |>
  summarise(codes = paste(raw, collapse = " + "), n = n(), .groups = "drop") |>
  filter(n > 1)

# ---- Calendar -----------------------------------------------------------------
# Monthly totals are not comparable, because the months have different numbers of
# each kind of day. August 2026 has 20 weekdays against July's 23: National Day
# falls on a Sunday, so Monday 10 August is a public holiday.
calendar <- tibble(date = seq(as.Date("2026-07-01"), as.Date("2026-09-30"), by = "day")) |>
  mutate(YEAR_MONTH = format(date, "%Y-%m"),
         dow = as.integer(format(date, "%u")),
         DAY_TYPE = ifelse(dow >= 6 | date == as.Date("2026-08-10"), "WEEKENDS/HOLIDAY", "WEEKDAY")) |>
  count(YEAR_MONTH, DAY_TYPE, name = "days")
audit$calendar <- calendar

audit$month_check <- node |>
  group_by(YEAR_MONTH, DAY_TYPE) |>
  summarise(tap_in = sum(TOTAL_TAP_IN_VOLUME), .groups = "drop") |>
  left_join(calendar, by = c("YEAR_MONTH", "DAY_TYPE")) |>
  mutate(per_day = tap_in / days)

# ---- Reference period ---------------------------------------------------------
# July is left out of the calibration: Circle Line Stage 6 opened on 12 July, so its
# stations have a partial month, and the rest of the network was adjusting to a
# closed circle. August and September describe one stable network.
REF <- c("2026-08", "2026-09")
ref_days <- calendar |>
  filter(YEAR_MONTH %in% REF) |>
  group_by(DAY_TYPE) |>
  summarise(days = sum(days), .groups = "drop")
n_wd <- ref_days$days[ref_days$DAY_TYPE == "WEEKDAY"]
n_we <- ref_days$days[ref_days$DAY_TYPE != "WEEKDAY"]
audit$ref_days <- ref_days

per <- function(h) case_when(h %in% 6:8 ~ "am", h %in% 17:19 ~ "pm", TRUE ~ "off")

rid <- node |>
  filter(YEAR_MONTH %in% REF) |>
  left_join(select(lookup, raw, station_id), by = c("PT_CODE" = "raw")) |>
  mutate(period = per(TIME_PER_HOUR)) |>
  group_by(station_id) |>
  summarise(
    # Average day across the reference period, weighting weekdays and weekends by how
    # many of each a typical week holds, so that the calendar does not matter.
    wd_in  = sum(TOTAL_TAP_IN_VOLUME[DAY_TYPE == "WEEKDAY"]) / n_wd,
    we_in  = sum(TOTAL_TAP_IN_VOLUME[DAY_TYPE != "WEEKDAY"]) / n_we,
    wd_out = sum(TOTAL_TAP_OUT_VOLUME[DAY_TYPE == "WEEKDAY"]) / n_wd,
    we_out = sum(TOTAL_TAP_OUT_VOLUME[DAY_TYPE != "WEEKDAY"]) / n_we,
    am_in  = sum(TOTAL_TAP_IN_VOLUME[DAY_TYPE == "WEEKDAY" & period == "am"]) / n_wd,
    pm_in  = sum(TOTAL_TAP_IN_VOLUME[DAY_TYPE == "WEEKDAY" & period == "pm"]) / n_wd,
    am_out = sum(TOTAL_TAP_OUT_VOLUME[DAY_TYPE == "WEEKDAY" & period == "am"]) / n_wd,
    pm_out = sum(TOTAL_TAP_OUT_VOLUME[DAY_TYPE == "WEEKDAY" & period == "pm"]) / n_wd,
    .groups = "drop") |>
  mutate(daily_in  = (5 * wd_in + 2 * we_in) / 7,
         daily_out = (5 * wd_out + 2 * we_out) / 7)
audit$daily_total <- sum(rid$daily_in)
audit$wd_total    <- sum(rid$wd_in)
audit$we_total    <- sum(rid$we_in)

# Whether the months differ at station level: a stability check on the reference.
stab <- node |>
  left_join(select(lookup, raw, station_id), by = c("PT_CODE" = "raw")) |>
  left_join(calendar, by = c("YEAR_MONTH", "DAY_TYPE")) |>
  group_by(YEAR_MONTH, DAY_TYPE, station_id) |>
  summarise(v = sum(TOTAL_TAP_IN_VOLUME) / first(days), .groups = "drop") |>
  filter(DAY_TYPE == "WEEKDAY") |>
  select(-DAY_TYPE) |>
  pivot_wider(names_from = YEAR_MONTH, values_from = v)
audit$stab_aug_sep_cor <- cor(log(stab$`2026-08`), log(stab$`2026-09`))
audit$stab_aug_sep_mape <- median(abs(stab$`2026-09` / stab$`2026-08` - 1))

# ---- OD against the station file ----------------------------------------------
# The OD file records completed journeys, so it should fall a little short of tap-ins.
od_o <- od |> group_by(YEAR_MONTH, raw = ORIGIN_PT_CODE) |> summarise(od_o = sum(TOTAL_TRIPS), .groups = "drop")
od_d <- od |> group_by(YEAR_MONTH, raw = DESTINATION_PT_CODE) |> summarise(od_d = sum(TOTAL_TRIPS), .groups = "drop")
rec <- node |>
  group_by(YEAR_MONTH, raw = PT_CODE) |>
  summarise(tin = sum(TOTAL_TAP_IN_VOLUME), tout = sum(TOTAL_TAP_OUT_VOLUME), .groups = "drop") |>
  full_join(od_o, by = c("YEAR_MONTH", "raw")) |>
  full_join(od_d, by = c("YEAR_MONTH", "raw")) |>
  mutate(r_in = od_o / tin, r_out = od_d / tout)
audit$od_to_node <- sum(rec$od_o) / sum(rec$tin)
audit$od_ratio_range <- range(rec$r_in)
audit$od_ratio_low <- rec |> arrange(r_in) |> head(5) |> select(YEAR_MONTH, raw, tin, od_o, r_in)
audit$od_selfloop <- sum(od$TOTAL_TRIPS[od$ORIGIN_PT_CODE == od$DESTINATION_PT_CODE])

# ---- Station-to-station flows -------------------------------------------------
# Used by the destination-choice model. Average daily trips over the reference
# period, by the same day weighting as above.
od_st <- od |>
  filter(YEAR_MONTH %in% REF) |>
  left_join(select(lookup, raw, o = station_id), by = c("ORIGIN_PT_CODE" = "raw")) |>
  left_join(select(lookup, raw, d = station_id), by = c("DESTINATION_PT_CODE" = "raw")) |>
  mutate(period = per(TIME_PER_HOUR)) |>
  group_by(o, d) |>
  summarise(wd = sum(TOTAL_TRIPS[DAY_TYPE == "WEEKDAY"]) / n_wd,
            we = sum(TOTAL_TRIPS[DAY_TYPE != "WEEKDAY"]) / n_we,
            am = sum(TOTAL_TRIPS[DAY_TYPE == "WEEKDAY" & period == "am"]) / n_wd,
            .groups = "drop") |>
  mutate(daily = (5 * wd + 2 * we) / 7) |>
  filter(o != d)
audit$od_pairs <- nrow(od_st)
audit$od_pairs_possible <- audit$stations * (audit$stations - 1)

saveRDS(rid, file.path(DER, "ridership.rds"))
saveRDS(od_st, file.path(DER, "od_station.rds"))
saveRDS(audit, file.path(DER, "audit_ridership.rds"))

str(audit, max.level = 1, give.attr = FALSE)
print(audit$month_check); print(audit$multi_node); print(audit$od_ratio_low)
