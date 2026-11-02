# =============================================================================
# R/rail.R -- the rail network as a graph
# Take-home Exercise 3 | ISSS626
#
# Lines are written out as ordered station codes. Station-to-station times are distance
# over a running speed that includes stops, and changing lines costs a fixed penalty.
# These are assumptions rather than measurements. Each station has separate entry and
# exit nodes, so that a path can only change lines by paying the transfer penalty, and
# never by passing through a station's gates.
# =============================================================================

SPEED <- c(MRT = 35, LRT = 20, JRL = 30)   # km/h, including dwell
TRANSFER_MIN <- 5

seqn <- function(p, v) paste0(p, v)
NETWORK <- list(
  NS = list(seqn("NS", c(1:5, 7:28)), "MRT"),
  EW = list(seqn("EW", 1:33), "MRT"),
  CG = list(c("EW4", "CG1", "CG2"), "MRT"),
  NE = list(seqn("NE", c(1, 3:18)), "MRT"),
  CC = list(c(seqn("CC", c(1:17, 19:32)), "CE2", "CE1", "CC4"), "MRT"),
  DT = list(seqn("DT", 1:35), "MRT"),
  TE = list(seqn("TE", c(1:9, 11:20, 22:29)), "MRT"),
  BP = list(c(seqn("BP", 1:13), "BP6"), "LRT"),
  SE = list(c("STC", seqn("SE", 1:5), "STC"), "LRT"),
  SW = list(c("STC", seqn("SW", 1:8), "STC"), "LRT"),
  PE = list(c("PTC", seqn("PE", 1:7), "PTC"), "LRT"),
  PW = list(c("PTC", seqn("PW", 1:7), "PTC"), "LRT")
)
# The Jurong Region Line's three arms meet at Tengah and Bahar Junction. Its medium
# capacity trains are given a lower speed than heavy rail.
PROPOSED <- list(
  JS = list(seqn("JS", 1:12), "JRL"),
  JW = list(c("JS7", seqn("JW", 1:5)), "JRL"),
  JE = list(c("JS3", seqn("JE", 1:7)), "JRL"),
  CR = list(seqn("CR", 2:13), "MRT")
)

# cx: one row per code with station_id, x, y
build_graph <- function(lines, cx, transfer = TRANSFER_MIN) {
  link <- purrr::map_dfr(names(lines), function(nm) {
    s <- lines[[nm]][[1]]
    tibble::tibble(from = s[-length(s)], to = s[-1], line = nm, mode = lines[[nm]][[2]])
  })
  link <- link |>
    dplyr::left_join(dplyr::select(cx, code, fx = x, fy = y), by = c("from" = "code")) |>
    dplyr::left_join(dplyr::select(cx, code, tx = x, ty = y), by = c("to" = "code")) |>
    dplyr::mutate(km = sqrt((fx - tx)^2 + (fy - ty)^2) / 1000,
                  minutes = km / SPEED[mode] * 60, type = "ride") |>
    dplyr::select(from, to, line, minutes, km, type)
  used_st <- dplyr::filter(cx, code %in% unique(c(link$from, link$to)))
  gates <- used_st |>
    dplyr::transmute(from = paste0("IN_", station_id), to = code, line = NA_character_,
                     minutes = 0, km = 0, type = "enter")
  exits <- used_st |>
    dplyr::transmute(from = code, to = paste0("OUT_", station_id), line = NA_character_,
                     minutes = 0, km = 0, type = "exit")
  xfer <- used_st |>
    dplyr::inner_join(used_st, by = "station_id", suffix = c("_a", "_b"), relationship = "many-to-many") |>
    dplyr::filter(code_a != code_b) |>
    dplyr::transmute(from = code_a, to = code_b, line = NA_character_, minutes = transfer,
                     km = 0, type = "transfer")
  rides <- dplyr::bind_rows(link, dplyr::rename(link, from = to, to = from))
  igraph::graph_from_data_frame(dplyr::bind_rows(rides, gates, exits, xfer), directed = TRUE)
}

rail_matrix <- function(g, ids) {
  m <- igraph::distances(g, v = paste0("IN_", ids), to = paste0("OUT_", ids), mode = "out",
                         weights = igraph::E(g)$minutes)
  dimnames(m) <- list(ids, ids)
  m
}
