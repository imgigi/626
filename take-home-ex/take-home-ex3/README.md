# Take-home Exercise 3 — Geospatial Analytics for Urban Sustainability

A probabilistic spatial interaction model of average daily tap-ins at Singapore's rail
stations, used to predict ridership at the Jurong Region Line and Cross Island Line
Phase 1 stations and to evaluate the impact on stations already open.

Part of the [ISSS626-GAA](../../index.qmd) coursework site. The rendered pages are
`take-home-ex3.html` (technical report) and `take-home-ex3-slides.html` (executive
summary); both carry an **EN / 中文** switch in the top-right corner.

## Layout

```
take-home-ex3.qmd            Technical report (pulls in _sections/)
take-home-ex3-slides.qmd     Executive summary, revealjs, 10 content slides
_sections/                   The report's fourteen sections
custom.scss, slides.scss     Report and slide styling
references.bib, apa.csl      Bibliography
R/utils.R                    Caching helper, palettes, plot themes, formatters
R/huff.R                     Huff model in matrix form, exact gradients, calibration
R/rail.R                     Rail graph, speeds, transfer penalty, travel times
scripts/                     The analysis pipeline (see below)
data/raw/                    Source data (large or keyed files are not committed)
data/derived/                Aggregated objects the report reads
data/cache/                  Slow results (git-ignored; seeded)
```

## Pipeline

Run from **this directory**:

```bash
Rscript scripts/01_markets.R      # HDB crosswalk, residents, hexagon grid
Rscript scripts/02_stations.R     # physical stations, attraction variables, activity markets
Rscript scripts/03_ridership.R    # calendar, reference period, OD against station volumes
Rscript scripts/04_frames.R       # the spatial interaction data frames
Rscript scripts/05_calibrate.R    # model ladder, hold-out, profiles, null, bootstrap, MCI check
Rscript scripts/06_scenario.R     # proposed network, impact, market areas, sensitivity
Rscript scripts/07_mci.R          # value-added: MCI on OD shares, transfers, link loads
Rscript scripts/08_report_data.R  # compact file the report reads
```

`05_calibrate.R` is the slow step (more than ten minutes from a cold cache). Then, from
the **site root**:

```bash
quarto render take-home-ex/take-home-ex3/take-home-ex3.qmd
quarto render take-home-ex/take-home-ex3/take-home-ex3-slides.qmd
```

## Data

| Data | Source | In the repository |
|---|---|---|
| Passenger Volume by Train Stations and by Origin Destination, Jul to Sep 2026 | LTA DataMall API | **No.** Needs your own AccountKey, and DataMall keeps only the latest three months. Run `scripts/00_fetch_lta.ps1` (hidden prompt) or `scripts/00_fetch_lta.sh`. |
| Train Station, Sep 2026; Bus Stop Location, Mar 2026 | LTA DataMall | **No.** Download from a browser into `data/raw/lta/`. |
| HDB Existing Building, Jun 2026 | data.gov.sg | **No** (57 MB). |
| OpenStreetMap, 31 Oct 2026 | Geofabrik | **No** (505 MB). |
| HDB Property Information, Oct 2026; residents by subzone, Jun 2024; Master Plan 2019 subzones; station coordinates | data.gov.sg, SingStat, URA, GitHub | Yes, in `data/raw/`. |
| Aggregated markets, frames and results | Built by the scripts | Yes, in `data/derived/`, except `model.rds` (53 MB of fitted matrices). The report renders from these alone. |

No account key is stored anywhere in this folder.

## Three things worth knowing before re-using the Huff code

1. The model needs a choice set. Without a maximum access distance, stations 10 km from
   anything new lose about 10% of their ridership, even though that specification has
   the best hold-out score.
2. Exact calibration, as `huff.fit()` in **MCI** does it, reads unexplained ridership as
   attractiveness. Kranji's mean catchment triples. Existing stations are therefore
   pivoted on observed ridership rather than calibrated.
3. `mci.fit()` and `huff.shares()` need complete rectangular tables. Station pairs with
   no trips and origins with no station in range stop them, so the log-centring is done
   by hand and checked against the package on a complete block.
