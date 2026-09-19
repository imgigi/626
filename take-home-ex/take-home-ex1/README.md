# Take-home Exercise 1 — Fire on Drained Peat

Spatial and spatio-temporal point pattern analysis of satellite active-fire
detections in **Pulang Pisau Regency, Central Kalimantan**, 4–11 September 2026.

Part of the [ISSS626-GAA](../../index.qmd) coursework site. The rendered pages
are `take-home-ex1.html` (technical report) and `take-home-ex1-slides.html`
(executive summary); both carry an **EN / 中文** switch in the top-right corner.

## Layout

```
take-home-ex1.qmd            Technical report (pulls in _sections/)
take-home-ex1-slides.qmd     Executive summary, revealjs, 10 content slides
_sections/                   The report's eleven sections
_include/                    EN/中文 language switch (head + body injections)
custom.scss, slides.scss     Report and slide styling
references.bib, apa.csl      Bibliography
R/utils.R                    Caching helper, plot theme, formatters
scripts/                     The analysis pipeline (see below)
data/raw/                    Source data as downloaded
data/derived/                Prepared objects the report reads
data/cache/                  Monte-Carlo results (git-ignored; seeded)
```

## Pipeline

Run from **this directory** — every path in the scripts and in the report is
relative to it, matching the convention used by the hands-on exercises.

```bash
Rscript scripts/03_prepare.R        # QA, CRS, point-event definitions, ppp objects
Rscript scripts/04_covariates.R     # distance-to-canal/river/road/settlement
Rscript scripts/05_firstorder.R     # quadrat, KDE, bandwidth ladder, rhohat
Rscript scripts/05b_stkde.R         # daily and spatio-temporal kernels
Rscript scripts/06_secondorder.R    # G/F/J, L, Linhom ladder, pcf, Knox
Rscript scripts/07_valueadded.R     # ppm and kppm point-process models
Rscript scripts/08_slim_for_repo.R  # shrink saved results to what the report reads
```

`06` and `07` are the expensive steps (a few tens of minutes of Monte-Carlo
simulation). Results are cached in `data/cache/`; `set.seed(1234)` is set at the
top of every script, so a recomputation reproduces the same values.

Then, from the **site root**:

```bash
quarto render take-home-ex/take-home-ex1/take-home-ex1.qmd
quarto render take-home-ex/take-home-ex1/take-home-ex1-slides.qmd
```

### Two things worth knowing before editing

- **`freeze` is off for this document.** Quarto's freeze cache does not notice
  edits inside files pulled in with `{{< include >}}`, so leaving it on would
  render stale output after a change to `_sections/`. The cost is that the
  report re-executes (~8 minutes) on every render.
- **`08_slim_for_repo.R` must be re-run after `06`/`07`.** Those scripts save
  their full working objects — `va_results.rds` comes out at ~121 MB, above
  GitHub's 100 MB per-file limit. The slimming step reduces it to ~1 MB by
  dropping the simulated curves kept by `savefuns = TRUE` and storing the two
  scalars the report needed out of the fitted models. It is lossless for every
  number and figure in the report.

## Data

The FIRMS open archive is a **rolling 7-day window**, so the exact records
analysed here cannot be recovered by re-running the download. The analysed CSV
(`data/raw/firms_pulangpisau_2026-09-04_2026-09-11.csv`) and the OSM extracts
are therefore committed. The bulk national boundary datasets are not — they are
re-downloadable from stable URLs, and the four small objects distilled from them
are committed instead (see `scripts/02b_context_layers.R`).

`scripts/01_download_firms.R` uses the authenticated FIRMS Area API and accepts
an arbitrary date range; supplying a free `MAP_KEY` (in `.firms_key` or the
`FIRMS_MAP_KEY` environment variable) and rerunning it over a full season
reproduces every result below without changing any other file. That is the
single highest-value extension of this work, since the eight-day window is what
forces the report to disclaim all temporal-trend findings.

| Data | Source | Licence |
|---|---|---|
| VIIRS / MODIS active fire | NASA FIRMS / LANCE | Open (NASA data policy) |
| Indonesia admin boundaries (COD-AB) | Badan Pusat Statistik / UN OCHA via HDX | CC BY 3.0 IGO |
| geoBoundaries ADM2 | William & Mary geoLab | CC BY 4.0 |
| Canals, rivers, roads, places | OpenStreetMap contributors | ODbL 1.0 |
