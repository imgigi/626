# Take-home Exercise 2 — Geospatial Analytics for Peace

Local measures of spatial association and emerging hot spot analysis of armed
conflict in **Myanmar**, by township, January 2021 to September 2025.

Part of the [ISSS626-GAA](../../index.qmd) coursework site. The rendered pages
are `take-home-ex2.html` (technical report) and `take-home-ex2-slides.html`
(executive summary); both carry an **EN / 中文** switch in the top-right corner.

## Layout

```
take-home-ex2.qmd            Technical report (pulls in _sections/)
take-home-ex2-slides.qmd     Executive summary, revealjs, 10 content slides
_sections/                   The report's twelve sections
custom.scss, slides.scss     Report and slide styling
references.bib, apa.csl      Bibliography
R/utils.R                    Caching helper, palettes, plot themes, formatters
R/ehsa.R                     Space-time Gi*, period-conditional null, EHSA classes
scripts/                     The analysis pipeline (see below)
data/raw/                    MIMU v9.4 boundary zip; the ACLED extract (not committed)
data/derived/                Township-level objects the report reads
data/cache/                  Permutation results (git-ignored; seeded)
```

## Pipeline

Run from **this directory**:

```bash
Rscript scripts/01_boundaries.R   # MIMU townships, CRS choice, islands, neighbours
Rscript scripts/02_acled.R        # ACLED audit, spatial join, screening, panels
Rscript scripts/03_lmsa.R         # global and local Moran's I, Gi*, sensitivity
Rscript scripts/04_ehsa.R         # space-time cube, EHSA, sfdep check, null calibration
Rscript scripts/05_markov.R       # spatial Markov chain, FE logit, onsets, back-test
```

`04_ehsa.R` is the slow step (about 20 minutes from a cold cache). Then, from
the **site root**:

```bash
quarto render take-home-ex/take-home-ex2/take-home-ex2.qmd
quarto render take-home-ex/take-home-ex2/take-home-ex2-slides.qmd
```

## Data

| Data | Source | In the repository |
|---|---|---|
| ACLED events, Myanmar, Jan 2021 – Sep 2025 | ACLED, distributed through the course | **No.** ACLED's terms of use do not allow the records to be republished. Place the CSV anywhere under `data/raw/` to re-run the scripts. |
| Township boundaries v9.4 | MIMU GeoNode (WFS shapefile export) | Yes, as the downloaded zip, because MIMU replaces its releases. |
| Township-level aggregates and results | Built by the scripts | Yes, in `data/derived/`; the report renders from these alone. |

## Three things worth knowing before re-using the EHSA code

`emerging_hotspot_analysis()` in **sfdep 0.2.5** was checked against an
independent implementation (`R/ehsa.R`), and three defects were found:

1. It labels its rows in geometry order but computes them in alphabetical
   order of location, so results land on the wrong township unless the
   geometry is sorted by the location column.
2. Its permutation null draws neighbours from every period of the cube, so a
   national trend turns into local significance.
3. Its "intensifying hotspot" rule counts significant periods of either sign,
   so a persistent cold spot with a rising trend is labelled a hot spot.

The report (§7.3) documents each one on these data. `R/ehsa.R` keeps sfdep's
statistic and class rules and changes only the null and the one rule.
