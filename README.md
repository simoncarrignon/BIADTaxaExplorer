## BIAD Taxa Explorer v2

A forked, optimised version of the BIAD Taxa Explorer Shiny app. It keeps the
original analysis workflow but speeds up data access, fixes a few correctness
issues, and rebuilds the interface around a clearer four-step workflow.

### What changed

**Performance**

- BIAD records are cached to `cache/biad_data.rds` and reused for 12 hours by
  default (`BIAD_CACHE_HOURS` overrides this). Startup drops from a database
  round-trip to a local read on subsequent launches.
- Grouping files are parsed once and cached; the cache is invalidated
  automatically when a file is uploaded, edited, or changed on disk.
- Eliminated repeated `st_intersects()` calls: site colours are memoised on the
  current data and polygon set.
- Removed the unused citation columns from the faunal query.
- The map no longer re-fits its bounds on every repaint; use **Zoom to data**
  or change the polygons to re-frame.

**Correctness**

- PCA explained-variance percentages now divide by the total variance of the
  input matrix, not just the retained components.
- Grouping discovery ignores non-grouping CSVs (for example site lists).
- The group editor now preserves the title/preamble rows and writes the file
  without inserting a new header row.
- Uploads are reduced to their basename before writing.
- NULL geometries/site tables are guarded in the map and analysis paths.
- Phase identifiers use a `||` separator so period labels containing `-` parse
  correctly (older hyphenated ids still load).

**Interface**

- Four-step accordion workflow with a persistent run card.
- Live data status bar (source, cache age, record counts).
- Metric tiles for dataset/sites/records/polygons/grouping.
- Spinners on maps and plots, clear empty states and inline validation.
- Quick time-range presets, period validation, and helper text throughout.
- Result tabs now include an Overview, Diagnostics split out from Dimensions,
  and a searchable/scrollable data preview.
- The map can show only sites inside the study areas, with the metric tile
  reporting "shown / total" sites.
- **Reference points (CSV)**: upload a CSV of coordinates under Study area.
  Longitude/latitude columns are detected from their names, can be overridden
  with two selectors, and valid points are drawn as red crosses. Invalid or
  out-of-range rows are reported and skipped.
- Culture composition is limited to the top cultures by total counts by default
  (adjustable), so the chart stays readable on rich datasets.
- The default grouping is a richer scheme (VIII All Mammal taxa / I ABot all
  groups) so the ordination map has at least two dimensions on first run.
- Download buttons only appear when there is something to download.
- Keyboard focus styling, responsive layout, and reduced-motion support.

### Run locally

```bash
Rscript -e "shiny::shinyAppDir('.', options = list(port = 1112))"
```

The app expects `BIADconnect` at `../BIADconnect/` (the fork lives beside the
original project so this resolves unchanged).

Useful options:

```bash
# Refresh the cached data on next launch
BIAD_CACHE_HOURS=0 Rscript -e "shiny::shinyAppDir('.', options = list(port = 1112))"

# Or remove the cache explicitly
rm -f cache/biad_data.rds
```

### Project structure

- `global.R` — libraries, BIADconnect loading, constants.
- `ui.R` — calls `app_ui()`.
- `server.R` — state, data loading, map syncing, analysis, downloads, grouping
  management.
- `R/cache.R` — small in-memory cache helpers.
- `R/data_store.R` — cached database access and data summaries.
- `R/tools.R` — grouping parsing/caching, taxon mapping, dataset helpers.
- `R/map_utils.R` — polygon normalisation and map helpers.
- `R/analysis.R` — period slicing, aggregation, CA/PCA preparation.
- `R/plots.R` — base and ggplot2 plot helpers.
- `R/ui_components.R` — theme, reusable widgets, and the full UI.
- `www/app.css` — interface styling.
