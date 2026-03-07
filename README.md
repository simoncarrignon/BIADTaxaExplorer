## BIAD Taxa Explorer

Shiny app for exploring faunal and botanical taxa inside custom map polygons, grouping them through predefined taxonomic schemes, and comparing spatial-temporal structure with correspondence analysis.

### Updated UI goals

- Guided workflow in the sidebar: dataset, study areas, time controls, then analysis/export.
- Map-first layout with live site recoloring when polygons change.
- Result tabs that only refresh after `Run analysis`, so the UI stays responsive.
- Clearer helper functions in `R/` so data prep, map logic, analysis, and UI building are separated.

### Project structure

- `ui.R` wires the app layout through `app_ui()`.
- `server.R` coordinates reactive state, map syncing, analysis runs, and downloads.
- `R/tools.R` contains grouping and small data utilities.
- `R/map_utils.R` contains polygon normalization and point-color logic.
- `R/analysis.R` contains period splitting, aggregation, and correspondence-analysis prep.
- `R/plots.R` contains plot helpers.
- `R/ui_components.R` contains the reusable UI cards and status widgets.
- `global.R` loads BIAD data through the sibling `../BIADconnect/` package.

### Run locally

```bash
Rscript -e "shiny::shinyAppDir('.', options = list(port = 1111))"
```

This app expects the `BIADconnect` package to be available at `../BIADconnect/`.
