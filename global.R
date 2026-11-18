# BIAD Taxa Explorer v2 -------------------------------------------------------
# Shared setup for the Shiny application.

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(leaflet)
  library(leaflet.extras)
  library(sf)
  library(geojsonsf)
  library(DT)
  library(ggplot2)
})

options(
  shiny.maxRequestSize = 64 * 1024^2,
  shiny.sanitize.errors = FALSE
)

# BIADconnect lives next to this project in the COREX workspace. Prefer the
# development copy so local changes are picked up, and fall back to the
# installed package if the sibling directory is not available.
if (dir.exists("../BIADconnect")) {
  if (requireNamespace("devtools", quietly = TRUE)) {
    suppressMessages(devtools::load_all("../BIADconnect/", quiet = TRUE))
  } else if (requireNamespace("BIADconnect", quietly = TRUE)) {
    library(BIADconnect)
  }
} else if (requireNamespace("BIADconnect", quietly = TRUE)) {
  library(BIADconnect)
}

# Constants -------------------------------------------------------------------
PHASE_SEPARATOR <- "||"
DEFAULT_POLYGON_PATH <- "groupings/spatial/polygon_18112026_SC_v2.gpkg"
DEFAULT_POLYGON_LAYER <- "polygon_18112026_SC_v2"

app_files <- c(
  "R/cache.R",
  "R/tools.R",
  "R/data_store.R",
  "R/reference_points.R",
  "R/map_utils.R",
  "R/analysis.R",
  "R/plots.R",
  "R/ui_components.R",
  "R/get.faunal.taxa.R",
  "R/get.botanical.taxa.R"
)

for (app_file in app_files) {
  if (!file.exists(app_file)) {
    stop("Missing required application file: ", app_file)
  }
  source(app_file, local = FALSE)
}

# Shared area palette. The first colour is reserved for sites outside any
# polygon so they always read as "unassigned" rather than a random area.
areapal <- c("#0f172a", grDevices::palette.colors(25, "Pastel 1", recycle = TRUE))
