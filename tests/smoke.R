# Lightweight smoke tests for the pure application helpers ---------------------
# Run with: Rscript tests/smoke.R

source("global.R")

fail <- function(...) stop(..., call. = FALSE)

# Grouping parsing -------------------------------------------------------------
sample_file <- "groupings/faunal_taxa/IV Domestic taxa & combined wild mammals.csv"
definition <- parse_grouping_file(sample_file)
if (!length(definition)) fail("Grouping parser returned no groups.")
if (!all(lengths(definition) > 0)) fail("Grouping parser produced empty groups.")

site_like <- tempfile(fileext = ".csv")
writeLines(c(
  "SiteID,SiteName,Region",
  "S1,Alpha,North",
  "S2,Beta,South"
), site_like)
if (is_valid_grouping_file(site_like)) fail("A site list was accepted as a grouping file.")

# Round-trip the editor format -------------------------------------------------
round_trip <- tempfile(fileext = ".csv")
writeLines(c(
  ",Example title,",
  ",,",
  "1,Group A,TAX1+TAX2",
  "2,Group B,TAX3"
), round_trip)
parsed_rows <- read_grouping_rows(round_trip)
groups <- parsed_rows$raw[parsed_rows$is_group, , drop = FALSE]
preamble <- parsed_rows$raw[!parsed_rows$is_group, , drop = FALSE]
groups[1, 2] <- "Renamed A"
write_grouping_rows(merge_grouping_rows(preamble, groups), round_trip)
reparsed <- parse_grouping_file(round_trip)
if (!"Renamed A" %in% names(reparsed)) fail("Grouping editor round-trip lost the edited label.")
if (length(reparsed) != 2) fail("Grouping editor round-trip changed the group count.")

# Phase labels -----------------------------------------------------------------
labels <- parse_phase_labels(c("5000||1", "Bronze Age||2", "Early Iron Age-3"))
if (!identical(labels$period, c("5000", "Bronze Age", "Early Iron Age"))) {
  fail("Phase label parsing returned the wrong periods.")
}
if (!identical(labels$area, c(1L, 2L, 3L))) fail("Phase label parsing returned the wrong areas.")

# PCA variance ----------------------------------------------------------------
set.seed(42)
counts <- matrix(
  sample(0:40, 48, replace = TRUE),
  nrow = 6,
  dimnames = list(paste0("p", 1:6), paste0("t", 1:8))
)
ordination <- run_ordination(counts, "pca", pca_transform = "raw", pca_scaling = "scale")
percentages <- ordination$eigenvalues[, "percentage of variance"]
if (any(!is.finite(percentages))) fail("PCA percentages are not finite.")
if (any(percentages < 0)) fail("PCA percentages are negative.")
if (sum(percentages) > 100 + 1e-8) fail("PCA percentages exceed 100%.")

# Small, fully explained case -------------------------------------------------
tiny <- matrix(
  c(1, 2, 3, 4, 5, 6, 2, 5, 1, 7, 4, 3, 9, 2, 6, 1, 8, 5),
  nrow = 6,
  byrow = TRUE
)
tiny_ordination <- run_ordination(tiny, "pca", pca_transform = "raw", pca_scaling = "none")
if (abs(sum(tiny_ordination$eigenvalues[, "percentage of variance"]) - 100) > 1e-6) {
  fail("PCA percentages should sum to 100 when all components are retained.")
}

# Reference points -------------------------------------------------------------
detected <- detect_reference_columns(c("id", "SiteName", "Longitude", "Latitude"))
if (!identical(detected$lon, "Longitude") || !identical(detected$lat, "Latitude")) {
  fail("Reference coordinate detection failed on Longitude/Latitude.")
}

detected <- detect_reference_columns(c("decimalLongitude", "decimalLatitude"))
if (!identical(detected$lon, "decimalLongitude") || !identical(detected$lat, "decimalLatitude")) {
  fail("Reference coordinate detection failed on decimalLongitude/decimalLatitude.")
}

reference <- build_reference_points(
  data.frame(id = 1:4, lon = c(10, NA, 200, 20), lat = c(50, 50, 50, 95)),
  "lon",
  "lat"
)
if (is.null(reference$points) || nrow(reference$points) != 1L) {
  fail("Reference point validation kept the wrong number of rows.")
}
if (reference$skipped != 3L) fail("Reference point skipped count is wrong.")
if (!identical(as.character(reference$points$.reference_label), "1")) {
  fail("Reference point label detection failed.")
}

same_column <- build_reference_points(data.frame(a = 1), "a", "a")
if (!is.null(same_column$points)) fail("Identical lon/lat columns were accepted.")

# UI construction --------------------------------------------------------------
ui <- app_ui()
if (!inherits(ui, "shiny.tag") && !inherits(ui, "shiny.tag.list")) {
  fail("app_ui() did not return a shiny tag.")
}

# Server / source parsing ------------------------------------------------------
parsed <- parse("server.R")
if (!length(parsed)) fail("server.R could not be parsed.")

cat("All smoke tests passed.\n")
