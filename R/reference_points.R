# Reference points from a user-supplied CSV -----------------------------------

lon_column_pattern <- paste(
  "(^|[^a-z])(lon(g(itude)?)?|lng|x|easting)([^a-z]|$)",
  "longitude",
  sep = "|"
)

lat_column_pattern <- paste(
  "(^|[^a-z])(lat(itude)?|y|northing)([^a-z]|$)",
  "latitude",
  sep = "|"
)

reference_label_pattern <- paste(
  "(^|[^a-z])(name|site|site_?id|label|id)([^a-z]|$)",
  "name",
  "site",
  sep = "|"
)

detect_reference_columns <- function(column_names) {
  column_names <- as.character(column_names)
  if (!length(column_names)) {
    return(list(lon = NULL, lat = NULL))
  }

  lon_hits <- grep(lon_column_pattern, column_names, ignore.case = TRUE, perl = TRUE)
  lat_hits <- grep(lat_column_pattern, column_names, ignore.case = TRUE, perl = TRUE)

  lon <- if (length(lon_hits)) column_names[lon_hits[[1]]] else NULL
  lat <- if (length(lat_hits)) column_names[lat_hits[[1]]] else NULL

  # A single column should never be used for both coordinates.
  if (!is.null(lon) && !is.null(lat) && identical(lon, lat)) {
    lat_hits <- setdiff(lat_hits, lon_hits[[1]])
    lat <- if (length(lat_hits)) column_names[lat_hits[[1]]] else NULL
  }

  list(lon = lon, lat = lat)
}

detect_reference_label_column <- function(column_names) {
  column_names <- as.character(column_names)
  hits <- grep(reference_label_pattern, column_names, ignore.case = TRUE, perl = TRUE)
  if (length(hits)) column_names[hits[[1]]] else NULL
}

build_reference_points <- function(data, lon_column, lat_column) {
  total <- nrow(data)

  if (is.null(lon_column) || !nzchar(lon_column) ||
      is.null(lat_column) || !nzchar(lat_column)) {
    return(list(
      points = NULL,
      total = total,
      skipped = 0L,
      message = "Choose a longitude and a latitude column."
    ))
  }

  if (identical(lon_column, lat_column)) {
    return(list(
      points = NULL,
      total = total,
      skipped = total,
      message = "Longitude and latitude must be different columns."
    ))
  }

  if (!lon_column %in% names(data) || !lat_column %in% names(data)) {
    return(list(
      points = NULL,
      total = total,
      skipped = total,
      message = "The selected columns are not present in this file."
    ))
  }

  lon <- suppressWarnings(as.numeric(as.character(data[[lon_column]])))
  lat <- suppressWarnings(as.numeric(as.character(data[[lat_column]])))
  valid <- is.finite(lon) & is.finite(lat) &
    lon >= -180 & lon <= 180 &
    lat >= -90 & lat <= 90
  skipped <- sum(!valid)

  if (!any(valid)) {
    return(list(
      points = NULL,
      total = total,
      skipped = skipped,
      message = "No rows contain valid longitude/latitude pairs."
    ))
  }

  filtered <- data[valid, , drop = FALSE]
  filtered[[lon_column]] <- lon[valid]
  filtered[[lat_column]] <- lat[valid]

  points <- sf::st_as_sf(
    filtered,
    coords = c(lon_column, lat_column),
    crs = 4326,
    remove = FALSE
  )

  label_column <- detect_reference_label_column(names(data))
  if (!is.null(label_column)) {
    labels <- as.character(filtered[[label_column]])
  } else {
    labels <- rep(NA_character_, nrow(filtered))
  }

  fallback <- paste("Point", which(valid))
  missing_labels <- is.na(labels) | !nzchar(labels)
  labels[missing_labels] <- fallback[missing_labels]
  points$.reference_label <- labels

  list(
    points = points,
    total = total,
    skipped = skipped,
    message = NULL
  )
}

reference_cross_icon <- function() {
  leaflet::makeIcon(
    iconUrl = "reference-cross.svg",
    iconWidth = 16,
    iconHeight = 16,
    iconAnchorX = 8,
    iconAnchorY = 8,
    popupAnchorX = 0,
    popupAnchorY = -8
  )
}

reference_popup_label <- function(points) {
  coordinates <- sf::st_coordinates(points)
  labels <- if (".reference_label" %in% names(points)) {
    as.character(points$.reference_label)
  } else {
    paste("Point", seq_len(nrow(points)))
  }

  missing_labels <- is.na(labels) | !nzchar(labels)
  labels[missing_labels] <- paste("Point", seq_along(labels))[missing_labels]

  paste0(
    "<strong>",
    htmltools::htmlEscape(labels),
    "</strong><br/>",
    sprintf("%.5f, %.5f", coordinates[, 1], coordinates[, 2])
  )
}
