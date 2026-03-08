normalize_polygon_data <- function(polygons, reference_data) {
  if (is.null(polygons)) {
    return(NULL)
  }

  if (!inherits(polygons, "sf")) {
    stop("Polygon data must be an sf object.")
  }

  if (nrow(polygons) == 0) {
    return(NULL)
  }

  geometry_types <- unique(as.character(sf::st_geometry_type(polygons, by_geometry = TRUE)))
  if (!all(grepl("POLYGON", geometry_types))) {
    stop("Polygon data must contain polygon geometries.")
  }

  normalized <- polygons
  reference_crs <- sf::st_crs(reference_data)
  polygon_crs <- sf::st_crs(normalized)
  same_crs <- isTRUE(all.equal(polygon_crs, reference_crs))

  if (!is.na(reference_crs) && is.na(polygon_crs)) {
    sf::st_crs(normalized) <- reference_crs
  } else if (!is.na(reference_crs) && !is.na(polygon_crs) && !same_crs) {
    normalized <- sf::st_transform(normalized, reference_crs)
  }

  geometry <- sf::st_geometry(normalized)
  attributes <- sf::st_drop_geometry(normalized)
  normalized <- sf::st_sf(attributes, geometry = geometry, crs = sf::st_crs(normalized))
  normalized$polygon_id <- seq_len(nrow(normalized))
  normalized
}

complete_polygon_attributes <- function(attribute_frame, all_names) {
  missing_names <- setdiff(all_names, names(attribute_frame))

  for (name in missing_names) {
    attribute_frame[[name]] <- NA
  }

  attribute_frame[, all_names, drop = FALSE]
}

bind_polygon_rows <- function(existing_polygons, new_polygons, reference_data) {
  existing_polygons <- normalize_polygon_data(existing_polygons, reference_data)
  new_polygons <- normalize_polygon_data(new_polygons, reference_data)

  if (is.null(existing_polygons)) {
    return(new_polygons)
  }
  if (is.null(new_polygons)) {
    return(existing_polygons)
  }

  existing_attributes <- sf::st_drop_geometry(existing_polygons)
  new_attributes <- sf::st_drop_geometry(new_polygons)
  all_attribute_names <- union(names(existing_attributes), names(new_attributes))

  existing_attributes <- complete_polygon_attributes(existing_attributes, all_attribute_names)
  new_attributes <- complete_polygon_attributes(new_attributes, all_attribute_names)

  combined_attributes <- rbind(existing_attributes, new_attributes)
  combined_geometry <- c(sf::st_geometry(existing_polygons), sf::st_geometry(new_polygons))
  combined_polygons <- sf::st_sf(
    combined_attributes,
    geometry = combined_geometry,
    crs = sf::st_crs(existing_polygons)
  )

  normalize_polygon_data(combined_polygons, reference_data)
}

append_drawn_polygon <- function(existing_polygons, feature, reference_data) {
  new_polygon <- geojsonsf::geojson_sf(jsonlite::toJSON(feature, auto_unbox = TRUE))
  bind_polygon_rows(existing_polygons, new_polygon, reference_data)
}

site_polygon_ids <- function(points, polygons, normalize = TRUE) {
  normalized_polygons <- if (normalize) {
    normalize_polygon_data(polygons, points)
  } else {
    polygons
  }

  if (is.null(normalized_polygons)) {
    return(rep(NA_integer_, nrow(points)))
  }

  first_intersection_id(sf::st_intersects(points, normalized_polygons))
}

site_colors <- function(points, polygons, palette = areapal, normalize = TRUE) {
  polygon_ids <- site_polygon_ids(points, polygons, normalize = normalize)
  color_index <- ifelse(is.na(polygon_ids), 1L, pmin(polygon_ids + 1L, length(palette)))
  palette[color_index]
}

site_popup_label <- function(points) {
  paste0(
    "<strong>",
    points$SiteName,
    "</strong><br/>",
    points$Country
  )
}

has_sf_features <- function(data) {
  inherits(data, "sf") && !is.null(data) && nrow(data) > 0
}

combined_map_geometry <- function(points = NULL, polygons = NULL) {
  geometries <- list()

  if (has_sf_features(points)) {
    geometries[[length(geometries) + 1]] <- sf::st_geometry(points)
  }

  if (has_sf_features(polygons)) {
    geometries[[length(geometries) + 1]] <- sf::st_geometry(polygons)
  }

  if (!length(geometries)) {
    return(NULL)
  }

  do.call(c, geometries)
}

fit_map_to_data <- function(map, points = NULL, polygons = NULL, pad_fraction = 0.05, min_span = 0.1) {
  geometry <- combined_map_geometry(points = points, polygons = polygons)
  if (is.null(geometry) || length(geometry) == 0) {
    return(map)
  }

  if (!is.na(sf::st_crs(geometry)) && !isTRUE(sf::st_is_longlat(geometry))) {
    geometry <- sf::st_transform(geometry, 4326)
  }

  bbox <- sf::st_bbox(geometry)
  x_span <- as.numeric(bbox[["xmax"]] - bbox[["xmin"]])
  y_span <- as.numeric(bbox[["ymax"]] - bbox[["ymin"]])
  x_pad <- max(x_span * pad_fraction, min_span)
  y_pad <- max(y_span * pad_fraction, min_span)

  leaflet::fitBounds(
    map,
    lng1 = as.numeric(bbox[["xmin"]] - x_pad),
    lat1 = as.numeric(bbox[["ymin"]] - y_pad),
    lng2 = as.numeric(bbox[["xmax"]] + x_pad),
    lat2 = as.numeric(bbox[["ymax"]] + y_pad)
  )
}
