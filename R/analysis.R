validate_auto_period_inputs <- function(start, end, len) {
  if (is.null(start) || is.null(end) || is.null(len)) {
    stop("Automatic time slicing needs a start, an end, and a duration.")
  }

  if (start <= end) {
    stop("For BP slices, the older bound must be greater than the younger bound.")
  }

  if (len <= 0) {
    stop("Time-slice duration must be greater than zero.")
  }
}

split_periods <- function(dataset, selection, start, end, len) {
  if (identical(selection, "auto")) {
    validate_auto_period_inputs(start, end, len)

    breaks <- seq(end, start, by = len)
    if (tail(breaks, 1) != start) {
      breaks <- c(breaks, start)
    }
    breaks <- sort(unique(breaks))

    if (length(breaks) < 2) {
      stop("Automatic time slicing needs at least two break points.")
    }

    labels <- as.character(breaks[-1])
    periods <- cut(
      dataset$GMM,
      breaks = breaks,
      labels = labels,
      include.lowest = TRUE,
      right = TRUE
    )

    return(factor(as.character(periods), levels = rev(labels), ordered = TRUE))
  }

  if (identical(selection, "groupings/periods/periods.csv")) {
    grouped_periods <- groupPeriod(dataset$Period, selection)
    return(
      factor(
        grouped_periods,
        levels = c("Neolithic", "Eneolithic", "Bronze Age", "Early Iron Age"),
        ordered = TRUE
      )
    )
  }

  stop("Unknown period grouping option selected.")
}

parse_phase_labels <- function(labels) {
  label_parts <- strsplit(labels, "-", fixed = TRUE)

  data.frame(
    period = vapply(
      label_parts,
      function(parts) {
        paste(parts[seq_len(length(parts) - 1L)], collapse = "-")
      },
      character(1)
    ),
    area = as.integer(vapply(label_parts, function(parts) parts[[length(parts)]], character(1))),
    stringsAsFactors = FALSE
  )
}

resolve_period_levels <- function(period_values) {
  if (is.factor(period_values)) {
    return(levels(droplevels(period_values)))
  }

  unique_periods <- unique(as.character(period_values))
  if (!length(unique_periods)) {
    return(character(0))
  }

  is_numeric_period <- all(grepl("^-?[0-9]+(?:\\.[0-9]+)?$", unique_periods))
  if (is_numeric_period) {
    return(unique_periods[order(as.numeric(unique_periods), decreasing = TRUE)])
  }

  unique_periods
}

phase_levels_from_data <- function(period_values, area_ids) {
  period_levels <- resolve_period_levels(period_values)
  if (!length(period_levels) || !length(area_ids)) {
    return(character(0))
  }

  as.vector(unlist(
    lapply(period_levels, function(period) paste(period, area_ids, sep = "-")),
    use.names = FALSE
  ))
}

prepare_analysis_data <- function(dataset, polygons) {
  normalized_polygons <- normalize_polygon_data(polygons, dataset)

  if (is.null(normalized_polygons)) {
    stop("Draw or upload at least one polygon before running the analysis.")
  }

  prepared <- dataset
  prepared$new_area <- first_intersection_id(sf::st_intersects(prepared, normalized_polygons))
  prepared <- prepared[
    !is.na(prepared$new_area) &
      !is.na(prepared$new_periods) &
      !is.na(prepared$new_txgroups),
    ,
    drop = FALSE
  ]

  if (nrow(prepared) == 0) {
    stop("No records fall inside the selected polygons for the current settings.")
  }

  area_ids <- sort(unique(prepared$new_area))
  phase_values <- paste(as.character(prepared$new_periods), prepared$new_area, sep = "-")
  phase_levels <- phase_levels_from_data(prepared$new_periods, area_ids)
  prepared$phase <- factor(phase_values, levels = phase_levels, ordered = TRUE)

  prepared
}

analysis_method_label <- function(analysis_method) {
  switch(
    analysis_method,
    "ca" = "Correspondence analysis",
    "pca" = "Principal component analysis",
    stop("Unknown analysis method selected.")
  )
}

ordination_diagnostic_label <- function(analysis_method) {
  switch(
    analysis_method,
    "ca" = "Percent inertia explained",
    "pca" = "Percent variance explained",
    stop("Unknown analysis method selected.")
  )
}

build_count_matrix <- function(subregions, count_column) {
  counts <- tapply(
    subregions[[count_column]],
    list(subregions$phase, subregions$new_txgroups),
    sum,
    na.rm = TRUE
  )
  counts <- as.matrix(counts)
  notes <- character(0)

  if (anyNA(counts)) {
    counts[is.na(counts)] <- 0
    notes <- c(notes, "Missing phase/taxon combinations were treated as zero.")
  }

  empty_rows <- rowSums(counts, na.rm = TRUE) == 0
  if (any(empty_rows)) {
    counts <- counts[!empty_rows, , drop = FALSE]
    notes <- c(notes, sprintf("%s empty phase rows were removed.", sum(empty_rows)))
  }

  empty_cols <- colSums(counts, na.rm = TRUE) == 0
  if (any(empty_cols)) {
    counts <- counts[, !empty_cols, drop = FALSE]
    notes <- c(notes, sprintf("%s empty taxon groups were removed.", sum(empty_cols)))
  }

  if (nrow(counts) < 2 || ncol(counts) < 2) {
    stop("Ordination needs at least two populated phase/area combinations and two populated taxon groups.")
  }

  list(counts = counts, notes = notes)
}

run_ordination <- function(count_matrix, analysis_method, pca_transform = "log1p", pca_scaling = "scale") {
  notes <- character(0)
  ordination_matrix <- count_matrix

  if (identical(analysis_method, "ca")) {
    model <- FactoMineR::CA(ordination_matrix, graph = FALSE)
    return(list(
      method = analysis_method,
      method_label = analysis_method_label(analysis_method),
      diagnostic_label = ordination_diagnostic_label(analysis_method),
      matrix = ordination_matrix,
      row_coordinates = as.matrix(model$row$coord),
      column_coordinates = as.matrix(model$col$coord),
      eigenvalues = model$eig,
      notes = notes
    ))
  }

  if (!identical(analysis_method, "pca")) {
    stop("Unknown analysis method selected.")
  }

  if (identical(pca_transform, "log1p")) {
    ordination_matrix <- log(ordination_matrix + 1)
    notes <- c(notes, "PCA used log(count + 1) transformed counts.")
  } else if (identical(pca_transform, "raw")) {
    notes <- c(notes, "PCA used raw counts.")
  } else {
    stop("Unknown PCA input transform selected.")
  }

  column_variance <- apply(ordination_matrix, 2, stats::var, na.rm = TRUE)
  zero_variance_columns <- !is.finite(column_variance) | column_variance == 0
  if (any(zero_variance_columns)) {
    ordination_matrix <- ordination_matrix[, !zero_variance_columns, drop = FALSE]
    notes <- c(notes, sprintf("%s taxon groups with zero variance were removed before PCA.", sum(zero_variance_columns)))
  }

  if (nrow(ordination_matrix) < 2 || ncol(ordination_matrix) < 2) {
    stop("PCA needs at least two populated phase/area combinations and two taxon groups with non-zero variance.")
  }

  if (identical(pca_scaling, "scale")) {
    center_columns <- TRUE
    scale_columns <- TRUE
    notes <- c(notes, "PCA was run on centered and scaled taxon-group columns.")
  } else if (identical(pca_scaling, "center")) {
    center_columns <- TRUE
    scale_columns <- FALSE
    notes <- c(notes, "PCA was run on centered taxon-group columns without scaling.")
  } else if (identical(pca_scaling, "none")) {
    center_columns <- FALSE
    scale_columns <- FALSE
    notes <- c(notes, "PCA was run without centering or scaling taxon-group columns.")
  } else {
    stop("Unknown PCA scaling option selected.")
  }

  pca_input <- as.matrix(ordination_matrix)
  model <- stats::prcomp(
    pca_input,
    center = center_columns,
    scale. = scale_columns,
    rank. = min(5, nrow(ordination_matrix) - 1L, ncol(ordination_matrix))
  )
  retained_dimensions <- seq_len(ncol(model$rotation))
  retained_sdev <- model$sdev[retained_dimensions]
  eigenvalues <- retained_sdev^2
  percentages <- eigenvalues / sum(eigenvalues) * 100
  cumulative <- cumsum(percentages)
  eig <- cbind(
    eigenvalue = eigenvalues,
    "percentage of variance" = percentages,
    "cumulative percentage of variance" = cumulative
  )
  rownames(eig) <- paste("comp", seq_along(eigenvalues))
  column_coordinates <- sweep(model$rotation, 2, retained_sdev, "*")

  if (identical(pca_scaling, "scale")) {
    notes <- c(notes, "Scaled PCA gives all taxon-group columns unit variance before decomposition.")
  } else if (identical(pca_scaling, "center")) {
    notes <- c(notes, "Centered-only PCA removes column means but keeps original variance differences.")
  } else if (identical(pca_scaling, "none")) {
    notes <- c(notes, "Uncentered PCA keeps both column means and variance differences.")
  }

  list(
    method = analysis_method,
    method_label = analysis_method_label(analysis_method),
    diagnostic_label = ordination_diagnostic_label(analysis_method),
    matrix = ordination_matrix,
    row_coordinates = as.matrix(model$x),
    column_coordinates = as.matrix(column_coordinates),
    eigenvalues = eig,
    notes = notes
  )
}

build_combined_data <- function(count_matrix, row_coordinates) {
  coordinates <- as.matrix(row_coordinates)
  if (is.null(colnames(coordinates))) {
    colnames(coordinates) <- paste("Dim", seq_len(ncol(coordinates)))
  }

  phase_info <- parse_phase_labels(rownames(count_matrix))

  base_frame <- data.frame(
    `Phase start` = phase_info$period,
    `Polygon ID` = phase_info$area,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  cbind(
    base_frame,
    as.data.frame(count_matrix, check.names = FALSE),
    as.data.frame(coordinates, check.names = FALSE)
  )
}

phase_assignment_key <- function(dataset) {
  key_columns <- intersect(
    c("SiteId", "PhaseId", "new_area", "new_periods", "phase"),
    names(dataset)
  )

  if (!length(key_columns)) {
    stop("Phase assignment export needs site/phase identifiers.")
  }

  do.call(
    paste,
    c(
      lapply(dataset[key_columns], function(column_values) {
        ifelse(is.na(column_values), "", as.character(column_values))
      }),
      sep = "||"
    )
  )
}

build_phase_assignment_taxon_counts <- function(subregions, count_column, phase_ids) {
  phase_frame <- sf::st_drop_geometry(subregions)
  taxon_columns <- as.character(phase_frame$TaxonCode)
  if ("Dataset" %in% names(phase_frame)) {
    taxon_columns <- paste(as.character(phase_frame$Dataset), taxon_columns, sep = ": ")
  }

  source_record_ids <- if ("SourceRecordId" %in% names(phase_frame)) {
    as.character(phase_frame$SourceRecordId)
  } else if ("FaunalSpeciesID" %in% names(phase_frame)) {
    paste("Faunal", as.character(phase_frame$FaunalSpeciesID), sep = ":")
  } else if ("SampleID" %in% names(phase_frame)) {
    paste("Botanical", as.character(phase_frame$SampleID), sep = ":")
  } else {
    paste("row", seq_len(nrow(phase_frame)), sep = ":")
  }

  count_frame <- unique(
    data.frame(
      phase_id = phase_assignment_key(phase_frame),
      taxon_column = taxon_columns,
      source_record_id = source_record_ids,
      count_value = phase_frame[[count_column]],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  )

  if (!nrow(count_frame)) {
    return(data.frame(row.names = seq_along(phase_ids)))
  }

  aggregated_counts <- stats::aggregate(
    count_value ~ phase_id + taxon_column,
    data = count_frame,
    FUN = sum
  )

  taxon_levels <- sort(unique(aggregated_counts$taxon_column))
  taxon_matrix <- matrix(
    NA_real_,
    nrow = length(phase_ids),
    ncol = length(taxon_levels),
    dimnames = list(NULL, taxon_levels)
  )

  taxon_matrix[
    cbind(
      match(aggregated_counts$phase_id, phase_ids),
      match(aggregated_counts$taxon_column, taxon_levels)
    )
  ] <- aggregated_counts$count_value

  as.data.frame(taxon_matrix, check.names = FALSE)
}

build_phase_assignment_data <- function(subregions, count_column) {
  coordinates <- sf::st_coordinates(sf::st_transform(subregions, 4326))
  phase_frame <- cbind(
    sf::st_drop_geometry(subregions),
    Longitude = coordinates[, "X"],
    Latitude = coordinates[, "Y"]
  )
  culture_column <- intersect(c("Culture", "Culture1"), names(phase_frame))
  selected_columns <- intersect(
    c(
      "SiteId",
      "SiteName",
      "Country",
      "Longitude",
      "Latitude",
      "PhaseId",
      "Period",
      "GMM",
      "GMS",
      "new_area",
      "new_periods",
      "phase"
    ),
    names(phase_frame)
  )
  selected_columns <- c(selected_columns, culture_column[1])
  selected_columns <- unique(selected_columns)

  phase_frame$phase_id <- phase_assignment_key(phase_frame)
  phase_frame <- unique(phase_frame[, c("phase_id", selected_columns), drop = FALSE])
  column_labels <- c(
    SiteId = "Site ID",
    SiteName = "Site name",
    Country = "Country",
    Longitude = "Longitude",
    Latitude = "Latitude",
    PhaseId = "Phase ID",
    Period = "Original period",
    Culture = "Culture",
    Culture1 = "Culture",
    GMM = "Raw time (GMM)",
    GMS = "Raw time uncertainty (GMS)",
    new_area = "Polygon ID",
    new_periods = "Time bin",
    phase = "Grouped phase"
  )
  metadata_columns <- selected_columns
  names(phase_frame)[match(metadata_columns, names(phase_frame))] <- unname(column_labels[metadata_columns])

  if ("Time bin" %in% names(phase_frame)) {
    phase_frame[["Time bin"]] <- as.character(phase_frame[["Time bin"]])
  }
  if ("Grouped phase" %in% names(phase_frame)) {
    phase_frame[["Grouped phase"]] <- as.character(phase_frame[["Grouped phase"]])
  }

  ordering_columns <- intersect(
    c("Polygon ID", "Raw time (GMM)", "Site ID", "Phase ID"),
    names(phase_frame)
  )

  if (length(ordering_columns)) {
    ordering_frame <- phase_frame[, ordering_columns, drop = FALSE]
    if ("Raw time (GMM)" %in% names(ordering_frame) && is.numeric(ordering_frame[["Raw time (GMM)"]])) {
      ordering_frame[["Raw time (GMM)"]] <- -ordering_frame[["Raw time (GMM)"]]
    }
    phase_frame <- phase_frame[do.call(order, ordering_frame), , drop = FALSE]
  }

  taxon_counts <- build_phase_assignment_taxon_counts(
    subregions = subregions,
    count_column = count_column,
    phase_ids = phase_frame$phase_id
  )
  phase_frame <- cbind(
    phase_frame[, setdiff(names(phase_frame), "phase_id"), drop = FALSE],
    taxon_counts
  )

  rownames(phase_frame) <- NULL
  phase_frame
}

compute_analysis <- function(dataset, polygons, data_type, analysis_method = "ca", pca_transform = "log1p", pca_scaling = "scale") {
  normalized_polygons <- normalize_polygon_data(polygons, dataset)
  selected_polygon_count <- if (is.null(normalized_polygons)) 0 else nrow(normalized_polygons)
  count_column <- count_column_for_type(data_type)
  dataset_notes <- attr(dataset, "analysis_notes", exact = TRUE)
  if (is.null(dataset_notes)) {
    dataset_notes <- character(0)
  }

  subregions <- prepare_analysis_data(dataset, normalized_polygons)
  count_bundle <- build_count_matrix(subregions, count_column)
  ordination_bundle <- run_ordination(
    count_matrix = count_bundle$counts,
    analysis_method = analysis_method,
    pca_transform = pca_transform,
    pca_scaling = pca_scaling
  )
  used_polygon_count <- length(unique(subregions$new_area))
  notes <- c(count_bundle$notes, dataset_notes, ordination_bundle$notes)

  if (identical(data_type, "Combined")) {
    notes <- c(notes, "Combined mode merges faunal NISP and botanical counts into the same matrix.")
  }

  if (selected_polygon_count > used_polygon_count) {
    notes <- c(
      notes,
      sprintf(
        "%s polygon(s) had no matching records for the current filters.",
        selected_polygon_count - used_polygon_count
      )
    )
  }

  list(
    subregions = subregions,
    count_matrix = ordination_bundle$matrix,
    ordination = ordination_bundle,
    phase_assignments = build_phase_assignment_data(subregions, count_column),
    combined_data = build_combined_data(ordination_bundle$matrix, ordination_bundle$row_coordinates),
    count_column = count_column,
    period_levels = resolve_period_levels(subregions$new_periods),
    summary = list(
      records = nrow(subregions),
      phase_count = nrow(ordination_bundle$matrix),
      polygon_count = used_polygon_count,
      taxon_group_count = ncol(ordination_bundle$matrix),
      analysis_method = ordination_bundle$method_label
    ),
    notes = notes
  )
}
