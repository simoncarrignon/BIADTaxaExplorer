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

build_count_matrix <- function(subregions, count_column, use_logs = FALSE) {
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
    stop("Correspondence analysis needs at least two populated phase/area combinations and two populated taxon groups.")
  }

  if (use_logs) {
    counts <- log(counts + 1)
    notes <- c(notes, "A log(count + 1) transform was applied before correspondence analysis.")
  }

  list(counts = counts, notes = notes)
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

compute_analysis <- function(dataset, polygons, data_type, use_logs = FALSE) {
  normalized_polygons <- normalize_polygon_data(polygons, dataset)
  selected_polygon_count <- if (is.null(normalized_polygons)) 0 else nrow(normalized_polygons)
  count_column <- count_column_for_type(data_type)

  subregions <- prepare_analysis_data(dataset, normalized_polygons)
  count_bundle <- build_count_matrix(subregions, count_column, use_logs)
  used_polygon_count <- length(unique(subregions$new_area))
  notes <- count_bundle$notes

  if (selected_polygon_count > used_polygon_count) {
    notes <- c(
      notes,
      sprintf(
        "%s polygon(s) had no matching records for the current filters.",
        selected_polygon_count - used_polygon_count
      )
    )
  }

  ca_result <- FactoMineR::CA(count_bundle$counts, graph = FALSE)

  list(
    subregions = subregions,
    count_matrix = count_bundle$counts,
    ca = ca_result,
    combined_data = build_combined_data(count_bundle$counts, ca_result$row$coord),
    count_column = count_column,
    period_levels = resolve_period_levels(subregions$new_periods),
    summary = list(
      records = nrow(subregions),
      phase_count = nrow(count_bundle$counts),
      polygon_count = used_polygon_count,
      taxon_group_count = ncol(count_bundle$counts)
    ),
    notes = notes
  )
}
