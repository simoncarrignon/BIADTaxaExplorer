empty_plot_message <- function(message) {
  plot.new()
  text(0.5, 0.5, message, cex = 1.05)
}

taxon_group_palette <- function(group_names) {
  group_names <- sort(unique(as.character(group_names)))
  colors <- palette.colors(max(length(group_names), 1), "Pastel 2", recycle = TRUE)
  stats::setNames(colors[seq_along(group_names)], group_names)
}

legend_text_size <- function(labels, base = 1.02, minimum = 0.82) {
  entry_count <- length(labels)
  longest_label <- max(nchar(as.character(labels)), 0)
  size <- base

  if (entry_count > 8) {
    size <- size - 0.04 * (entry_count - 8)
  }
  if (longest_label > 18) {
    size <- size - 0.01 * (longest_label - 18)
  }

  max(minimum, min(base, size))
}

legend_column_count <- function(labels) {
  entry_count <- length(labels)
  if (entry_count > 18) {
    return(3)
  }
  if (entry_count > 9) {
    return(2)
  }
  1
}

separate_legend_column_count <- function(labels) {
  entry_count <- length(labels)
  if (entry_count > 15) {
    return(4)
  }
  if (entry_count > 7) {
    return(3)
  }
  if (entry_count > 3) {
    return(2)
  }
  1
}

clear_plot_legend <- function(position = "topright", labels, ..., cex = NULL, ncol = NULL, inset = 0) {
  labels <- as.character(labels)
  if (is.null(cex)) {
    cex <- legend_text_size(labels)
  }
  if (is.null(ncol)) {
    ncol <- legend_column_count(labels)
  }

  legend(
    position,
    legend = labels,
    ...,
    bg = "#ffffff",
    box.col = "#cbd5e1",
    text.col = "#111827",
    border = "#cbd5e1",
    bty = "o",
    cex = cex,
    ncol = ncol,
    x.intersp = 0.8,
    y.intersp = 1.05,
    inset = inset
  )
}

plotTaxonGroupLegend <- function(dataset) {
  group_names <- sort(unique(as.character(stats::na.omit(dataset$new_txgroups))))

  if (!length(group_names)) {
    empty_plot_message("No taxon groups are available for the current selection.")
    return(invisible(NULL))
  }

  group_colors <- taxon_group_palette(group_names)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  par(mar = c(0.5, 0.5, 1.2, 0.5), xpd = NA)
  plot.new()
  title("Taxon group legend", cex.main = 1.08, line = 0.2)
  clear_plot_legend(
    "center",
    labels = names(group_colors),
    fill = unname(group_colors),
    cex = legend_text_size(names(group_colors), base = 1.08, minimum = 0.9),
    ncol = separate_legend_column_count(names(group_colors))
  )
}

countTotal <- function(dataset, count_column) {
  counts <- tapply(
    dataset[[count_column]],
    list(dataset$new_periods, dataset$new_area),
    sum,
    na.rm = TRUE
  )
  counts <- as.matrix(counts)
  counts[is.na(counts)] <- 0

  if (nrow(counts) == 0 || ncol(counts) == 0) {
    empty_plot_message("No data available for the selected configuration.")
    return(invisible(NULL))
  }

  area_ids <- colnames(counts)
  area_colors <- palette.colors(max(length(area_ids), 1), "Pastel 1", recycle = TRUE)
  x_positions <- seq_len(nrow(counts))
  y_range <- range(counts, na.rm = TRUE)

  if (!all(is.finite(y_range))) {
    y_range <- c(0, 1)
  }
  if (diff(y_range) == 0) {
    y_range <- c(0, y_range[2] + 1)
  }

  plot(
    x_positions,
    counts[, 1],
    type = "n",
    xlab = "",
    ylab = count_column,
    xaxt = "n",
    ylim = y_range
  )

  for (column_index in seq_len(ncol(counts))) {
    lines(
      x_positions,
      counts[, column_index],
      type = "o",
      lwd = 2.5,
      pch = 16,
      col = area_colors[column_index]
    )
  }

  axis(1, at = x_positions, labels = rownames(counts), las = 2)
  clear_plot_legend(
    "topright",
    labels = paste("Area", area_ids),
    col = area_colors,
    lwd = 2.5,
    pch = 16,
    cex = 1.08
  )
}

bySpeciesComposition <- function(dataset, count_column) {
  area_ids <- sort(unique(na.omit(dataset$new_area)))
  all_groups <- sort(unique(as.character(stats::na.omit(dataset$new_txgroups))))

  if (length(area_ids) == 0 || !length(all_groups)) {
    empty_plot_message("No selected areas contain grouped data.")
    return(invisible(NULL))
  }

  group_colors <- taxon_group_palette(all_groups)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  par(
    mfrow = c(1, length(area_ids)),
    oma = c(2, 0, 0, 5),
    mar = c(8, 4, 3, 1),
    xpd = NA
  )

  for (area_index in seq_along(area_ids)) {
    area_id <- area_ids[area_index]
    region <- dataset[dataset$new_area == area_id & !is.na(dataset$new_txgroups), , drop = FALSE]
    composition <- tapply(
      region[[count_column]],
      list(region$new_txgroups, region$new_periods),
      sum,
      na.rm = TRUE
    )
    composition <- as.matrix(composition)
    composition[is.na(composition)] <- 0
    composition <- composition[rowSums(composition) > 0, , drop = FALSE]

    if (nrow(composition) == 0 || ncol(composition) == 0) {
      empty_plot_message(paste("No grouped taxa for Area", area_id))
      next
    }

    proportions <- sweep(composition, 2, colSums(composition), FUN = "/")
    proportions[is.na(proportions)] <- 0

    barplot(
      proportions,
      col = unname(group_colors[rownames(proportions)]),
      main = paste("Area", area_id),
      names.arg = colnames(proportions),
      las = 2,
      ylab = "Share of counts",
      border = NA
    )

  }
}

byCultureComposition <- function(dataset, count_column) {
  culture_column <- intersect(c("Culture", "Culture1"), names(dataset))[1]

  if (is.na(culture_column) || is.null(culture_column)) {
    empty_plot_message("No culture field is available for the selected dataset.")
    return(invisible(NULL))
  }

  culture_values <- trimws(as.character(dataset[[culture_column]]))
  valid_rows <- !is.na(dataset$new_txgroups) & nzchar(culture_values)
  dataset <- dataset[valid_rows, , drop = FALSE]
  culture_values <- culture_values[valid_rows]

  if (nrow(dataset) == 0) {
    empty_plot_message("No grouped taxa with culture information are available.")
    return(invisible(NULL))
  }

  composition <- tapply(
    dataset[[count_column]],
    list(dataset$new_txgroups, culture_values),
    sum,
    na.rm = TRUE
  )
  composition <- as.matrix(composition)
  composition[is.na(composition)] <- 0
  composition <- composition[rowSums(composition) > 0, , drop = FALSE]
  composition <- composition[, colSums(composition) > 0, drop = FALSE]

  if (nrow(composition) == 0 || ncol(composition) == 0) {
    empty_plot_message("No grouped taxa with culture information are available.")
    return(invisible(NULL))
  }

  proportions <- sweep(composition, 2, colSums(composition), FUN = "/")
  proportions[is.na(proportions)] <- 0
  group_colors <- taxon_group_palette(rownames(proportions))

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  par(mar = c(9, 4, 4, 1), xpd = TRUE)
  barplot(
    proportions,
    col = unname(group_colors[rownames(proportions)]),
    las = 2,
    ylab = "Share of counts",
    main = "Taxon composition by culture",
    border = NA
  )

  invisible(NULL)
}

shorten_segments <- function(coordinates, shorten_length = 0.03) {
  if (nrow(coordinates) < 2) {
    return(NULL)
  }

  start_points <- coordinates[-nrow(coordinates), , drop = FALSE]
  end_points <- coordinates[-1, , drop = FALSE]
  delta <- end_points - start_points
  distance <- sqrt(rowSums(delta^2))
  valid_rows <- distance > 0

  if (!any(valid_rows)) {
    return(NULL)
  }

  start_points <- start_points[valid_rows, , drop = FALSE]
  end_points <- end_points[valid_rows, , drop = FALSE]
  delta <- delta[valid_rows, , drop = FALSE]
  distance <- distance[valid_rows]
  scale <- shorten_length / distance

  cbind(
    x0 = start_points[, 1] + delta[, 1] * scale,
    y0 = start_points[, 2] + delta[, 2] * scale,
    x1 = end_points[, 1] - delta[, 1] * scale,
    y1 = end_points[, 2] - delta[, 2] * scale
  )
}

present_period_levels <- function(period_values, period_levels = NULL) {
  if (!is.null(period_levels) && length(period_levels)) {
    used_periods <- unique(as.character(period_values))
    period_levels <- period_levels[period_levels %in% used_periods]
    if (length(period_levels)) {
      return(period_levels)
    }
  }

  resolve_period_levels(period_values)
}

plotOrdinationArrows <- function(ordination_result, period_levels = NULL, culture_coordinates = NULL) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  par(xpd = TRUE)  
  row_coordinates <- as.matrix(ordination_result$row_coordinates)
  if (is.null(colnames(row_coordinates))) {
    colnames(row_coordinates) <- paste('Dim', seq_len(ncol(row_coordinates)))
  }
  if (ncol(row_coordinates) < 2) {
    empty_plot_message(sprintf("At least two dimensions are needed for the %s map.", tolower(ordination_result$method_label)))
    return(invisible(NULL))
  }

  x_name <- colnames(row_coordinates)[1]
  y_name <- colnames(row_coordinates)[2]
  coordinate_frame <- cbind(
    parse_phase_labels(rownames(row_coordinates)),
    label = rownames(row_coordinates),
    as.data.frame(row_coordinates[, 1:2, drop = FALSE], check.names = FALSE)
  )

  period_levels <- present_period_levels(coordinate_frame$period, period_levels)
  coordinate_frame$period <- factor(coordinate_frame$period, levels = period_levels, ordered = TRUE)
  area_ids <- sort(unique(coordinate_frame$area))
  area_colors <- palette.colors(max(length(area_ids), 1), "Pastel 1", recycle = TRUE)

  plot(
    row_coordinates[, 1],
    row_coordinates[, 2],
    type = "n",
    xlab = x_name,
    ylab = y_name
  )
  abline(v = 0, h = 0, lty = 2, col = "#9ca3af")

  for (area_index in seq_along(area_ids)) {
    area_data <- coordinate_frame[coordinate_frame$area == area_ids[area_index], , drop = FALSE]
    area_data <- area_data[order(area_data$period), , drop = FALSE]

    points( area_data[[x_name]], area_data[[y_name]], pch = 16, cex = 1.2, col = area_colors[area_index])
    text(
      area_data[[x_name]],
      area_data[[y_name]],
      labels = area_data$label,
      pos = 3,
      cex = 0.8,
      col = area_colors[area_index]
    )

    arrow_segments <- shorten_segments(as.matrix(area_data[, c(x_name, y_name)]), 0.03)
    if (!is.null(arrow_segments)) {
      arrows(
        x0 = arrow_segments[, 1],
        y0 = arrow_segments[, 2],
        x1 = arrow_segments[, 3],
        y1 = arrow_segments[, 4],
        col = area_colors[area_index],
        length = 0.08,
        lwd = 2.2
      )
    }
  }

  column_coordinates <- as.matrix(ordination_result$column_coordinates)
  if (!is.null(column_coordinates) && ncol(column_coordinates) >= 2) {
    text(
      column_coordinates[, 1],
      column_coordinates[, 2],
      rownames(column_coordinates),
      col = "darkgreen",
      font = 3
    )
  }

  culture_color <- "#b45309"
  if (!is.null(culture_coordinates)) {
    culture_coordinates <- as.matrix(culture_coordinates)
    if (nrow(culture_coordinates) && ncol(culture_coordinates) >= 2) {
      points(
        culture_coordinates[, 1],
        culture_coordinates[, 2],
        pch = 17,
        cex = 1.25,
        col = culture_color
      )
      text(
        culture_coordinates[, 1],
        culture_coordinates[, 2],
        labels = rownames(culture_coordinates),
        pos = 4,
        cex = 0.78,
        col = culture_color,
        font = 2
      )
    }
  }

  legend_labels <- paste("Area", area_ids)
  legend_colors <- area_colors
  legend_symbols <- rep(16, length(area_ids))
  if (!is.null(culture_coordinates) && nrow(as.matrix(culture_coordinates))) {
    legend_labels <- c(legend_labels, "Culture")
    legend_colors <- c(legend_colors, culture_color)
    legend_symbols <- c(legend_symbols, 17)
  }

  clear_plot_legend(
    "topright",
    labels = legend_labels,
    col = legend_colors,
    pch = legend_symbols,
    lwd = 2.2,
    cex = 1.05
  )
}

plotOrdinationDimensions <- function(ordination_result, period_levels = NULL) {
  row_coordinates <- as.matrix(ordination_result$row_coordinates)
  if (is.null(dim(row_coordinates))) {
    row_coordinates <- matrix(
      row_coordinates,
      ncol = 1,
      dimnames = list(names(ordination_result$row_coordinates), "Dim 1")
    )
  }
  if (is.null(colnames(row_coordinates))) {
    colnames(row_coordinates) <- paste('Dim', seq_len(ncol(row_coordinates)))
  }

  coordinate_frame <- cbind(
    parse_phase_labels(rownames(row_coordinates)),
    as.data.frame(row_coordinates, check.names = FALSE)
  )

  period_levels <- present_period_levels(coordinate_frame$period, period_levels)
  coordinate_frame$period <- factor(coordinate_frame$period, levels = period_levels, ordered = TRUE)
  area_ids <- sort(unique(coordinate_frame$area))
  area_colors <- palette.colors(max(length(area_ids), 1), "Pastel 1", recycle = TRUE)
  dimensions_to_plot <- seq_len(min(2, ncol(row_coordinates)))

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  par(mfrow = c(length(dimensions_to_plot), 1), oma = c(4, 4, 1, 1), mar = c(2, 4, 2, 1))

  for (dimension_index in dimensions_to_plot) {
    dimension_name <- colnames(row_coordinates)[dimension_index]
    y_values <- coordinate_frame[[dimension_name]]
    y_range <- range(y_values, na.rm = TRUE)

    if (!all(is.finite(y_range))) {
      y_range <- c(-1, 1)
    }
    if (diff(y_range) == 0) {
      y_range <- y_range + c(-1, 1)
    }

    plot(
      seq_along(period_levels),
      rep(0, length(period_levels)),
      type = "n",
      xlab = "",
      ylab = dimension_name,
      xaxt = "n",
      xlim = c(1, length(period_levels)),
      ylim = y_range
    )
    axis(1, at = seq_along(period_levels), labels = period_levels, las = 2)

    for (area_index in seq_along(area_ids)) {
      area_data <- coordinate_frame[coordinate_frame$area == area_ids[area_index], , drop = FALSE]
      area_data <- area_data[order(area_data$period), , drop = FALSE]
      x_positions <- match(as.character(area_data$period), period_levels)

      lines(
        x_positions,
        area_data[[dimension_name]],
        type = "o",
        pch = 16,
        lwd = 2.2,
        col = area_colors[area_index]
      )
    }

    clear_plot_legend(
      "topright",
      labels = paste("Area", area_ids),
      col = area_colors,
      lwd = 2.2,
      pch = 16,
      cex = 1.05
    )
  }
}

plotOrdinationDiagnostics <- function(ordination_result, max_dimensions = 10) {
  eig <- ordination_result$eigenvalues
  if (is.null(eig) || !nrow(eig)) {
    empty_plot_message("No ordination diagnostics are available.")
    return(invisible(NULL))
  }

  dimension_count <- min(nrow(eig), max_dimensions)
  percentages <- eig[seq_len(dimension_count), 2]
  dimension_labels <- paste("Dim", seq_len(dimension_count))
  y_max <- max(percentages, na.rm = TRUE)
  if (!is.finite(y_max) || y_max <= 0) {
    y_max <- 1
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  par(mar = c(6, 4, 3, 1))
  mids <- barplot(
    percentages,
    names.arg = dimension_labels,
    las = 2,
    ylab = ordination_result$diagnostic_label,
    main = sprintf("%s diagnostics", ordination_result$method_label),
    ylim = c(0, y_max * 1.15),
    col = "#9ecae1",
    border = "#4b5563"
  )

  text(
    x = mids,
    y = percentages,
    labels = sprintf("%.1f%%", percentages),
    pos = 3,
    cex = 0.85
  )
}
