empty_plot_message <- function(message) {
  plot.new()
  text(0.5, 0.5, message, cex = 1.05)
}

taxon_group_palette <- function(group_names) {
  group_names <- sort(unique(as.character(group_names)))
  colors <- palette.colors(max(length(group_names), 1), "Dark 2", recycle = TRUE)
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
  shorten <- pmin(shorten_length, distance * 0.42)
  scale <- shorten / distance

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

plotOrdinationArrowsBase <- function(ordination_result, period_levels = NULL, culture_coordinates = NULL) {
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
  column_coordinates <- if (is.null(ordination_result$column_coordinates)) NULL else as.matrix(ordination_result$column_coordinates)
  culture_matrix <- if (is.null(culture_coordinates)) NULL else as.matrix(culture_coordinates)

  plot_x <- row_coordinates[, 1]
  plot_y <- row_coordinates[, 2]
  if (!is.null(column_coordinates) && nrow(column_coordinates) && ncol(column_coordinates) >= 2) {
    plot_x <- c(plot_x, column_coordinates[, 1])
    plot_y <- c(plot_y, column_coordinates[, 2])
  }
  if (!is.null(culture_matrix) && nrow(culture_matrix) && ncol(culture_matrix) >= 2) {
    plot_x <- c(plot_x, culture_matrix[, 1])
    plot_y <- c(plot_y, culture_matrix[, 2])
  }

  plot_x <- plot_x[is.finite(plot_x)]
  plot_y <- plot_y[is.finite(plot_y)]
  x_range <- range(plot_x, na.rm = TRUE)
  y_range <- range(plot_y, na.rm = TRUE)
  if (!all(is.finite(x_range)) || diff(x_range) == 0) {
    x_range <- x_range + c(-1, 1)
  }
  if (!all(is.finite(y_range)) || diff(y_range) == 0) {
    y_range <- y_range + c(-1, 1)
  }

  x_padding <- diff(x_range) * 0.22
  y_padding <- diff(y_range) * 0.20
  x_limits <- x_range + c(-x_padding, x_padding * 1.45)
  y_limits <- y_range + c(-y_padding, y_padding * 1.35)

  plot(
    row_coordinates[, 1],
    row_coordinates[, 2],
    type = "n",
    xlab = x_name,
    ylab = y_name,
    xlim = x_limits,
    ylim = y_limits
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

  if (!is.null(column_coordinates) && ncol(column_coordinates) >= 2) {
    text(
      column_coordinates[, 1],
      column_coordinates[, 2],
      rownames(column_coordinates),
      pos = 4,
      col = "#111827",
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

plotOrdinationArrows <- function(
  ordination_result,
  period_levels = NULL,
  culture_coordinates = NULL,
  show_year_labels = FALSE,
  show_taxa_labels = FALSE,
  show_culture_labels = FALSE
) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("ggrepel", quietly = TRUE)) {
    plotOrdinationArrowsBase(ordination_result, period_levels, culture_coordinates)
    return(invisible(NULL))
  }

  row_coordinates <- as.matrix(ordination_result$row_coordinates)
  if (is.null(colnames(row_coordinates))) {
    colnames(row_coordinates) <- paste("Dim", seq_len(ncol(row_coordinates)))
  }
  if (ncol(row_coordinates) < 2) {
    empty_plot_message(sprintf("At least two dimensions are needed for the %s map.", tolower(ordination_result$method_label)))
    return(invisible(NULL))
  }

  x_name <- colnames(row_coordinates)[1]
  y_name <- colnames(row_coordinates)[2]
  phase_labels <- parse_phase_labels(rownames(row_coordinates))
  coordinate_frame <- data.frame(
    phase_labels,
    label = paste(phase_labels$period, phase_labels$area, sep = "-"),
    x = row_coordinates[, 1],
    y = row_coordinates[, 2],
    stringsAsFactors = FALSE
  )
  coordinate_frame <- coordinate_frame[
    is.finite(coordinate_frame$x) & is.finite(coordinate_frame$y),
    ,
    drop = FALSE
  ]
  if (!nrow(coordinate_frame)) {
    empty_plot_message("No finite ordination coordinates are available for the map.")
    return(invisible(NULL))
  }

  period_levels <- present_period_levels(coordinate_frame$period, period_levels)
  coordinate_frame$period <- factor(coordinate_frame$period, levels = period_levels, ordered = TRUE)
  area_ids <- sort(unique(coordinate_frame$area))
  area_keys <- paste("Area", area_ids)
  area_colors <- palette.colors(max(length(area_ids), 1), "Pastel 1", recycle = TRUE)
  names(area_colors) <- area_keys
  coordinate_frame$legend_key <- paste("Area", coordinate_frame$area)

  column_frame <- data.frame()
  column_coordinates <- if (is.null(ordination_result$column_coordinates)) NULL else as.matrix(ordination_result$column_coordinates)
  if (isTRUE(show_taxa_labels) && !is.null(column_coordinates) && nrow(column_coordinates) && ncol(column_coordinates) >= 2) {
    column_labels <- rownames(column_coordinates)
    if (is.null(column_labels)) {
      column_labels <- paste("Taxon group", seq_len(nrow(column_coordinates)))
    }
    column_frame <- data.frame(
      label = column_labels,
      x = column_coordinates[, 1],
      y = column_coordinates[, 2],
      legend_key = "Taxon group",
      stringsAsFactors = FALSE
    )
    column_frame <- column_frame[
      is.finite(column_frame$x) & is.finite(column_frame$y),
      ,
      drop = FALSE
    ]
  }

  culture_frame <- data.frame()
  culture_matrix <- if (is.null(culture_coordinates)) NULL else as.matrix(culture_coordinates)
  if (isTRUE(show_culture_labels) && !is.null(culture_matrix) && nrow(culture_matrix) && ncol(culture_matrix) >= 2) {
    culture_labels <- rownames(culture_matrix)
    if (is.null(culture_labels)) {
      culture_labels <- paste("Culture", seq_len(nrow(culture_matrix)))
    }
    culture_frame <- data.frame(
      label = culture_labels,
      x = culture_matrix[, 1],
      y = culture_matrix[, 2],
      legend_key = "Culture",
      stringsAsFactors = FALSE
    )
    culture_frame <- culture_frame[
      is.finite(culture_frame$x) & is.finite(culture_frame$y),
      ,
      drop = FALSE
    ]
  }

  all_x <- c(coordinate_frame$x, column_frame$x, culture_frame$x)
  all_y <- c(coordinate_frame$y, column_frame$y, culture_frame$y)
  x_range <- range(all_x, na.rm = TRUE)
  y_range <- range(all_y, na.rm = TRUE)
  if (!all(is.finite(x_range)) || diff(x_range) == 0) {
    x_range <- x_range + c(-1, 1)
  }
  if (!all(is.finite(y_range)) || diff(y_range) == 0) {
    y_range <- y_range + c(-1, 1)
  }
  map_diagonal <- sqrt(diff(x_range)^2 + diff(y_range)^2)
  arrow_shorten <- if (is.finite(map_diagonal) && map_diagonal > 0) map_diagonal * 0.016 else 0.03

  arrow_frames <- lapply(area_ids, function(area_id) {
    area_data <- coordinate_frame[coordinate_frame$area == area_id, , drop = FALSE]
    area_data <- area_data[order(area_data$period), , drop = FALSE]
    arrow_segments <- shorten_segments(as.matrix(area_data[, c("x", "y")]), arrow_shorten)
    if (is.null(arrow_segments)) {
      return(NULL)
    }

    data.frame(
      arrow_segments,
      legend_key = paste("Area", area_id),
      stringsAsFactors = FALSE
    )
  })
  arrow_frame <- do.call(rbind, arrow_frames[!vapply(arrow_frames, is.null, logical(1))])
  if (is.null(arrow_frame)) {
    arrow_frame <- data.frame()
  }

  label_frame <- rbind(
    if (isTRUE(show_year_labels)) {
      data.frame(
        label = coordinate_frame$label,
        x = coordinate_frame$x,
        y = coordinate_frame$y,
        legend_key = coordinate_frame$legend_key,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame()
    },
    column_frame,
    culture_frame
  )
  label_frame$label <- as.character(label_frame$label)
  label_frame <- label_frame[nzchar(label_frame$label), , drop = FALSE]
  label_frame$label_fill_key <- label_frame$legend_key

  label_fill_alpha <- 0.66
  color_values <- area_colors
  fill_values <- grDevices::adjustcolor(area_colors, alpha.f = label_fill_alpha)
  names(fill_values) <- names(area_colors)
  legend_breaks <- area_keys
  if (nrow(column_frame)) {
    color_values <- c(color_values, "Taxon group" = "#111827")
    fill_values <- c(
      fill_values,
      "Taxon group" = grDevices::adjustcolor("#f8fafc", alpha.f = 0.76)
    )
    legend_breaks <- c(legend_breaks, "Taxon group")
  }
  if (nrow(culture_frame)) {
    color_values <- c(color_values, "Culture" = "#b45309")
    fill_values <- c(
      fill_values,
      "Culture" = grDevices::adjustcolor("#ffedd5", alpha.f = 0.76)
    )
    legend_breaks <- c(legend_breaks, "Culture")
  }

  plot <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "#cbd5e1", linewidth = 0.35) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "#cbd5e1", linewidth = 0.35)

  if (nrow(arrow_frame)) {
    plot <- plot +
      ggplot2::geom_segment(
        data = arrow_frame,
        ggplot2::aes(x = x0, y = y0, xend = x1, yend = y1, color = legend_key),
        arrow = grid::arrow(type = "closed", length = grid::unit(0.11, "inches")),
        lineend = "round",
        linewidth = 0.8,
        alpha = 0.88
      )
  }

  if (isTRUE(show_year_labels)) {
    plot <- plot +
      ggplot2::geom_point(
        data = coordinate_frame,
        ggplot2::aes(x = x, y = y, color = legend_key),
        shape = 16,
        size = 2.8,
        alpha = 0.95,
        show.legend = FALSE
      )
  }

  if (nrow(column_frame)) {
    plot <- plot +
      ggplot2::geom_point(
        data = column_frame,
        ggplot2::aes(x = x, y = y, color = legend_key),
        shape = 4,
        size = 2.2,
        stroke = 0.8
      )
  }

  if (nrow(culture_frame)) {
    plot <- plot +
      ggplot2::geom_point(
        data = culture_frame,
        ggplot2::aes(x = x, y = y, color = legend_key),
        shape = 17,
        size = 3.1
      )
  }

  if (nrow(label_frame)) {
    plot <- plot +
      ggrepel::geom_label_repel(
        data = label_frame,
        ggplot2::aes(x = x, y = y, label = label, fill = label_fill_key),
        color = "#111827",
        label.size = 0.14,
        label.padding = 0.14,
        label.r = 0.06,
        box.padding = 0.36,
        point.padding = 0.20,
        min.segment.length = 0,
        max.overlaps = Inf,
        max.iter = 30000,
        max.time = 2.5,
        force = 2.3,
        force_pull = 0.16,
        seed = 42,
        segment.color = "#475569",
        segment.alpha = 0.62,
        segment.linetype = "dotted",
        segment.size = 0.25,
        size = 3.25,
        show.legend = FALSE
      )
  }

  plot <- plot +
    ggplot2::scale_color_manual(
      values = color_values,
      breaks = legend_breaks,
      name = NULL,
      drop = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = fill_values,
      guide = "none",
      drop = FALSE
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.18, 0.30))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.18, 0.28))) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = x_name, y = y_name) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        nrow = 2,
        byrow = TRUE,
        override.aes = list(size = 3, alpha = 1)
      )
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.key.height = grid::unit(0.16, "inches"),
      legend.key.width = grid::unit(0.26, "inches"),
      legend.text = ggplot2::element_text(size = 9),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "#e5e7eb", linewidth = 0.35),
      axis.title = ggplot2::element_text(color = "#111827"),
      axis.text = ggplot2::element_text(color = "#374151"),
      plot.margin = ggplot2::margin(12, 18, 10, 12)
    )

  print(plot)
  invisible(plot)
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
