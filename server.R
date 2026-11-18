# Server -----------------------------------------------------------------------

server <- function(input, output, session) {
  # State ---------------------------------------------------------------------
  data_bundle <- reactiveVal(NULL)
  data_loading <- reactiveVal(TRUE)
  data_error <- reactiveVal(NULL)
  analysis_result <- reactiveVal(NULL)
  analysis_stale <- reactiveVal(FALSE)
  current_polygons <- reactiveVal(NULL)
  default_polygons_attempted <- reactiveVal(FALSE)
  grouping_revision <- reactiveVal(0)
  group_catalog_version <- reactiveVal(0)
  reference_csv <- reactiveVal(NULL)
  reference_points <- reactiveVal(NULL)
  reference_info <- reactiveVal(NULL)
  reference_signature <- reactiveVal(NULL)

  group_edit_store <- new.env(parent = emptyenv())

  # Prefer the current study-area file even if the deployment still has an
  # older global.R defining DEFAULT_POLYGON_PATH.
  default_polygon_path <- "groupings/spatial/polygon_18112026_SC_v2.gpkg"
  if (!file.exists(default_polygon_path) && exists("DEFAULT_POLYGON_PATH", inherits = TRUE)) {
    default_polygon_path <- DEFAULT_POLYGON_PATH
  }

  # Data loading --------------------------------------------------------------
  load_data <- function(refresh = FALSE) {
    data_loading(TRUE)
    data_error(NULL)

    on.exit(data_loading(FALSE), add = TRUE)

    shiny::withProgress(
      message = if (refresh) "Refreshing BIAD data" else "Loading BIAD data",
      value = 0,
      {
        shiny::setProgress(0.15, detail = "Reading cache or querying the database")
        bundle <- tryCatch(
          load_biad_data(refresh = refresh),
          error = function(error) error
        )
        shiny::setProgress(0.85, detail = "Preparing the interface")

        if (inherits(bundle, "error")) {
          data_error(conditionMessage(bundle))
          data_bundle(NULL)
          return(invisible(NULL))
        }

        data_bundle(bundle)
        data_error(NULL)

        # A new data version invalidates the spatial and analytical state.
        current_polygons(NULL)
        default_polygons_attempted(FALSE)
        analysis_result(NULL)
        analysis_stale(FALSE)
      }
    )

    invisible(NULL)
  }

  observe({
    load_data(refresh = FALSE)
  })

  observeEvent(input$refresh_data, {
    load_data(refresh = TRUE)
    if (is.null(data_error())) {
      showNotification("BIAD data refreshed.", type = "message", duration = 4)
    } else {
      showNotification(data_error(), type = "error", duration = 10)
    }
  }, ignoreInit = TRUE)

  bundle <- reactive({
    shiny::req(data_bundle())
  })

  # Data status ---------------------------------------------------------------
  output$data_status <- renderUI({
    if (data_loading()) {
      return(div(
        class = "data-status data-status--loading",
        tags$span(class = "data-status-dot"),
        "Loading BIAD data — this can take a moment on the first run."
      ))
    }

    if (!is.null(data_error())) {
      return(div(
        class = "data-status data-status--error",
        app_icon("triangle-exclamation"),
        tags$span(data_error())
      ))
    }

    if (is.null(data_bundle())) return(NULL)

    summary <- biad_data_summary(data_bundle())
    source_label <- if (identical(summary$source, "cache")) {
      paste("cache", format_cache_age(summary$cache_age))
    } else {
      "live database"
    }

    div(
      class = "data-status",
      tags$span(class = "data-status-item", app_icon("circle-check"), "Data ready"),
      tags$span(class = "data-status-item", app_icon("hard-drive"), source_label),
      tags$span(class = "data-status-item", app_icon("paw"), sprintf("%s faunal records", format(summary$faunal_records, big.mark = ","))),
      tags$span(class = "data-status-item", app_icon("seedling"), sprintf("%s botanical records", format(summary$botanical_records, big.mark = ","))),
      tags$span(class = "data-status-item", app_icon("location-dot"), sprintf("%s faunal sites", format(summary$faunal_sites, big.mark = ","))),
      tags$span(class = "data-status-item", app_icon("location-dot"), sprintf("%s botanical sites", format(summary$botanical_sites, big.mark = ",")))
    )
  })

  # Dataset and site selection ------------------------------------------------
  selected_dataset <- reactive({
    data <- bundle()
    switch(
      input$data_type_selector,
      "Botanical" = data$botanical,
      "Combined" = combine_taxa_datasets(data$faunal, data$botanical),
      data$faunal
    )
  })

  selected_sites <- reactive({
    data <- data_bundle()
    if (is.null(data)) return(NULL)
    switch(
      input$data_type_selector,
      "Botanical" = data$botanical_sites,
      "Combined" = combine_site_points(data$faunal_sites, data$botanical_sites),
      data$faunal_sites
    )
  })

  record_count <- reactive({
    data <- data_bundle()
    if (is.null(data)) return(0L)
    switch(
      input$data_type_selector,
      "Botanical" = nrow(data$botanical),
      "Combined" = nrow(data$faunal) + nrow(data$botanical),
      nrow(data$faunal)
    )
  })

  # Grouping selection ---------------------------------------------------------
  active_grouping_path <- reactive({
    if (identical(input$data_type_selector, "Combined")) {
      input$file_selector_tx_combined_faunal
    } else {
      input$file_selector_tx
    }
  })

  refresh_grouping_inputs <- function(selected_path = NULL) {
    if (identical(input$data_type_selector, "Combined")) {
      choices <- available_groupings("Faunal")
      current <- selected_path %||% isolate(input$file_selector_tx_combined_faunal)
      if (is.null(current) || !current %in% unname(choices)) current <- preferred_grouping_path("Faunal")
      updateSelectInput(session, "file_selector_tx_combined_faunal", choices = choices, selected = current)
    } else {
      choices <- available_groupings(input$data_type_selector)
      current <- selected_path %||% isolate(input$file_selector_tx)
      if (is.null(current) || !current %in% unname(choices)) current <- preferred_grouping_path(input$data_type_selector)
      updateSelectInput(session, "file_selector_tx", choices = choices, selected = current)
    }
  }

  observeEvent(input$data_type_selector, {
    if (identical(input$data_type_selector, "Combined")) {
      choices <- available_groupings("Faunal")
      current <- isolate(input$file_selector_tx_combined_faunal)
      if (is.null(current) || !current %in% unname(choices)) current <- preferred_grouping_path("Faunal")
      updateSelectInput(session, "file_selector_tx_combined_faunal", choices = choices, selected = current)
    } else {
      choices <- available_groupings(input$data_type_selector)
      current <- isolate(input$file_selector_tx)
      if (is.null(current) || !current %in% unname(choices)) current <- preferred_grouping_path(input$data_type_selector)
      updateSelectInput(session, "file_selector_tx", choices = choices, selected = current)
    }
  }, ignoreInit = TRUE)

  output$combined_botanical_groups_ui <- renderUI({
    group_catalog_version()
    botanical_choices <- group_catalog_choices("Botanical")
    selected_groups <- isolate(input$combined_botanical_groups)
    if (is.null(selected_groups)) selected_groups <- character(0)

    catalog <- group_catalog("Botanical")
    selected_groups <- intersect(selected_groups, catalog$option_id)

    selectizeInput(
      "combined_botanical_groups",
      "Botanical subgroups to include",
      choices = botanical_choices,
      selected = selected_groups,
      multiple = TRUE,
      width = "100%",
      options = list(
        plugins = list("remove_button"),
        placeholder = "Choose one or more botanical subgroups"
      )
    )
  })

  selected_grouping_label <- reactive({
    if (!identical(input$data_type_selector, "Combined")) {
      if (is.null(input$file_selector_tx)) return("—")
      return(format_grouping_label(input$file_selector_tx))
    }

    faunal_label <- if (is.null(input$file_selector_tx_combined_faunal)) {
      "—"
    } else {
      format_grouping_label(input$file_selector_tx_combined_faunal)
    }

    botanical_groups <- input$combined_botanical_groups
    if (is.null(botanical_groups)) botanical_groups <- character(0)

    botanical_catalog <- group_catalog("Botanical")
    selected_botanical <- botanical_catalog[match(botanical_groups, botanical_catalog$option_id), , drop = FALSE]
    selected_botanical <- selected_botanical[!is.na(selected_botanical$option_id), , drop = FALSE]

    botanical_suffix <- if (!length(botanical_groups)) {
      "no botanical subgroups selected"
    } else if (length(botanical_groups) <= 2) {
      paste(selected_botanical$output_label, collapse = ", ")
    } else {
      sprintf("%s botanical subgroups", length(botanical_groups))
    }

    sprintf("Faunal: %s | Botanical: %s", faunal_label, botanical_suffix)
  })

  grouped_dataset <- reactive({
    grouping_revision()
    data <- bundle()

    if (identical(input$data_type_selector, "Combined")) {
      shiny::req(input$file_selector_tx_combined_faunal)
      shiny::req(input$combined_botanical_groups)
      return(combine_grouped_taxa_datasets(
        faunal_dataset = data$faunal,
        botanical_dataset = data$botanical,
        faunal_grouping_file = input$file_selector_tx_combined_faunal,
        botanical_group_ids = input$combined_botanical_groups
      ))
    }

    shiny::req(input$file_selector_tx)
    groupTaxons(selected_dataset(), input$file_selector_tx, logging = FALSE)
  })

  prepared_dataset <- reactive({
    dataset <- grouped_dataset()
    notes <- attr(dataset, "analysis_notes", exact = TRUE)
    dataset$new_periods <- split_periods(
      dataset = dataset,
      selection = input$file_selector_per,
      start = input$start_value,
      end = input$end_value,
      len = input$duration
    )
    if (!is.null(notes)) attr(dataset, "analysis_notes") <- notes
    dataset
  })

  # Staleness -----------------------------------------------------------------
  observe({
    data_bundle()
    input$data_type_selector
    input$file_selector_tx
    input$file_selector_tx_combined_faunal
    input$combined_botanical_groups
    input$file_selector_per
    input$start_value
    input$end_value
    input$duration
    input$analysis_method
    input$pca_transform
    input$pca_scaling
    current_polygons()

    if (!is.null(isolate(analysis_result()))) {
      analysis_stale(TRUE)
    }
  })

  # Sidebar outputs -----------------------------------------------------------
  output$taxon_table <- renderUI({
    if (data_loading()) return(helper_text("Loading data ..."))
    if (is.null(data_bundle())) return(NULL)

    if (identical(input$data_type_selector, "Combined")) {
      if (is.null(input$combined_botanical_groups) || !length(input$combined_botanical_groups)) {
        return(helper_text("Choose at least one botanical subgroup to build the combined dataset."))
      }
    } else if (is.null(input$file_selector_tx) || !file.exists(input$file_selector_tx)) {
      return(helper_text("Choose a grouping file to preview its composition."))
    }

    stats <- tryCatch(group_summary_stats(grouped_dataset()), error = function(error) NULL)
    if (is.null(stats)) return(helper_text("Select a valid grouping file to preview its composition."))

    build_group_summary_ui(stats)
  })

  output$selection_summary <- renderUI({
    if (data_loading()) {
      return(build_selection_summary_ui(
        data_type = input$data_type_selector,
        site_count = 0,
        record_count = 0,
        polygon_count = 0,
        grouping_label = "—",
        stale = FALSE,
        has_results = FALSE,
        loading = TRUE
      ))
    }

    if (is.null(data_bundle())) {
      return(empty_state(
        "Data unavailable",
        data_error() %||% "The BIAD data could not be loaded.",
        "triangle-exclamation"
      ))
    }

    points <- selected_sites()
    polygons <- current_polygons()
    total_sites <- if (is.null(points)) 0L else nrow(points)

    site_count_label <- NULL
    if (isTRUE(input$only_selected_sites)) {
      visible <- map_site_layer()$points
      visible_sites <- if (is.null(visible)) 0L else nrow(visible)
      site_count_label <- sprintf(
        "%s / %s",
        format(visible_sites, big.mark = ","),
        format(total_sites, big.mark = ",")
      )
    }

    build_selection_summary_ui(
      data_type = input$data_type_selector,
      site_count = total_sites,
      record_count = record_count(),
      polygon_count = if (is.null(polygons)) 0L else nrow(polygons),
      grouping_label = selected_grouping_label(),
      stale = analysis_stale(),
      has_results = !is.null(analysis_result()),
      site_count_label = site_count_label
    )
  })

  output$time_validation <- renderUI({
    if (!identical(input$file_selector_per, "auto")) return(NULL)

    start <- input$start_value
    end <- input$end_value
    duration <- input$duration
    if (is.null(start) || is.null(end) || is.null(duration)) return(NULL)

    if (start <= end) {
      return(div(
        class = "inline-alert inline-alert--error",
        app_icon("triangle-exclamation"),
        "Start (older bound) must be greater than End (younger bound) for BP slices."
      ))
    }
    if (duration <= 0) {
      return(div(
        class = "inline-alert inline-alert--error",
        app_icon("triangle-exclamation"),
        "Slice width must be greater than zero."
      ))
    }

    slices <- max(1, ceiling((start - end) / duration))
    div(
      class = "inline-alert inline-alert--info",
      app_icon("circle-info"),
      sprintf("About %s time slices of up to %s years.", slices, duration)
    )
  })

  output$method_note <- renderUI({
    if (identical(input$analysis_method, "pca")) {
      helper_text(
        "PCA is a linear ordination. Log-transforming counts and scaling columns ",
        "usually makes the structure across rare and common groups easier to compare."
      )
    } else {
      helper_text(
        "Correspondence analysis is designed for count tables and gives equal ",
        "weight to rare and common taxon groups through the chi-square metric."
      )
    }
  })

  output$step1_hint <- renderUI({
    if (data_loading()) return("Loading data")
    if (is.null(data_bundle())) return("Data unavailable")

    if (identical(input$data_type_selector, "Combined")) {
      count <- if (is.null(input$combined_botanical_groups)) 0L else length(input$combined_botanical_groups)
      if (count == 0) return("Combined · choose botanical subgroups")
      return(sprintf("Combined · %s botanical subgroup%s", count, if (count == 1) "" else "s"))
    }

    grouping <- if (is.null(input$file_selector_tx)) {
      "no grouping"
    } else {
      format_grouping_label(input$file_selector_tx)
    }
    sprintf("%s · %s", input$data_type_selector, grouping)
  })

  output$step2_hint <- renderUI({
    if (identical(input$file_selector_per, "auto")) {
      start <- input$start_value
      end <- input$end_value
      if (is.null(start) || is.null(end)) return("Automatic BP slices")
      return(sprintf("Automatic · %s–%s BP", start, end))
    }
    "Named archaeological periods"
  })

  output$step3_hint <- renderUI({
    if (identical(input$analysis_method, "pca")) {
      transform <- if (identical(input$pca_transform, "log1p")) "log(count+1)" else "raw counts"
      scaling <- switch(
        input$pca_scaling,
        "scale" = "center + scale",
        "center" = "center only",
        "none" = "no scaling",
        "scaled"
      )
      return(sprintf("PCA · %s, %s", transform, scaling))
    }
    "Correspondence analysis"
  })

  output$step4_hint <- renderUI({
    if (data_loading()) return("Loading data")
    if (is.null(data_bundle())) return("Data unavailable")

    polygons <- current_polygons()
    count <- if (is.null(polygons)) 0L else nrow(polygons)
    if (count == 0) return("Add a polygon")
    if (!is.null(analysis_result()) && analysis_stale()) return("Inputs changed")
    if (!is.null(analysis_result())) return("Results up to date")
    sprintf("Ready · %s polygon%s", count, if (count == 1) "" else "s")
  })

  output$run_hint <- renderUI({
    if (data_loading()) return(helper_text("Loading data ..."))
    if (is.null(data_bundle())) return(status_pill("Data unavailable", "warning"))

    polygons <- current_polygons()
    if (is.null(polygons) || nrow(polygons) == 0) {
      return(helper_text("Add at least one polygon to enable analysis."))
    }
    if (!is.null(analysis_result()) && analysis_stale()) {
      return(status_pill("Inputs changed — rerun when ready.", "stale"))
    }
    status_pill("Ready to run.", "neutral")
  })

  # Reference points (CSV) ----------------------------------------------------
  output$reference_meta <- renderUI({
    info <- reference_info()
    if (is.null(info) || is.null(info$count) || info$count == 0) return("none")
    sprintf("%s point%s", info$count, if (info$count == 1) "" else "s")
  })

  output$reference_columns_ui <- renderUI({
    state <- reference_csv()
    if (is.null(state)) return(NULL)

    columns <- names(state$data)
    choices <- c("Select a column" = "", stats::setNames(columns, columns))

    tagList(
      selectInput(
        "reference_lon",
        "Longitude column",
        choices = choices,
        selected = state$lon %||% "",
        width = "100%"
      ),
      selectInput(
        "reference_lat",
        "Latitude column",
        choices = choices,
        selected = state$lat %||% "",
        width = "100%"
      )
    )
  })

  output$reference_status_ui <- renderUI({
    state <- reference_csv()
    if (is.null(state)) return(NULL)

    info <- reference_info()
    if (is.null(info)) return(helper_text("Choose the longitude and latitude columns."))

    if (is.null(info$count) || info$count == 0) {
      return(helper_text(info$message %||% "No valid points yet."))
    }

    detail <- sprintf("%s of %s rows have valid coordinates", info$count, info$total)
    if (info$skipped > 0) {
      detail <- paste0(detail, sprintf(" · %s skipped", info$skipped))
    }
    helper_text(detail)
  })

  output$reference_controls_ui <- renderUI({
    info <- reference_info()
    if (is.null(info) || is.null(info$count) || info$count == 0) return(NULL)

    div(
      class = "reference-controls",
      checkboxInput("show_reference_points", "Show on map", TRUE),
      actionButton(
        "clear_reference_points",
        tagList(app_icon("trash"), "Remove"),
        class = "btn-sm btn-danger-soft"
      )
    )
  })

  observeEvent(input$reference_csv, {
    upload <- input$reference_csv
    req(upload$datapath)

    tryCatch({
      data <- utils::read.csv(
        upload$datapath,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        na.strings = c("", "NA", "NaN", "null")
      )

      if (!ncol(data) || !nrow(data)) {
        stop("The CSV is empty.")
      }

      names(data)[1] <- sub("^\ufeff", "", names(data)[1])
      names(data) <- make.unique(names(data))

      detected <- detect_reference_columns(names(data))
      state <- list(
        data = data,
        name = basename(upload$name),
        lon = detected$lon,
        lat = detected$lat
      )
      reference_csv(state)

      info <- build_reference_points(data, detected$lon, detected$lat)
      reference_points(info$points)
      reference_info(list(
        count = if (is.null(info$points)) 0L else nrow(info$points),
        skipped = info$skipped,
        total = info$total,
        message = info$message
      ))
      reference_signature(paste(
        state$name,
        nrow(data),
        detected$lon %||% "",
        detected$lat %||% "",
        sep = "|"
      ))

      if (!is.null(info$points)) {
        msg <- sprintf(
          "Loaded %s: %s point%s",
          state$name,
          nrow(info$points),
          if (nrow(info$points) == 1) "" else "s"
        )
        if (info$skipped > 0) {
          msg <- paste0(msg, sprintf(" (%s skipped)", info$skipped))
        }
        if (!is.null(detected$lon) && !is.null(detected$lat)) {
          msg <- paste0(msg, sprintf(". Using %s / %s.", detected$lon, detected$lat))
        }
        showNotification(msg, type = "message", duration = 6)
      } else {
        showNotification(
          info$message %||% "Choose the longitude and latitude columns.",
          type = "warning",
          duration = 6
        )
      }
    }, error = function(error) {
      showNotification(
        paste("Could not read CSV:", conditionMessage(error)),
        type = "error",
        duration = 8
      )
    })
  }, ignoreInit = TRUE)

  observeEvent(list(input$reference_lon, input$reference_lat), {
    state <- reference_csv()
    if (is.null(state)) return()

    if (!is.null(input$reference_lon) && nzchar(input$reference_lon) &&
        !input$reference_lon %in% names(state$data)) {
      return()
    }
    if (!is.null(input$reference_lat) && nzchar(input$reference_lat) &&
        !input$reference_lat %in% names(state$data)) {
      return()
    }

    info <- build_reference_points(state$data, input$reference_lon, input$reference_lat)
    reference_points(info$points)
    reference_info(list(
      count = if (is.null(info$points)) 0L else nrow(info$points),
      skipped = info$skipped,
      total = info$total,
      message = info$message
    ))

    signature <- paste(
      state$name,
      nrow(state$data),
      input$reference_lon %||% "",
      input$reference_lat %||% "",
      sep = "|"
    )
    if (identical(signature, reference_signature())) return()
    reference_signature(signature)

    if (!is.null(info$points)) {
      showNotification(
        sprintf(
          "%s reference point%s loaded.",
          nrow(info$points),
          if (nrow(info$points) == 1) "" else "s"
        ),
        type = "message",
        duration = 3
      )
      if (info$skipped > 0) {
        showNotification(
          sprintf(
            "%s row%s skipped because coordinates were missing or out of range.",
            info$skipped,
            if (info$skipped == 1) "" else "s"
          ),
          type = "warning",
          duration = 6
        )
      }
    } else if (!is.null(info$message)) {
      showNotification(info$message, type = "warning", duration = 5)
    }
  }, ignoreInit = FALSE)

  observeEvent(input$clear_reference_points, {
    reference_csv(NULL)
    reference_points(NULL)
    reference_info(NULL)
    reference_signature(NULL)
    leafletProxy("map") |> clearGroup("reference-points")
    showNotification("Reference points removed.", type = "message", duration = 3)
  }, ignoreInit = TRUE)

  output$analysis_status <- renderUI({
    build_analysis_status_ui(analysis_result(), analysis_stale())
  })

  output$analysis_overview <- renderUI({
    build_overview_ui(analysis_result(), analysis_stale())
  })

  # Map -----------------------------------------------------------------------
  map_draw_toolbar <- function(proxy, clear_features = FALSE) {
    if (clear_features) proxy <- leaflet.extras::removeDrawToolbar(proxy, clearFeatures = TRUE)
    leaflet.extras::addDrawToolbar(
      proxy,
      targetGroup = "selected-polygons",
      polygonOptions = leaflet.extras::drawPolygonOptions(showArea = FALSE),
      editOptions = leaflet.extras::editToolbarOptions(edit = FALSE, remove = FALSE),
      polylineOptions = FALSE,
      rectangleOptions = FALSE,
      circleOptions = FALSE,
      markerOptions = FALSE,
      circleMarkerOptions = FALSE
    )
  }

  output$map <- renderLeaflet({
    leaflet(options = leafletOptions(preferCanvas = TRUE)) |>
      addProviderTiles(providers$Esri.WorldGrayCanvas) |>
      addScaleBar(position = "bottomleft") |>
      map_draw_toolbar()
  })

  output$map_exports <- renderUI({
    polygons <- current_polygons()
    has_polygons <- !is.null(polygons) && nrow(polygons) > 0
    has_results <- !is.null(analysis_result())

    div(
      class = "export-actions",
      if (has_polygons) {
        downloadButton(
          "download_polygons",
          tagList(app_icon("download"), "Polygons"),
          class = "btn-sm btn-soft"
        )
      },
      if (has_results) {
        downloadButton(
          "download_raw_table",
          tagList(app_icon("download"), "Raw phases"),
          class = "btn-sm btn-soft"
        )
      },
      if (has_results) {
        downloadButton(
          "download_table",
          tagList(app_icon("download"), "Aggregated"),
          class = "btn-sm btn-soft"
        )
      }
    )
  })

  normalized_polygons <- reactive({
    polygons <- current_polygons()
    if (is.null(polygons) || !nrow(polygons)) return(NULL)

    points <- selected_sites()
    if (is.null(points) || !nrow(points)) return(NULL)

    tryCatch(
      normalize_polygon_data(polygons, points),
      error = function(error) NULL
    )
  })

  map_site_layer <- reactive({
    points <- selected_sites()
    if (is.null(points) || !nrow(points)) {
      return(list(points = NULL, colors = character(0)))
    }

    polygons <- normalized_polygons()
    colors <- site_colors(points, polygons, areapal, normalize = FALSE)

    if (isTRUE(input$only_selected_sites)) {
      ids <- site_polygon_ids(points, polygons, normalize = FALSE)
      keep <- !is.na(ids)
      return(list(points = points[keep, , drop = FALSE], colors = colors[keep]))
    }

    list(points = points, colors = colors)
  })

  # Default polygons ----------------------------------------------------------
  observe({
    if (default_polygons_attempted()) return()

    points <- selected_sites()
    if (is.null(points) || !nrow(points)) return()

    default_polygons_attempted(TRUE)
    default_polygons <- tryCatch(
      load_polygon_file(default_polygon_path, points),
      error = function(error) {
        warning(sprintf(
          "Failed to load default polygons from '%s': %s",
          default_polygon_path,
          conditionMessage(error)
        ))
        NULL
      }
    )

    if (!is.null(default_polygons)) {
      current_polygons(default_polygons)
    }
  })

  # Redraw the map ------------------------------------------------------------
  observe({
    layer <- map_site_layer()
    points <- layer$points
    colors <- layer$colors
    polygons <- normalized_polygons()

    proxy <- leafletProxy("map") |>
      clearGroup("sites") |>
      clearGroup("selected-polygons") |>
      clearGroup("reference-points")

    if (!is.null(points) && nrow(points)) {
      proxy <- proxy |>
        addCircleMarkers(
          data = points,
          radius = 4,
          stroke = TRUE,
          weight = 1,
          color = "#1f2933",
          opacity = 0.8,
          fillOpacity = 0.9,
          fillColor = colors,
          popup = site_popup_label(points),
          group = "sites"
        )
    }

    if (!is.null(polygons) && nrow(polygons)) {
      proxy <- proxy |>
        addPolygons(
          data = polygons,
          group = "selected-polygons",
          fillColor = "#2f5d5a",
          fillOpacity = 0.10,
          color = "#1f2933",
          weight = 2,
          label = paste("Polygon", polygons$polygon_id)
        )
    }

    reference <- reference_points()
    show_reference <- if (is.null(input$show_reference_points)) {
      TRUE
    } else {
      isTRUE(input$show_reference_points)
    }
    if (!is.null(reference) && nrow(reference) && show_reference) {
      proxy <- proxy |>
        addMarkers(
          data = reference,
          icon = reference_cross_icon(),
          popup = reference_popup_label(reference),
          group = "reference-points"
        )
    }

    invisible(proxy)
  })

  # Fit the view only when the data, polygons, or the explicit button change.
  observeEvent(list(selected_sites(), normalized_polygons(), input$zoom_to_data), {
    points <- selected_sites()
    polygons <- normalized_polygons()

    if ((is.null(points) || !nrow(points)) && (is.null(polygons) || !nrow(polygons))) {
      return()
    }

    # When a study area exists, frame that rather than the global site extent.
    proxy <- leafletProxy("map")
    if (!is.null(polygons) && nrow(polygons)) {
      proxy |> fit_map_to_data(polygons = polygons)
    } else {
      proxy |> fit_map_to_data(points = points)
    }
  }, ignoreInit = FALSE)

  observeEvent(reference_points(), {
    reference <- reference_points()
    if (is.null(reference) || !nrow(reference)) return()
    leafletProxy("map") |> fit_map_to_data(points = reference)
  }, ignoreInit = TRUE)

  observeEvent(input$map_draw_new_feature, {
    current_polygons(append_drawn_polygon(
      current_polygons(),
      input$map_draw_new_feature,
      selected_sites()
    ))
  }, ignoreInit = TRUE)

  observeEvent(input$shapefile, {
    upload <- input$shapefile
    req(upload$datapath)

    tryCatch({
      uploaded_polygons <- sf::st_read(upload$datapath, quiet = TRUE)
      normalized <- normalize_polygon_data(uploaded_polygons, selected_sites())

      if (is.null(normalized)) {
        stop("The uploaded file does not contain any polygon features.")
      }

      current_polygons(normalized)
      showNotification(
        sprintf("Loaded %s polygon(s).", nrow(normalized)),
        type = "message",
        duration = 4
      )
    }, error = function(error) {
      showNotification(conditionMessage(error), type = "error", duration = 8)
    })
  }, ignoreInit = TRUE)

  observeEvent(input$clear_polygons_btn, {
    current_polygons(NULL)
    leafletProxy("map") |>
      clearGroup("selected-polygons") |>
      clearGroup("reference-points") |>
      map_draw_toolbar(clear_features = TRUE) |>
      map_draw_toolbar()
    showNotification("Polygons cleared.", type = "message", duration = 3)
  }, ignoreInit = TRUE)

  observeEvent(input$undo_polygon_btn, {
    polygons <- current_polygons()
    if (is.null(polygons) || !nrow(polygons)) {
      showNotification("There are no polygons to remove.", type = "warning", duration = 3)
      return()
    }

    if (nrow(polygons) > 1) {
      current_polygons(polygons[-nrow(polygons), , drop = FALSE])
    } else {
      current_polygons(NULL)
    }
    showNotification("Removed the last polygon.", type = "message", duration = 3)
  }, ignoreInit = TRUE)

  # Analysis ------------------------------------------------------------------
  observeEvent(input$run_btn, {
    tryCatch({
      if (is.null(data_bundle())) stop("BIAD data is still loading. Please try again in a moment.")

      polygons <- current_polygons()
      if (is.null(polygons) || !nrow(polygons)) {
        stop("Draw or upload at least one polygon before running the analysis.")
      }

      result <- withProgress(
        message = "Running ordination analysis",
        value = 0,
        {
          incProgress(0.2, detail = "Preparing grouped records")
          dataset <- prepared_dataset()

          incProgress(0.4, detail = "Aggregating counts and fitting the selected method")
          result <- compute_analysis(
            dataset = dataset,
            polygons = polygons,
            data_type = input$data_type_selector,
            analysis_method = input$analysis_method,
            pca_transform = input$pca_transform,
            pca_scaling = input$pca_scaling
          )

          incProgress(0.3, detail = "Packaging plots and exports")
          result
        }
      )

      analysis_result(result)
      analysis_stale(FALSE)

      overlap_warnings <- grep("^Selected botanical subgroups overlap", result$notes, value = TRUE)
      if (length(overlap_warnings)) {
        showNotification(overlap_warnings[[1]], type = "warning", duration = 10)
      }

      showNotification("Analysis updated.", type = "message", duration = 3)
    }, error = function(error) {
      showNotification(conditionMessage(error), type = "error", duration = 8)
    })
  }, ignoreInit = TRUE)

  # Result outputs ------------------------------------------------------------
  output$plot1 <- renderPlot({
    result <- analysis_result()
    shiny::validate(shiny::need(!is.null(result), "Run analysis to view raw counts."))
    countTotal(result$subregions, result$count_column)
  })

  output$plot2 <- renderPlot({
    result <- analysis_result()
    shiny::validate(shiny::need(!is.null(result), "Run analysis to view taxon composition."))
    bySpeciesComposition(result$subregions, result$count_column)
  })

  output$plot2_legend <- renderPlot({
    result <- analysis_result()
    shiny::validate(shiny::need(!is.null(result), "Run analysis to view taxon composition."))
    plotTaxonGroupLegend(result$subregions)
  })

  output$plot2_culture <- renderPlot({
    result <- analysis_result()
    shiny::validate(shiny::need(!is.null(result), "Run analysis to view culture composition."))
    byCultureComposition(
      result$subregions,
      result$count_column,
      culture_limit = input$plot_culture_composition_limit
    )
  })

  output$plot_culture_composition_limit_ui <- renderUI({
    result <- analysis_result()
    if (is.null(result) || is.null(result$culture_projection)) return(NULL)

    culture_summary <- result$culture_projection$summary
    if (is.null(culture_summary) || nrow(culture_summary) <= 1) return(NULL)

    max_cultures <- nrow(culture_summary)
    selected_limit <- isolate(input$plot_culture_composition_limit)
    if (is.null(selected_limit) || !is.numeric(selected_limit) || selected_limit < 1) {
      selected_limit <- min(20, max_cultures)
    }
    selected_limit <- min(selected_limit, max_cultures)

    sliderInput(
      "plot_culture_composition_limit",
      "Cultures shown (top by total counts)",
      min = 1,
      max = max_cultures,
      value = selected_limit,
      step = 1,
      width = "100%"
    )
  })

  output$plot_culture_limit_ui <- renderUI({
    result <- analysis_result()
    if (!isTRUE(input$plot_label_culture) || is.null(result) || is.null(result$culture_projection)) {
      return(NULL)
    }

    culture_summary <- result$culture_projection$summary
    if (is.null(culture_summary) || nrow(culture_summary) <= 1) return(NULL)

    max_cultures <- nrow(culture_summary)
    selected_limit <- isolate(input$plot_culture_limit)
    if (is.null(selected_limit) || !is.numeric(selected_limit) || selected_limit < 1) {
      selected_limit <- min(20, max_cultures)
    }
    selected_limit <- min(selected_limit, max_cultures)

    sliderInput(
      "plot_culture_limit",
      "Top cultures by number of grouped records",
      min = 1,
      max = max_cultures,
      value = selected_limit,
      step = 1,
      width = "100%"
    )
  })

  output$plot3 <- renderPlot({
    result <- analysis_result()
    shiny::validate(shiny::need(!is.null(result), "Run analysis to view the ordination map."))

    culture_coordinates <- NULL
    if (isTRUE(input$plot_label_culture) && !is.null(result$culture_projection)) {
      culture_summary <- result$culture_projection$summary
      culture_coordinates <- result$culture_projection$coordinates

      if (!is.null(culture_summary) && nrow(culture_summary)) {
        culture_limit <- input$plot_culture_limit
        if (is.null(culture_limit) || !is.numeric(culture_limit) || culture_limit < 1) {
          culture_limit <- min(20, nrow(culture_summary))
        }
        selected_cultures <- head(rownames(culture_summary), min(culture_limit, nrow(culture_summary)))
        culture_coordinates <- culture_coordinates[selected_cultures, , drop = FALSE]
      }
    }

    plotOrdinationArrows(
      result$ordination,
      period_levels = result$period_levels,
      culture_coordinates = culture_coordinates,
      show_year_labels = isTRUE(input$plot_label_year),
      show_taxa_labels = isTRUE(input$plot_label_taxa),
      show_culture_labels = isTRUE(input$plot_label_culture)
    )
  })

  output$plot4 <- renderPlot({
    result <- analysis_result()
    shiny::validate(shiny::need(!is.null(result), "Run analysis to view ordination dimensions."))
    plotOrdinationDimensions(result$ordination, period_levels = result$period_levels)
  })

  output$plot4_diag <- renderPlot({
    result <- analysis_result()
    shiny::validate(shiny::need(!is.null(result), "Run analysis to view ordination diagnostics."))
    plotOrdinationDiagnostics(result$ordination)
  })

  output$analysis_preview <- DT::renderDT({
    result <- analysis_result()
    shiny::validate(shiny::need(!is.null(result), "Run analysis to preview the export table."))

    DT::datatable(
      utils::head(result$combined_data, 30),
      rownames = FALSE,
      class = "compact stripe nowrap",
      options = list(scrollX = TRUE, pageLength = 10, dom = "tip")
    )
  })

  # Downloads -----------------------------------------------------------------
  output$download_polygons <- downloadHandler(
    filename = function() "selected_polygons.gpkg",
    content = function(file) {
      polygons <- current_polygons()
      if (is.null(polygons) || !nrow(polygons)) {
        stop("No polygons are available to download.")
      }
      sf::st_write(polygons, file, driver = "GPKG", quiet = TRUE)
    }
  )

  output$download_table <- downloadHandler(
    filename = function() sprintf("aggregated_%s_taxa.csv", tolower(input$data_type_selector)),
    content = function(file) {
      result <- analysis_result()
      if (is.null(result)) stop("Run the analysis before downloading results.")
      utils::write.csv(result$combined_data, file, row.names = FALSE, na = "")
    }
  )

  output$download_raw_table <- downloadHandler(
    filename = function() sprintf("raw_phase_assignments_%s_taxa.csv", tolower(input$data_type_selector)),
    content = function(file) {
      result <- analysis_result()
      if (is.null(result)) stop("Run the analysis before downloading phase assignments.")
      utils::write.csv(result$phase_assignments, file, row.names = FALSE, na = "")
    }
  )

  # Grouping management -------------------------------------------------------
  observeEvent(input$upload_group, {
    upload <- input$upload_group
    req(upload)

    tryCatch({
      file_name <- basename(upload$name)
      if (!grepl("\\.csv$", file_name, ignore.case = TRUE)) {
        stop("Please upload a .csv grouping file.")
      }

      parsed <- parse_grouping_file(upload$datapath)
      if (!length(parsed)) {
        stop(paste0(
          "No groups found. The file needs at least one row with a group number, ",
          "a label, and taxa codes separated by '+'."
        ))
      }

      type_dir <- grouping_directory(input$data_type_selector)
      dir.create(type_dir, recursive = TRUE, showWarnings = FALSE)
      destination <- file.path(type_dir, file_name)

      if (!file.copy(upload$datapath, destination, overwrite = TRUE)) {
        stop("The uploaded file could not be saved.")
      }

      invalidate_grouping_cache()
      grouping_revision(grouping_revision() + 1)
      group_catalog_version(group_catalog_version() + 1)
      refresh_grouping_inputs(destination)

      showNotification(
        sprintf("Uploaded '%s' — %d groups found.", file_name, length(parsed)),
        type = "message",
        duration = 5
      )
    }, error = function(error) {
      showNotification(conditionMessage(error), type = "error", duration = 8)
    })
  }, ignoreInit = TRUE)

  output$download_group <- downloadHandler(
    filename = function() {
      path <- active_grouping_path()
      if (is.null(path)) "group.csv" else basename(path)
    },
    content = function(file) {
      path <- active_grouping_path()
      if (is.null(path) || !file.exists(path)) {
        stop("No grouping file is currently selected or available.")
      }
      if (!file.copy(path, file, overwrite = TRUE)) {
        stop("The grouping file could not be copied.")
      }
    }
  )

  output$download_template <- downloadHandler(
    filename = function() "template_grouping.csv",
    content = function(file) {
      template <- "groupings/template_grouping.csv"
      if (!file.exists(template)) stop("The grouping template is missing.")
      if (!file.copy(template, file, overwrite = TRUE)) {
        stop("The grouping template could not be copied.")
      }
    }
  )

  observeEvent(input$edit_group, {
    path <- active_grouping_path()
    if (is.null(path) || !file.exists(path)) {
      showNotification("Grouping file not found.", type = "error", duration = 4)
      return()
    }

    parsed <- tryCatch(read_grouping_rows(path), error = function(error) {
      showNotification(conditionMessage(error), type = "error", duration = 6)
      NULL
    })
    if (is.null(parsed)) return()

    group_edit_store$path <- path
    group_edit_store$preamble <- parsed$raw[!parsed$is_group, , drop = FALSE]
    group_edit_store$groups <- parsed$raw[parsed$is_group, , drop = FALSE]

    showModal(modalDialog(
      title = paste("Edit grouping:", basename(path)),
      size = "l",
      easyClose = TRUE,
      div(
        class = "modal-help",
        "Edit group labels and taxon codes. Separate taxon codes with '+'. ",
        "Title and note rows at the top of the file are preserved."
      ),
      DT::DTOutput("edit_group_table"),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("save_group_edits", "Save changes", class = "btn-primary")
      )
    ))
  }, ignoreInit = TRUE)

  output$edit_group_table <- DT::renderDT({
    groups <- group_edit_store$groups
    if (is.null(groups) || !nrow(groups)) {
      return(DT::datatable(data.frame(Message = "No editable group rows found.")))
    }

    extra_names <- if (ncol(groups) > 3) paste0("Extra ", seq_len(ncol(groups) - 3)) else character(0)

    DT::datatable(
      groups,
      editable = list(target = "cell"),
      rownames = FALSE,
      colnames = c("Group", "Label", "Taxa (+ separated)", extra_names),
      options = list(scrollX = TRUE, pageLength = 15, dom = "ftip")
    )
  })

  # Accumulate edits without re-rendering the table (which would reset the
  # user's scroll position and focus).
  observeEvent(input$edit_group_table_cell_edit, {
    groups <- group_edit_store$groups
    req(groups)
    info <- input$edit_group_table_cell_edit

    groups[info$row, info$col + 1] <- DT::coerceValue(
      info$value,
      groups[info$row, info$col + 1]
    )
    group_edit_store$groups <- groups
  }, ignoreInit = TRUE)

  observeEvent(input$save_group_edits, {
    path <- group_edit_store$path
    if (is.null(path) || !file.exists(path)) {
      showNotification("Cannot save: grouping file not found.", type = "error", duration = 5)
      return()
    }

    tryCatch({
      group_numbers <- trimws(as.character(group_edit_store$groups[[1]]))
      if (any(!grepl("^[0-9]+[.)]?$", group_numbers))) {
        stop("Every group row needs a numeric group number in the first column.")
      }

      merged <- merge_grouping_rows(group_edit_store$preamble, group_edit_store$groups)
      write_grouping_rows(merged, path)

      invalidate_grouping_cache()
      grouping_revision(grouping_revision() + 1)
      group_catalog_version(group_catalog_version() + 1)
      refresh_grouping_inputs(path)

      removeModal()
      showNotification("Grouping file saved.", type = "message", duration = 4)
    }, error = function(error) {
      showNotification(conditionMessage(error), type = "error", duration = 8)
    })
  }, ignoreInit = TRUE)

  # Reset ---------------------------------------------------------------------
  observeEvent(input$reset_app, {
    current_polygons(NULL)
    reference_csv(NULL)
    reference_points(NULL)
    reference_info(NULL)
    reference_signature(NULL)
    analysis_result(NULL)
    analysis_stale(FALSE)
    default_polygons_attempted(FALSE)

    updateSelectInput(session, "data_type_selector", selected = "Faunal")
    faunal_choices <- available_groupings("Faunal")
    faunal_default <- preferred_grouping_path("Faunal")
    updateSelectInput(
      session, "file_selector_tx",
      choices = faunal_choices,
      selected = faunal_default
    )
    updateSelectInput(
      session, "file_selector_tx_combined_faunal",
      choices = faunal_choices,
      selected = faunal_default
    )
    updateSelectInput(session, "file_selector_per", selected = "auto")
    updateSelectInput(session, "time_preset", selected = "custom")
    updateNumericInput(session, "start_value", value = 8500)
    updateNumericInput(session, "end_value", value = 3000)
    updateSliderInput(session, "duration", value = 500)
    updateSelectInput(session, "analysis_method", selected = "ca")
    updateSelectInput(session, "pca_transform", selected = "log1p")
    updateSelectInput(session, "pca_scaling", selected = "scale")
    updateCheckboxInput(session, "plot_label_year", value = FALSE)
    updateCheckboxInput(session, "plot_label_taxa", value = FALSE)
    updateCheckboxInput(session, "plot_label_culture", value = FALSE)
    updateCheckboxInput(session, "only_selected_sites", value = FALSE)
    updateSelectizeInput(session, "combined_botanical_groups", selected = character(0))

    leafletProxy("map") |>
      clearGroup("selected-polygons") |>
      clearGroup("reference-points") |>
      map_draw_toolbar(clear_features = TRUE) |>
      map_draw_toolbar()

    showNotification("Workspace reset.", type = "message", duration = 3)
  }, ignoreInit = TRUE)

  observeEvent(input$time_preset, {
    if (identical(input$time_preset, "custom")) return()

    parts <- strsplit(input$time_preset, "_", fixed = TRUE)[[1]]
    if (length(parts) != 3) return()

    updateNumericInput(session, "start_value", value = as.numeric(parts[[1]]))
    updateNumericInput(session, "end_value", value = as.numeric(parts[[2]]))
    updateSliderInput(session, "duration", value = as.numeric(parts[[3]]))
  }, ignoreInit = TRUE)
}
