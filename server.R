server <- function(input, output, session) {
  if (exists("check.conn", mode = "function")) {
    check.conn(conn)
  }

  current_polygons <- reactiveVal(NULL)
  analysis_result <- reactiveVal(NULL)
  analysis_stale <- reactiveVal(FALSE)

  selected_dataset <- reactive({
    if (identical(input$data_type_selector, "Botanical")) {
      allBTaxa
    } else {
      allFTaxa
    }
  })

  selected_sites <- reactive({
    if (identical(input$data_type_selector, "Botanical")) {
      botanicalSites
    } else {
      faunalSites
    }
  })

  observe({
    grouping_choices <- available_groupings(input$data_type_selector)
    if (!length(grouping_choices)) {
      return()
    }

    selected_grouping <- isolate(input$file_selector_tx)
    if (is.null(selected_grouping) || !selected_grouping %in% unname(grouping_choices)) {
      selected_grouping <- unname(grouping_choices)[1]
    }

    updateSelectInput(
      session,
      "file_selector_tx",
      choices = grouping_choices,
      selected = selected_grouping
    )
  })

  grouped_dataset <- reactive({
    req(input$file_selector_tx)
    groupTaxons(selected_dataset(), input$file_selector_tx, logging = FALSE)
  })

  prepared_dataset <- reactive({
    dataset <- grouped_dataset()
    dataset$new_periods <- split_periods(
      dataset = dataset,
      selection = input$file_selector_per,
      start = input$start_value,
      end = input$end_value,
      len = input$duration
    )
    dataset
  })

  observe({
    input$data_type_selector
    input$file_selector_tx
    input$file_selector_per
    input$start_value
    input$end_value
    input$duration
    input$use_logs
    current_polygons()

    if (!is.null(isolate(analysis_result()))) {
      analysis_stale(TRUE)
    }
  })

  output$taxon_table <- renderUI({
    build_group_summary_ui(group_summary_stats(grouped_dataset()))
  })

  output$selection_summary <- renderUI({
    polygon_count <- if (is.null(current_polygons())) 0 else nrow(current_polygons())
    grouping_label <- if (is.null(input$file_selector_tx)) {
      "—"
    } else {
      format_grouping_label(input$file_selector_tx)
    }

    build_selection_summary_ui(
      data_type = input$data_type_selector,
      site_count = nrow(selected_sites()),
      record_count = nrow(selected_dataset()),
      polygon_count = polygon_count,
      grouping_label = grouping_label,
      stale = analysis_stale(),
      has_results = !is.null(analysis_result())
    )
  })

  output$analysis_status <- renderUI({
    build_analysis_status_ui(analysis_result(), analysis_stale())
  })

  output$map <- renderLeaflet({
    leaflet() |>
      addTiles() |>
      fit_map_to_data(points = selected_sites()) |>
      addDrawToolbar(
        targetGroup = "selected-polygons",
        polygonOptions = drawPolygonOptions(showArea = TRUE),
        editOptions = editToolbarOptions(edit = FALSE, remove = FALSE),
        polylineOptions = FALSE,
        rectangleOptions = FALSE,
        circleOptions = FALSE,
        markerOptions = FALSE,
        circleMarkerOptions = FALSE
      )
  })

  observeEvent(selected_sites(), {
    leafletProxy("map") |>
      fit_map_to_data(points = selected_sites())
  }, ignoreInit = TRUE)

  observe({
    points <- selected_sites()
    polygons <- normalize_polygon_data(current_polygons(), points)
    point_colors <- site_colors(points, polygons, areapal, normalize = FALSE)

    proxy <- leafletProxy("map") |>
      clearGroup("sites") |>
      clearGroup("selected-polygons")

    if (nrow(points) > 0) {
      proxy <- proxy |>
        addCircleMarkers(
          data = points,
          radius = 4,
          stroke = FALSE,
          fillOpacity = 0.9,
          color = point_colors,
          popup = site_popup_label(points),
          group = "sites"
        )
    }

    if (!is.null(polygons) && nrow(polygons) > 0) {
      proxy |>
        addPolygons(
          data = polygons,
          group = "selected-polygons",
          fillColor = "#1d4ed8",
          fillOpacity = 0.08,
          color = "#1f2937",
          weight = 2
        )
    }
  })

  observeEvent(input$map_draw_new_feature, {
    updated_polygons <- append_drawn_polygon(
      current_polygons(),
      input$map_draw_new_feature,
      selected_sites()
    )
    current_polygons(updated_polygons)
  }, ignoreInit = TRUE)

  observeEvent(input$shapefile, {
    req(input$shapefile$datapath)

    tryCatch({
      uploaded_polygons <- sf::st_read(input$shapefile$datapath, quiet = TRUE)
      normalized_polygons <- normalize_polygon_data(uploaded_polygons, selected_sites())

      if (is.null(normalized_polygons)) {
        stop("The uploaded file does not contain any polygon features.")
      }

      current_polygons(normalized_polygons)
      showNotification(
        sprintf("Loaded %s polygon(s).", nrow(normalized_polygons)),
        type = "message",
        duration = 3
      )
    }, error = function(error) {
      showNotification(conditionMessage(error), type = "error", duration = 6)
    })
  }, ignoreInit = TRUE)

  observeEvent(input$clear_polygons_btn, {
    current_polygons(NULL)
    showNotification("Polygons cleared.", type = "message", duration = 2)
  }, ignoreInit = TRUE)

  observeEvent(input$run_btn, {
    tryCatch({
      result <- withProgress(message = "Running correspondence analysis", value = 0, {
        incProgress(0.25, detail = "Preparing grouped records")
        dataset <- prepared_dataset()
        polygons <- current_polygons()

        incProgress(0.45, detail = "Aggregating counts and fitting CA")
        result <- compute_analysis(
          dataset = dataset,
          polygons = polygons,
          data_type = input$data_type_selector,
          use_logs = input$use_logs
        )

        incProgress(0.30, detail = "Packaging plots and exports")

        result
      })

      analysis_result(result)
      analysis_stale(FALSE)
      showNotification("Analysis updated.", type = "message", duration = 3)
    }, error = function(error) {
      showNotification(conditionMessage(error), type = "error", duration = 8)
    })
  }, ignoreInit = TRUE)

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

  output$plot3 <- renderPlot({
    result <- analysis_result()
    shiny::validate(shiny::need(!is.null(result), "Run analysis to view the CA map."))
    plotCAarrows(result$ca, period_levels = result$period_levels)
  })

  output$plot4 <- renderPlot({
    result <- analysis_result()
    shiny::validate(shiny::need(!is.null(result), "Run analysis to view CA dimensions."))
    plot2dim(result$ca, period_levels = result$period_levels)
  })

  output$analysis_preview <- renderTable({
    result <- analysis_result()
    shiny::validate(shiny::need(!is.null(result), "Run analysis to preview the export table."))
    utils::head(result$combined_data, 12)
  }, rownames = FALSE, striped = TRUE, bordered = TRUE)

  output$download_polygons <- downloadHandler(
    filename = function() {
      "selected_polygons.gpkg"
    },
    content = function(file) {
      polygons <- current_polygons()
      if (is.null(polygons) || nrow(polygons) == 0) {
        stop("No polygons are available to download.")
      }

      sf::st_write(polygons, file, driver = "GPKG", quiet = TRUE)
    }
  )

  output$download_table <- downloadHandler(
    filename = function() {
      sprintf("aggregated_%s_taxa.csv", tolower(input$data_type_selector))
    },
    content = function(file) {
      result <- analysis_result()
      if (is.null(result)) {
        stop("Run the analysis before downloading results.")
      }

      utils::write.csv(result$combined_data, file, row.names = FALSE)
    }
  )

  # ── Group management: upload ──────────────────────────────────────────────────
  observeEvent(input$upload_group, {
    req(input$upload_group$datapath)

    tryCatch({
      # Validate: must parse as a usable grouping file
      grouping <- read_taxon_grouping(input$upload_group$datapath)
      if (length(grouping) == 0) {
        stop("No groups found. The file must have at least one row with a group number, label, and taxa codes separated by '+'.")
      }

      type_dir <- if (identical(input$data_type_selector, "Botanical")) {
        "groupings/botanical_taxa"
      } else {
        "groupings/faunal_taxa"
      }

      if (!dir.exists(type_dir)) {
        dir.create(type_dir, recursive = TRUE)
      }

      # Ensure filename starts with "group_" so it appears in the selector
      dest_name <- input$upload_group$name
      if (!grepl("^group_", dest_name)) {
        dest_name <- paste0("group_", dest_name)
      }

      dest_path <- file.path(type_dir, dest_name)
      file.copy(input$upload_group$datapath, dest_path, overwrite = TRUE)

      # Refresh the grouping selector
      grouping_choices <- available_groupings(input$data_type_selector)
      updateSelectInput(session, "file_selector_tx",
                        choices = grouping_choices,
                        selected = dest_path)

      showNotification(
        sprintf("Uploaded '%s' — %d groups found.", dest_name, length(grouping)),
        type = "message", duration = 4
      )
    }, error = function(error) {
      showNotification(conditionMessage(error), type = "error", duration = 6)
    })
  }, ignoreInit = TRUE)

  # ── Group management: download ────────────────────────────────────────────────
  output$download_group <- downloadHandler(
    filename = function() {
      if (is.null(input$file_selector_tx)) "group.csv" else basename(input$file_selector_tx)
    },
    content = function(file) {
      path <- input$file_selector_tx
      if (is.null(path) || !file.exists(path)) {
        stop("No grouping file is currently selected or available.")
      }
      file.copy(path, file)
    }
  )

  # ── Group management: download template ──────────────────────────────────────
  output$download_template <- downloadHandler(
    filename = function() "template_grouping.csv",
    content = function(file) {
      file.copy("groupings/template_grouping.csv", file)
    }
  )

  # ── Group management: edit modal ──────────────────────────────────────────────
  observeEvent(input$edit_group, {
    req(input$file_selector_tx)

    path <- input$file_selector_tx
    if (!file.exists(path)) {
      showNotification("Grouping file not found.", type = "error", duration = 4)
      return()
    }

    group_data <- utils::read.csv(path, header = FALSE, fill = TRUE)

    showModal(modalDialog(
      title = paste("Edit:", basename(path)),
      size = "l",
      DT::DTOutput("edit_group_table"),
      footer = tagList(
        actionButton("save_group_edits", "Save", class = "btn-primary"),
        modalButton("Cancel")
      )
    ))

    output$edit_group_table <- DT::renderDT({
      DT::datatable(
        group_data,
        editable = list(target = "cell"),
        rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 20)
      )
    })
  }, ignoreInit = TRUE)

  observeEvent(input$save_group_edits, {
    req(input$edit_group_table_cell_edit)
    req(input$file_selector_tx)

    path <- input$file_selector_tx
    if (!file.exists(path)) {
      showNotification("Cannot save: grouping file not found.", type = "error", duration = 4)
      return()
    }

    tryCatch({
      group_data <- utils::read.csv(path, header = FALSE, fill = TRUE)
      info <- input$edit_group_table_cell_edit
      group_data[info$row, info$col + 1] <- DT::coerceValue(info$value, group_data[info$row, info$col + 1])
      utils::write.csv(group_data, path, row.names = FALSE)

      removeModal()
      showNotification("Group file saved.", type = "message", duration = 3)
    }, error = function(error) {
      showNotification(conditionMessage(error), type = "error", duration = 6)
    })
  }, ignoreInit = TRUE)
}
