section_card <- function(title, subtitle = NULL, ..., step = NULL, class_name = NULL) {
  div(
    class = paste("card", class_name),
    if (!is.null(step)) div(class = "step-label", step),
    tags$h3(class = "card-title", title),
    if (!is.null(subtitle)) tags$p(class = "helper-text", subtitle),
    ...
  )
}

metric_chip <- function(label, value) {
  div(
    class = "metric-chip",
    tags$span(class = "metric-label", label),
    tags$span(class = "metric-value", as.character(value))
  )
}

metrics_grid <- function(...) {
  div(class = "metrics-grid", ...)
}

status_pill <- function(text, tone = c("neutral", "fresh", "stale", "warning")) {
  tone <- match.arg(tone)
  div(class = paste("status-pill", paste0("status-pill--", tone)), text)
}

build_group_summary_ui <- function(summary_stats) {
  groups <- summary_stats$groups

  if (nrow(groups) == 0) {
    return(
      div(
        class = "group-summary",
        tags$p(
          class = "helper-text helper-text--compact",
          "This grouping file does not map any taxa in the current dataset."
        )
      )
    )
  }

  group_items <- lapply(seq_len(nrow(groups)), function(index) {
    taxa <- groups$taxa[[index]]

    tags$details(
      class = "group-summary-item",
      tags$summary(
        class = "group-summary-toggle",
        tags$span(class = "group-summary-name", groups$group_name[[index]]),
        tags$span(class = "group-summary-count", paste(groups$taxon_count[[index]], "taxa"))
      ),
      div(
        class = "group-summary-taxa",
        if (length(taxa) == 0) {
          tags$span(class = "group-summary-empty", "No taxa listed.")
        } else {
          do.call(
            div,
            c(
              list(class = "group-summary-taxa-list"),
              lapply(taxa, function(taxon_label) {
                tags$span(class = "group-summary-taxon", taxon_label)
              })
            )
          )
        }
      )
    )
  })

  note <- if (summary_stats$unmatched > 0) {
    tags$p(
      class = "helper-text helper-text--compact",
      sprintf("%s taxa remain unmatched by this grouping file.", summary_stats$unmatched)
    )
  } else {
    tags$p(
      class = "helper-text helper-text--compact",
      "All taxa in the current dataset are covered by this grouping file."
    )
  }

  div(
    class = "group-summary",
    div(
      class = "group-summary-header",
      tags$span(class = "group-summary-title", "Grouping composition"),
      tags$span(class = "group-summary-meta", paste(nrow(groups), "groups"))
    ),
    tags$p(
      class = "group-summary-help",
      "Click a taxa count to expand the exact taxa list."
    ),
    do.call(div, c(list(class = "group-summary-list"), group_items)),
    note
  )
}

build_selection_summary_ui <- function(
  data_type,
  site_count,
  record_count,
  polygon_count,
  grouping_label,
  stale,
  has_results
) {
  status <- if (polygon_count == 0) {
    status_pill("Add at least one polygon to define a study area.", "warning")
  } else if (has_results && stale) {
    status_pill("Inputs changed — click Run analysis to refresh the charts.", "stale")
  } else if (has_results) {
    status_pill("Results are synced with the current map and controls.", "fresh")
  } else {
    status_pill("Ready to run once the study area looks right.", "neutral")
  }

  metrics <- metrics_grid(
    metric_chip("Dataset", data_type),
    metric_chip("Sites", format(site_count, big.mark = ",")),
    metric_chip("Records", format(record_count, big.mark = ",")),
    metric_chip("Polygons", polygon_count),
    metric_chip("Grouping", grouping_label)
  )

  div(class = "summary-block", status, metrics)
}

build_analysis_status_ui <- function(result, stale) {
  if (is.null(result)) {
    return(
      div(
        class = "summary-block",
        status_pill("No analysis has been run yet.", "neutral"),
        tags$p(
          class = "helper-text",
          "Choose a dataset, define polygons, and click Run analysis to compute the correspondence analysis."
        )
      )
    )
  }

  status <- if (stale) {
    status_pill("Showing the last successful run. Rerun to sync with current inputs.", "stale")
  } else {
    status_pill("Results are up to date.", "fresh")
  }

  metrics <- metrics_grid(
    metric_chip("Used polygons", result$summary$polygon_count),
    metric_chip("Phases", result$summary$phase_count),
    metric_chip("Taxon groups", result$summary$taxon_group_count),
    metric_chip("Records analysed", format(result$summary$records, big.mark = ","))
  )

  notes <- NULL
  if (length(result$notes) > 0) {
    notes <- do.call(
      tags$ul,
      c(
        list(class = "note-list"),
        lapply(result$notes, tags$li)
      )
    )
  }

  div(class = "summary-block", status, metrics, notes)
}

app_ui <- function() {
  bslib::page_fluid(
    title = "BIAD Taxa Explorer",
    theme = bslib::bs_theme(
      version = 5,
      base_font = bslib::font_google("Inter"),
      bg = "#eef3f8",
      fg = "#142133",
      primary = "#1976D2",
      secondary = "#111827",
      success = "#2e7d32",
      info = "#1976D2",
      warning = "#f59e0b",
      border_radius = "0.35rem",
      btn_border_radius = "0.25rem"
    ),
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "app.css")
    ),
    div(
      class = "page-wrap",
      div(
        class = "hero",
        div(
          class = "hero-brand",
          tags$a(
            href = "https://biadwiki.org",
            target = "_blank",
            rel = "noopener noreferrer",
            class = "hero-logo-link",
            tags$img(
              src = "biad.logo.png",
              alt = "BIAD logo",
              class = "hero-logo"
            )
          ),
          div(
            class = "hero-copy",
            tags$h1("BIAD Taxa Explorer"),
            tags$p(
              "Select a dataset, define study areas, group taxa, and compare spatial-temporal patterns with correspondence analysis."
            ),
            div(
              class = "warning-banner",
              "Experimental internal interface for research and testing."
            )
          )
        )
      ),
      fluidRow(
        column(
          width = 4,
          div(
            class = "app-sidebar-shell",
            div(
              class = "sticky-column",
              section_card(
                title = "Data and study areas",
                subtitle = "Choose the dataset, taxon grouping, and polygons to analyse.",
                step = "Step 1",
                class_name = "sidebar-card",
                div(
                  class = "sidebar-subsection",
                  selectInput(
                    "data_type_selector",
                    "Dataset",
                    choices = c("Faunal", "Botanical"),
                    selected = "Faunal"
                  ),
                  selectInput("file_selector_tx", "Taxon grouping", choices = character(0)),
                  tags$details(
                    style = "margin-top: 4px; font-size: 11px;",
                    tags$summary(
                      style = "color: #999; cursor: pointer; user-select: none;",
                      "group file options"
                    ),
                    div(
                      style = "margin-top: 6px; padding: 6px; background: #f9f9f9; border-radius: 4px; border: 1px solid #eee;",
                      fileInput("upload_group", NULL,
                        accept = ".csv",
                        buttonLabel = "Upload .csv",
                        placeholder = "no file selected",
                        width = "100%"
                      ),
                      div(
                        style = "display: flex; gap: 8px; margin-top: -10px;",
                        downloadButton("download_group", "↓ download",
                          style = "font-size: 11px; padding: 2px 8px; height: auto; color: #555; background: #fff; border: 1px solid #ccc;"
                        ),
                        actionButton("edit_group", "✎ edit labels",
                          style = "font-size: 11px; padding: 2px 8px; height: auto; color: #555; background: #fff; border: 1px solid #ccc;"
                        )
                      )
                    )
                  ),
                  uiOutput("taxon_table")
                ),
                div(
                  class = "sidebar-subsection sidebar-action-surface",
                  tags$div(class = "sidebar-subsection-title", "Study area polygons"),
                  tags$p(
                    class = "sidebar-subsection-help",
                    "Upload a GeoPackage here, or draw polygons directly on the map."
                  ),
                  fileInput("shapefile", "Upload polygons (.gpkg)", accept = ".gpkg"),
                  actionButton(
                    "clear_polygons_btn",
                    "Clear polygons",
                    class = "btn-default btn-block sidebar-secondary-btn sidebar-danger-btn"
                  )
                )
              ),
              section_card(
                title = "Time controls",
                subtitle = "Pick predefined archaeological periods or build fixed BP slices.",
                step = "Step 2",
                class_name = "sidebar-card",
                div(
                  class = "sidebar-subsection",
                  selectInput(
                    "file_selector_per",
                    "Period grouping",
                    choices = c(
                      "Automatic time slices" = "auto",
                      "Named archaeological periods" = "groupings/periods/periods.csv"
                    ),
                    selected = "auto"
                  )
                ),
                conditionalPanel(
                  condition = "input.file_selector_per == 'auto'",
                  div(
                    class = "sidebar-subsection sidebar-time-surface",
                    div(
                      class = "time-controls-panel",
                      div(
                        class = "time-boundary-grid",
                        div(
                          class = "time-boundary-field",
                          numericInput("start_value", "Start (BP)", value = 8500, min = 0, width = "100%")
                        ),
                        div(
                          class = "time-boundary-field",
                          numericInput("end_value", "End (BP)", value = 3000, min = 0, width = "100%")
                        )
                      ),
                      div(
                        class = "slider-panel slider-panel--highlight",
                        div(
                          class = "slider-panel-heading",
                          tags$div(class = "slider-panel-label", "Time-slice width"),
                          tags$p(
                            class = "slider-panel-help",
                            "Move the slider to choose the width of each BP slice."
                          )
                        ),
                        sliderInput(
                          "duration",
                          NULL,
                          min = 100,
                          max = 1000,
                          value = 500,
                          step = 50,
                          width = "100%"
                        )
                      )
                    )
                  )
                ),
                div(
                  class = "sidebar-option-row",
                  checkboxInput("use_logs", "Apply log(count + 1) before CA", FALSE)
                )
              ),
              section_card(
                title = "Run and export",
                subtitle = "Analysis only reruns when you ask for it, so the interface stays responsive.",
                step = "Step 3",
                class_name = "sidebar-card sidebar-actions-card",
                div(
                  class = "control-actions",
                  div(
                    class = "sidebar-subsection sidebar-run-surface",
                    tags$div(class = "sidebar-subsection-title", "Run analysis"),
                    tags$p(
                      class = "sidebar-subsection-help",
                      "Charts and tables refresh only when you click Run analysis."
                    ),
                    actionButton(
                      "run_btn",
                      "Run analysis",
                      class = "btn-primary btn-lg btn-block sidebar-primary-btn"
                    )
                  ),
                  div(
                    class = "sidebar-subsection sidebar-export-surface",
                    tags$div(class = "sidebar-subsection-title", "Export"),
                    tags$p(
                      class = "sidebar-subsection-help",
                      "Download the current study polygons or the aggregated results table."
                    ),
                    div(
                      class = "export-actions",
                      downloadButton(
                        "download_polygons",
                        "Download polygons",
                        class = "btn-default btn-block sidebar-download-btn"
                      ),
                      downloadButton(
                        "download_table",
                        "Download results",
                        class = "btn-default btn-block sidebar-download-btn"
                      )
                    )
                  )
                )
              )
            )
          )
        ),
        column(
          width = 8,
          div(
            class = "app-main-shell",
            section_card(
              title = "Map and coverage",
              subtitle = "Draw one or more polygons directly on the map or upload a GeoPackage.",
              uiOutput("selection_summary"),
              leafletOutput("map", height = "62vh")
            ),
            section_card(
              title = "Analysis results",
              subtitle = "Compare raw counts, taxon composition, and correspondence-analysis structure.",
              uiOutput("analysis_status"),
              div(
                class = "results-tabs",
                tabsetPanel(
                  id = "results_tabs",
                  tabPanel(
                    "Raw counts",
                    tags$p(
                      class = "plot-caption",
                      "Total counts per period for each selected polygon."
                    ),
                    plotOutput("plot1", height = "320px")
                  ),
                  tabPanel(
                    "Composition",
                    tags$p(
                      class = "plot-caption",
                      "Relative share of taxon groups across periods within each polygon."
                    ),
                    plotOutput("plot2", height = "340px")
                  ),
                  tabPanel(
                    "CA map",
                    tags$p(
                      class = "plot-caption",
                      "Correspondence analysis with arrows that track phase order inside each polygon."
                    ),
                    plotOutput("plot3", height = "520px")
                  ),
                  tabPanel(
                    "Dimensions",
                    tags$p(
                      class = "plot-caption",
                      "How the first two CA dimensions evolve across time for each polygon."
                    ),
                    plotOutput("plot4", height = "520px")
                  ),
                  tabPanel(
                    "Data preview",
                    tags$p(
                      class = "plot-caption",
                      "Preview of the exported aggregated table."
                    ),
                    tableOutput("analysis_preview")
                  )
                )
              )
            )
          )
        )
      )
    )
  )
}
