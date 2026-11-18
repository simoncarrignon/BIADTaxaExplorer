# UI construction --------------------------------------------------------------

app_theme <- function() {
  bslib::bs_theme(
    version = 5,
    base_font = "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif",
    heading_font = "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif",
    bg = "#f4f4f2",
    fg = "#1f2933",
    primary = "#2f5d5a",
    secondary = "#4b5563",
    success = "#2f6b3f",
    info = "#2f5d5a",
    warning = "#8a5a12",
    danger = "#9b2c2c",
    border_radius = "0.25rem",
    btn_border_radius = "0.25rem"
  )
}

# Small building blocks --------------------------------------------------------

app_icon <- function(name) shiny::icon(name)

status_pill <- function(text, tone = c("neutral", "fresh", "stale", "warning", "loading")) {
  tone <- match.arg(tone)
  div(
    class = paste("status-line", paste0("status-line--", tone)),
    tags$span(class = "status-dot"),
    tags$span(text)
  )
}

metric_chip <- function(label, value, icon_name = NULL) {
  div(
    class = "metric",
    tags$span(class = "metric-label", label),
    tags$span(class = "metric-value", as.character(value))
  )
}

metrics_grid <- function(...) {
  div(class = "stats-row", ...)
}

helper_text <- function(...) {
  tags$p(class = "helper-text", ...)
}

empty_state <- function(title, body, icon_name = "circle-info") {
  div(
    class = "empty-state",
    tags$h4(title),
    tags$p(body)
  )
}

panel_title <- function(text, icon_name = NULL) {
  tags$div(class = "panel-title",
    if (!is.null(icon_name)) app_icon(icon_name),
    tags$span(text)
  )
}

subsection_label <- function(text, help = NULL) {
  tagList(
    tags$div(class = "subsection-label", text),
    if (!is.null(help)) tags$p(class = "subsection-help", help)
  )
}

step_title <- function(number, title, hint_id) {
  tags$div(
    class = "step-head",
    tags$span(class = "step-num", number),
    tags$div(
      class = "step-head-text",
      tags$span(class = "step-title", title),
      uiOutput(hint_id, inline = TRUE)
    )
  )
}

# Grouping and status summaries ------------------------------------------------

build_group_summary_ui <- function(summary_stats) {
  groups <- summary_stats$groups

  if (is.null(groups) || nrow(groups) == 0) {
    return(
      div(
        class = "group-preview",
        tags$details(
          class = "group-preview-details",
          tags$summary(
            class = "group-preview-summary",
            tags$span(class = "group-preview-title", "Grouping preview"),
            tags$span(class = "group-preview-meta", "no matches")
          ),
          div(
            class = "group-preview-body",
            tags$p(
              class = "helper-text helper-text--compact",
              "This grouping file does not map any taxa in the current dataset."
            )
          )
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
    class = "group-preview",
    tags$details(
      class = "group-preview-details",
      tags$summary(
        class = "group-preview-summary",
        tags$span(class = "group-preview-title", "Grouping preview"),
        tags$span(
          class = "group-preview-meta",
          sprintf("%s groups · %s unmatched", nrow(groups), summary_stats$unmatched)
        )
      ),
      div(
        class = "group-preview-body",
        do.call(div, c(list(class = "group-summary-list"), group_items)),
        note
      )
    )
  )
}

build_selection_summary_ui <- function(
  data_type,
  site_count,
  record_count,
  polygon_count,
  grouping_label,
  stale,
  has_results,
  loading = FALSE,
  site_count_label = NULL
) {
  status <- if (loading) {
    status_pill("Loading BIAD data ...", "loading")
  } else if (polygon_count == 0) {
    status_pill("Add at least one polygon to define a study area.", "warning")
  } else if (has_results && stale) {
    status_pill("Inputs changed — rerun the analysis to refresh the charts.", "stale")
  } else if (has_results) {
    status_pill("Results are in sync with the current map and controls.", "fresh")
  } else {
    status_pill("Ready to run once the study area looks right.", "neutral")
  }

  metrics <- metrics_grid(
    metric_chip("Dataset", data_type, "database"),
    metric_chip("Sites", site_count_label %||% format(site_count, big.mark = ","), "location-dot"),
    metric_chip("Records", format(record_count, big.mark = ","), "layer-group"),
    metric_chip("Polygons", polygon_count, "draw-polygon"),
    metric_chip("Grouping", grouping_label, "sitemap")
  )

  div(class = "summary-block", status, metrics)
}

build_analysis_status_ui <- function(result, stale) {
  if (is.null(result)) {
    return(
      div(
        class = "summary-block",
        status_pill("No analysis has been run yet.", "neutral"),
        helper_text(
          "Choose a dataset, define polygons, select an analysis method, ",
          "then click Run analysis to compute the ordination."
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
    metric_chip("Method", result$summary$analysis_method, "chart-line"),
    metric_chip("Used polygons", result$summary$polygon_count, "draw-polygon"),
    metric_chip("Phases", result$summary$phase_count, "clock"),
    metric_chip("Taxon groups", result$summary$taxon_group_count, "tags"),
    metric_chip("Records analysed", format(result$summary$records, big.mark = ","), "list-check")
  )

  notes <- NULL
  if (length(result$notes) > 0) {
    notes <- do.call(
      tags$ul,
      c(list(class = "note-list"), lapply(result$notes, tags$li))
    )
  }

  div(class = "summary-block", status, metrics, notes)
}

build_overview_ui <- function(result, stale) {
  if (is.null(result)) {
    return(empty_state(
      "No results yet",
      "Run an analysis to see a summary of the ordination, notes, and export contents.",
      "flask"
    ))
  }

  div(
    class = "overview-stack",
    build_analysis_status_ui(result, stale)
  )
}

# Main UI ----------------------------------------------------------------------

app_ui <- function() {
  faunal_groupings <- available_groupings("Faunal")
  faunal_default <- preferred_grouping_path("Faunal")

  # The non-Combined selector initially shows the Faunal schemes. The server
  # swaps the choices when the dataset selector changes.
  bslib::page_fluid(
    title = "BIAD Taxa Explorer",
    theme = app_theme(),
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "app.css"),
      tags$link(rel = "icon", type = "image/png", href = "biad.logo.png"),
      tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
      tags$meta(name = "color-scheme", content = "light")
    ),
    div(
      class = "app-shell",

      # Header ----------------------------------------------------------------
      tags$header(
        class = "app-header",
        div(
          class = "app-header-brand",
          tags$a(
            href = "https://biadwiki.org",
            target = "_blank",
            rel = "noopener noreferrer",
            class = "brand-logo-link",
            tags$img(src = "biad.logo.png", alt = "BIAD logo", class = "brand-logo")
          ),
          div(
            class = "brand-text",
            tags$h1("BIAD Taxa Explorer"),
            tags$p(
              "Faunal and botanical taxa in custom study areas, grouped through ",
              "published schemes and compared with CA or PCA."
            )
          )
        ),
        div(
          class = "app-header-actions",
          actionButton(
            "refresh_data",
            tagList(app_icon("arrows-rotate"), "Refresh data"),
            class = "btn btn-sm btn-outline-secondary"
          ),
          actionButton(
            "reset_app",
            tagList(app_icon("rotate-left"), "Reset"),
            class = "btn btn-sm btn-outline-secondary"
          )
        )
      ),

      uiOutput("data_status"),

      bslib::layout_columns(
        col_widths = c(4, 8),
        gap = "1.1rem",

        # Sidebar -------------------------------------------------------------
        div(
          class = "sidebar-stack",
          div(
            class = "sidebar-panel",
          bslib::accordion(
            id = "workflow_accordion",
            open = c("panel-data", "panel-time"),
            multiple = TRUE,
            bslib::accordion_panel(
              title = step_title("1", "Data and study area", "step1_hint"),
              value = "panel-data",
              selectInput(
                "data_type_selector",
                "Dataset",
                choices = c("Faunal" = "Faunal", "Botanical" = "Botanical", "Combined" = "Combined"),
                selected = "Faunal",
                width = "100%"
              ),
              conditionalPanel(
                condition = "input.data_type_selector != 'Combined'",
                selectInput(
                  "file_selector_tx",
                  "Taxon grouping",
                  choices = faunal_groupings,
                  selected = faunal_default,
                  width = "100%"
                )
              ),
              conditionalPanel(
                condition = "input.data_type_selector == 'Combined'",
                selectInput(
                  "file_selector_tx_combined_faunal",
                  "Faunal grouping",
                  choices = faunal_groupings,
                  selected = faunal_default,
                  width = "100%"
                ),
                uiOutput("combined_botanical_groups_ui"),
                helper_text(
                  "Faunal groups are kept and the chosen botanical subgroups are added. ",
                  "Taxa used by several subgroups are counted in each and flagged after the analysis."
                )
              ),
              div(class = "sidebar-divider"),
              subsection_label(
                "Study area",
                "Upload a GeoPackage or draw polygons directly on the map."
              ),
              fileInput(
                "shapefile",
                NULL,
                accept = c(".gpkg", ".geojson", ".json"),
                buttonLabel = tagList(app_icon("upload"), "Upload GeoPackage"),
                placeholder = "No file selected",
                width = "100%"
              ),
              div(
                class = "button-row",
                actionButton(
                  "clear_polygons_btn",
                  tagList(app_icon("trash"), "Clear polygons"),
                  class = "btn-sm btn-danger-soft"
                )
              ),
              tags$details(
                class = "manage-details reference-details",
                tags$summary(
                  tags$span("Reference points (CSV)"),
                  uiOutput("reference_meta", inline = TRUE)
                ),
                div(
                  class = "manage-body",
                  fileInput(
                    "reference_csv",
                    NULL,
                    accept = ".csv",
                    buttonLabel = tagList(app_icon("upload"), "Upload CSV"),
                    placeholder = "No file selected",
                    width = "100%"
                  ),
                  uiOutput("reference_columns_ui"),
                  uiOutput("reference_status_ui"),
                  uiOutput("reference_controls_ui"),
                  helper_text(
                    "Coordinates are treated as WGS84 decimal degrees. ",
                    "Longitude and latitude columns are detected from their names; ",
                    "adjust the chooser if needed."
                  )
                )
              ),
              uiOutput("taxon_table"),
              tags$details(
                class = "manage-details",
                tags$summary("Manage grouping files"),
                div(
                  class = "manage-body",
                  fileInput(
                    "upload_group",
                    NULL,
                    accept = ".csv",
                    buttonLabel = tagList(app_icon("upload"), "Upload CSV"),
                    placeholder = "No file selected",
                    width = "100%"
                  ),
                  div(
                    class = "manage-actions",
                    downloadButton("download_group", tagList(app_icon("download"), "Download"), class = "btn-sm btn-soft"),
                    actionButton("edit_group", tagList(app_icon("pen"), "Edit labels"), class = "btn-sm btn-soft"),
                    downloadButton("download_template", tagList(app_icon("file-lines"), "Template"), class = "btn-sm btn-soft")
                  ),
                  helper_text("Format: group number, group label, taxa codes separated by +.")
                )
              )
            ),
            bslib::accordion_panel(
              title = step_title("2", "Time window", "step2_hint"),
              value = "panel-time",
              helper_text("Use the named archaeological periods or build fixed BP slices."),
              selectInput(
                "file_selector_per",
                "Period grouping",
                choices = c(
                  "Automatic BP slices" = "auto",
                  "Named archaeological periods" = "groupings/periods/periods.csv"
                ),
                selected = "auto",
                width = "100%"
              ),
              conditionalPanel(
                condition = "input.file_selector_per == 'auto'",
                selectInput(
                  "time_preset",
                  "Quick range",
                  choices = c(
                    "Custom" = "custom",
                    "8500–3000 BP · 500 yr" = "8500_3000_500",
                    "6000–2000 BP · 500 yr" = "6000_2000_500",
                    "4000–0 BP · 250 yr" = "4000_0_250"
                  ),
                  selected = "custom",
                  width = "100%"
                ),
                div(
                  class = "field-grid",
                  numericInput("start_value", "Start (BP)", value = 8500, min = 0, width = "100%"),
                  numericInput("end_value", "End (BP)", value = 3000, min = 0, width = "100%")
                ),
                sliderInput(
                  "duration",
                  "Time-slice width (years)",
                  min = 100,
                  max = 1000,
                  value = 500,
                  step = 50,
                  width = "100%"
                ),
                uiOutput("time_validation")
              )
            ),
            bslib::accordion_panel(
              title = step_title("3", "Analysis settings", "step3_hint"),
              value = "panel-analysis",
              selectInput(
                "analysis_method",
                "Ordination method",
                choices = c(
                  "Correspondence analysis (CA)" = "ca",
                  "Principal component analysis (PCA)" = "pca"
                ),
                selected = "ca",
                width = "100%"
              ),
              conditionalPanel(
                condition = "input.analysis_method == 'pca'",
                selectInput(
                  "pca_transform",
                  "PCA input matrix",
                  choices = c("log(count + 1)" = "log1p", "Raw counts" = "raw"),
                  selected = "log1p",
                  width = "100%"
                ),
                selectInput(
                  "pca_scaling",
                  "PCA column scaling",
                  choices = c(
                    "Center + scale" = "scale",
                    "Center only" = "center",
                    "No centering or scaling" = "none"
                  ),
                  selected = "scale",
                  width = "100%"
                )
              ),
              uiOutput("method_note")
            )
          ),

          div(
            class = "run-section",
            step_title("4", "Run analysis", "step4_hint"),
            tags$p(class = "run-section-help", "Charts refresh only when you click the button."),
            actionButton(
              "run_btn",
              tagList(app_icon("play"), "Run analysis"),
              class = "btn-primary btn-block btn-run"
            ),
            uiOutput("run_hint")
          )
          )
        ),

        # Main ----------------------------------------------------------------
        div(
          class = "main-stack",
          bslib::card(
            class = "app-card map-card",
            bslib::card_header(
              class = "app-card-header",
              div(
                class = "app-card-title-wrap",
                div(
                  tags$h3(class = "app-card-title", "Map and coverage"),
                  tags$p(class = "app-card-subtitle", "Draw or upload polygons to define the study area.")
                )
              ),
              div(class = "app-card-actions", uiOutput("map_exports"))
            ),
            bslib::card_body(
              class = "app-card-body",
              uiOutput("selection_summary"),
              div(
                class = "map-toolbar",
                actionButton("zoom_to_data", tagList(app_icon("magnifying-glass-plus"), "Zoom to data"), class = "btn-sm btn-soft"),
                actionButton("undo_polygon_btn", tagList(app_icon("rotate-left"), "Undo polygon"), class = "btn-sm btn-soft"),
                div(
                  class = "map-toggle",
                  checkboxInput("only_selected_sites", "Only sites in study areas", FALSE)
                )
              ),
              shinycssloaders::withSpinner(
                leafletOutput("map", height = "60vh"),
                type = 8,
                color = "#2f5d5a",
                size = 0.7
              )
            )
          ),

          div(
            class = "results-wrap",
            bslib::navset_card_tab(
              id = "results_tabs",
              title = "Analysis results",
              full_screen = TRUE,
              bslib::nav_panel(
                "Overview",
                uiOutput("analysis_overview")
              ),
              bslib::nav_panel(
                "Raw counts",
                shinycssloaders::withSpinner(plotOutput("plot1", height = "360px"), type = 8, color = "#2f5d5a", size = 0.7),
                tags$p(class = "plot-caption", "Total counts per period for each selected polygon.")
              ),
              bslib::nav_panel(
                "Composition",
                shinycssloaders::withSpinner(plotOutput("plot2_legend", height = "150px"), type = 8, color = "#2f5d5a", size = 0.6),
                shinycssloaders::withSpinner(plotOutput("plot2", height = "360px"), type = 8, color = "#2f5d5a", size = 0.7),
                tags$p(class = "plot-caption", "Relative share of taxon groups across periods within each polygon."),
                uiOutput("plot_culture_composition_limit_ui"),
                shinycssloaders::withSpinner(plotOutput("plot2_culture", height = "430px"), type = 8, color = "#2f5d5a", size = 0.7),
                tags$p(class = "plot-caption", "Relative share of taxon groups within each culture, aggregated across the current selection.")
              ),
              bslib::nav_panel(
                "Ordination",
                div(
                  class = "label-controls",
                  checkboxInput("plot_label_year", "Year labels", FALSE),
                  checkboxInput("plot_label_taxa", "Taxon labels", FALSE),
                  checkboxInput("plot_label_culture", "Culture labels", FALSE)
                ),
                uiOutput("plot_culture_limit_ui"),
                div(
                  class = "overview-note",
                  tags$strong("Reading the ordination"),
                  tags$p(
                    "Points are phase × polygon combinations. Arrows connect phases in time ",
                    "order inside each polygon, so arrow direction shows how composition ",
                    "changes through time. Use the controls above to toggle year, taxon, and ",
                    "culture labels."
                  )
                ),
                shinycssloaders::withSpinner(plotOutput("plot3", height = "680px"), type = 8, color = "#2f5d5a", size = 0.7),
                tags$p(class = "plot-caption", "Arrows follow phase order inside each polygon."),
                shinycssloaders::withSpinner(plotOutput("plot4_diag", height = "300px"), type = 8, color = "#2f5d5a", size = 0.6),
                tags$p(class = "plot-caption", "Variance explained (PCA) or inertia explained (CA) per dimension."),
                shinycssloaders::withSpinner(plotOutput("plot4", height = "520px"), type = 8, color = "#2f5d5a", size = 0.7),
                tags$p(class = "plot-caption", "How the first two ordination dimensions evolve through time for each polygon.")
              ),
              bslib::nav_panel(
                "Data preview",
                DT::DTOutput("analysis_preview"),
                tags$p(class = "plot-caption", "Preview of the aggregated export table.")
              )
            )
          ),

          tags$footer(
            class = "app-footer",
            tags$span("BIAD Taxa Explorer"),
            tags$span(class = "footer-dot", "·"),
            tags$span("Internal research tool")
          )
        )
      )
    )
  )
}
