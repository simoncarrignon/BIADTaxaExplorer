library(sf)
library(leaflet)
library(leaflet.extras)
library(shiny)
library(jsonlite)
library(geojsonsf)


ui <- fluidPage(
  leafletOutput("map", height='60vh'),
  verbatimTextOutput("groups"),
  actionButton("run_btn", "Run Analysis"),
  actionButton("clear_polygons_btn", "Clear Polygons"),
  selectInput("file_selector_tx", "Select File to group taxon:", choices = list.files(path = "groupings/faunal_taxa/",pattern = "group_.*.csv",full.names = TRUE),selected="groupings/faunal_taxa/group_I.csv"),
  selectInput("file_selector_per", "How to group periods", choices = list(Auto="auto","Period Name"="groupings/periods/periods.csv",selected="auto"),width="150px"),
  conditionalPanel(
    condition = "input.file_selector_per == 'auto'",
    fluidRow(
      column(width = 2, numericInput("start_value", "Start (BP):", value = 8500, width = "80px")),
      column(width = 2, numericInput("end_value", "End (BP):", value = 2500, width = "80px")),
      column(width = 2, sliderInput("duration", "Time slice duration:", min = 100, max = 1000, value = 500, step = 50, width = "100%"))
    )
  ),
  fluidRow(
    column(4, plotOutput("plot1", height = "300px")),
    column(4, plotOutput("plot2", height = "300px"))
    ),
  fluidRow(
    column(width=6, plotOutput("plot3", height = "500px",width="100%"))
  ),
  fluidRow(
    column(width=6, plotOutput("plot4", height = "500px",width="100%"))
  )
)
