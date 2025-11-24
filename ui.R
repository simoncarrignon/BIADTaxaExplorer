library(sf)
library(leaflet)
library(leaflet.extras)
library(shiny)
library(jsonlite)
library(geojsonsf)


ui <- fluidPage(
  tags$div(
    style = "background-color: #f8d7da; color: #721c24; padding: 10px; border-radius: 5px; margin-bottom: 20px; position:fixed;z-index:1000;width:100%;top:0;left:0",
    tags$h2("Warning: This is an experimental website for internal and testing use only.")
  ),
  leafletOutput("map", height='60vh'),
  verbatimTextOutput("groups"),
  actionButton("run_btn", "Run Analysis"),
  actionButton("clear_polygons_btn", "Clear Polygons"),
  # Create a fluidRow for the select input and tableOutput
  fluidRow(
    column(width = 6,
           selectInput("file_selector_tx", "Select File to group taxon:", 
                       choices = list.files(path = "groupings/faunal_taxa/", 
                                            pattern = "group_.*.csv", 
                                            full.names = TRUE),
                       selected = "groupings/faunal_taxa/group_IV.csv"),
             div(style = "font-size: 12px; margin-bottom: 10px;margin-left: 10px; margin-top: -10px;",  # Small text and close position
           htmlOutput("taxon_table")  # Output for the table
           )
    ),
  ),

  selectInput("file_selector_per", "How to group periods", choices = list("Auto"="auto","Period Name"="groupings/periods/periods.csv"),selected="auto",width="150px"),
  conditionalPanel(
    condition = "input.file_selector_per == 'auto'",
    fluidRow(
      column(width = 2, numericInput("start_value", "Start (BP):", value = 8500, width = "80px")),
      column(width = 2, numericInput("end_value", "End (BP):", value = 2500, width = "80px")),
      column(width = 2, sliderInput("duration", "Time slice duration:", min = 100, max = 1000, value = 500, step = 50, width = "100%"))
    )
  ),
HTML("<p><b>Figure 1:</b> <br/> Raw data. Left: total number of NISP for each phase in each selected region. Right: percentage of each group in each phase for each area </p>"),
fluidRow(
  column(4, plotOutput("plot1", height = "300px")),
  column(4, plotOutput("plot2", height = "300px"))
  ),
HTML("<p><b>Figure 2:</b> <br/> Result of the Correspondence Analysis, with arrow representing the order of the phases</p>"),
fluidRow(
  column(width=6, plotOutput("plot3", height = "500px",width="100%"))
),
HTML("<p><b>Figure 3:</b> <br/> Split of the two dimensions of the CA.</p>"),
fluidRow(
  column(width=6, plotOutput("plot4", height = "500px",width="100%"))
  )
)
