library(sf)
library(shiny)
library(bslib)
library(leaflet)
library(leaflet.extras)
library(geojsonsf)

devtools::load_all("../BIADconnect/")

source("R/tools.R")
source("R/map_utils.R")
source("R/analysis.R")
source("R/get.faunal.taxa.R")
source("R/get.botanical.taxa.R")
source("R/plots.R")
source("R/ui_components.R")

conn <- BIADconnect::init.conn()

area_of_interest <- c(
  "France", "Germany", "Netherlands", "Sweden", "Great Britain", "Belgium",
  "Luxembourg", "Switzerland", "Ireland", "Denmark", "Czechia / Czech Republic",
  "Poland", "Austria", "Andorra", "Isle of Man", "Slovakia", "Liechtenstein",
  "Serbia", "Croatia", "Hungary", "Romania", "Bosnia and Herzegovina", "Bulgaria",
  "Latvia", "North Macedonia", "Ukraine", "Slovenia", "Montenegro", "Guernsey",
  "Russia", "Estonia", "Lithuania", "Albania", "Kosovo", "Moldova", "Norway",
  "Finland", "Belarus"
)

allFTaxa <- getFaunalTaxa(conn = conn)
allBTaxa <- getABotTaxa(conn = conn)
allFTaxa <- allFTaxa[allFTaxa$Country %in% area_of_interest, ]
allBTaxa <- allBTaxa[allBTaxa$Country %in% area_of_interest, ]

faunalSites <- unique_site_points(allFTaxa)
botanicalSites <- unique_site_points(allBTaxa)

areapal <- c("black", palette.colors(25, "Pastel 1", recycle = TRUE))
