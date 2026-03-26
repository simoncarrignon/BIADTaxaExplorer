library(sf)
library(shiny)
library(bslib)
library(leaflet)
library(leaflet.extras)
library(geojsonsf)
library(DT)

devtools::load_all("../BIADconnect/")

source("R/tools.R")
source("R/map_utils.R")
source("R/analysis.R")
source("R/get.faunal.taxa.R")
source("R/get.botanical.taxa.R")
source("R/plots.R")
source("R/ui_components.R")

conn <- BIADconnect::init.conn()

allFTaxa <- getFaunalTaxa(conn = conn)
allBTaxa <- getABotTaxa(conn = conn)

faunalSites <- unique_site_points(allFTaxa)
botanicalSites <- unique_site_points(allBTaxa)

areapal <- c("black", palette.colors(25, "Pastel 1", recycle = TRUE))
