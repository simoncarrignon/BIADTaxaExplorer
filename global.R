require(sf)
devtools::load_all("../BIADconnect/")
source("R/tools.R")
source("R/get.faunal.taxa.R")
source("R/get.botanical.taxa.R")
source("R/plots.R")


conn <- BIADconnect::init.conn()

area_of_interest  <- c( "France","Germany","Netherlands","Sweden", "Great Britain","Belgium","Luxembourg","Switzerland", "Ireland","Denmark","Czechia / Czech Republic","Poland", "Austria" ,"Andorra","Isle of Man", "Slovakia","Liechtenstein","Serbia","Croatia", "Hungary","Romania","Bosnia and Herzegovina","Bulgaria", "Latvia","North Macedonia","Ukraine","Slovenia", "Montenegro", "Guernsey","Russia","Estonia","Lithuania", "Albania","Kosovo","Moldova", "Norway","Finland", "Belarus")

alltaxons <- getFaunalTaxa()
taxons_input  <-  alltaxons[ alltaxons$Country %in%  area_of_interest,]

points_sf <- st_sfc(unique(st_geometry(alltaxons)))
st_crs(points_sf)  <-  st_crs(taxons_input)

areapal <- c("black",palette.colors(25,"Pastel 1",recycle=T))

if (file.exists("groupings/spatial/groups.RDS")) {
    previous_polygon <- readRDS(file = "groupings/spatial/groups.RDS")
} else {
    previous_polygon <- list()
}
