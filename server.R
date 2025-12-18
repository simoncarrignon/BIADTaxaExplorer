library(leaflet)
library(leaflet.extras)
library(shiny)
library(jsonlite)
library(geojsonsf)

shinyServer(function(input, output, session) {
  print(Sys.time())
  check.conn(conn)
  taxon_dataset <- reactiveVal()
  combinedData <- reactiveVal()  # Define a reactiveVal to store the combined data

  drawnPolygons <- reactiveValues(polygons = previous_polygon)
  groupInfo <- reactiveValues(ngroup = length(previous_polygon), group = rep(areapal[1], length(points_sf)))
  taxon_dataset(allBTaxa)  # You should load the initial data here

  observeEvent(input$data_type_selector, {
    sel <- input$data_type_selector
    type <- ifelse(sel == "Faunal", "faunal_taxa/", "botanical_taxa/")
    path <- file.path("groupings",type)
    updateSelectInput(
      session, "file_selector_tx",
      choices = list.files(path , pattern = "group_.*.csv", full.names = TRUE),
    )
    if( sel == "Faunal")
        dataset <- allFTaxa
    else
        dataset <- allBTaxa
    head(dataset)
    taxon_dataset(dataset)
     output$plot1 <- renderPlot({})
     output$plot2 <- renderPlot({})
     output$plot3 <- renderPlot({})
     output$plot4 <- renderPlot({})
  })
    
  output$map <- renderLeaflet({
    maps  <- leaflet() |>
      addTiles() |>
      addCircleMarkers(data = points_sf, radius = 3, color = groupInfo$group, group = "sites") 
     
    if(!is.null(drawnPolygons$polygons) && length(drawnPolygons$polygons)>0){
      maps  <- maps |>  
        addPolygons(data = do.call(rbind, drawnPolygons$polygons), fillColor = "transparent", stroke = TRUE) 
    }
    maps <- maps |>
      addDrawToolbar(
        targetGroup = 'drawn',
        polygonOptions = drawPolygonOptions(),
        editOptions = F,
        polylineOptions = F, rectangleOptions = F, circleOptions = F,
          markerOptions = F, circleMarkerOptions = F
      )
    maps
  })

  output$download_polygons <- downloadHandler(
    filename = "exported_polygons.gpkg", 
    content = function(file) {
      if (!is.null(drawnPolygons$polygons) && length(drawnPolygons$polygons) > 0) {
        all_polygons <- do.call(rbind, drawnPolygons$polygons)
        st_write(all_polygons, file, driver = "gpkg")
      } else {
        showNotification("No polygons to download!", type = "warning")
      }
    }
  )

  output$download_table <- downloadHandler(
    filename = "aggregated_data.csv", 
    content = function(file) {
      data <- combinedData()
      if (!is.null(data)) {
        write.csv(data, file, row.names=F)
      } else {
        showNotification("No data available for download.", type = "warning",duration=5)
        stop("No data available for download.")
      }
    }
  )

  
  observeEvent(input$file_selector_tx, {
    selected_file <- input$file_selector_tx
    taxon_new <-taxon_dataset()

    taxon_new <- groupTaxons(taxon_new,selected_file)
    taxon_dataset(taxon_new)  # Update the reactive value
    output$taxon_table <- renderText({
         ns=lengths(tapply(taxon_new$TaxonCode,taxon_new$new_txgroups,unique))
         paste("\t\t<i> groups are: ", paste0(paste(names(ns),":",ns,"tx"),collapse=", "),"</i>  ")
        })

  })
  observeEvent(input$file_selector_per, {
     taxon_new <- taxon_dataset()
     if(input$file_selector_per == "auto"){
         abs_periods=seq(input$start_value,input$end_value,-(input$duration))
         taxon_new$new_periods <- cut(taxon_new$GMM,breaks=abs_periods,label=rev(abs_periods[-length(abs_periods)]))
     }
     else if(input$file_selector_per == "groupings/periods/periods.csv"){
         taxon_new$new_periods <- groupPeriod(taxon_new$Period,input$file_selector_per)
         taxon_new$new_periods <- factor(taxon_new$new_periods,levels=c("Neolithic","Eneolithic","Bronze Age","Early Iron Age"),ordered=T)
     }
     taxon_dataset(taxon_new)  
  })

  observeEvent(input$map_draw_new_feature, {
    feature <- input$map_draw_new_feature
    pol <- geojson_sf(jsonlite::toJSON(feature,auto_unbox=T))
    drawnPolygons$polygons <- append(drawnPolygons$polygons, list(pol))
    groupInfo$ngroup <- groupInfo$ngroup+1
    groupInfo$group <- ifelse(lengths(st_intersects(points_sf,pol)),areapal[groupInfo$ngroup],groupInfo$group)

    leafletProxy("map") |>
      clearGroup("sites") |>
      addPolygons(data = do.call(rbind, drawnPolygons$polygons), fillColor = "transparent", stroke = TRUE) |>
      addCircleMarkers(data = points_sf, radius = 3, color = groupInfo$group, group = "sites")
  })
  
  
  observeEvent(input$run_btn, { 
     if(is.null(drawnPolygons$polygons) || length(drawnPolygons$polygons)==0)
     {
         showNotification("Please draw at least one area before proceeding.", type = "warning");
         return()
     }
     saveRDS(drawnPolygons$polygons, file = "groupings/spatial/groups.RDS");
     groups <- do.call("rbind",drawnPolygons$polygons) #that should be a multipolygon with ach area manually selected
     cur_data=taxon_dataset()
     st_crs(cur_data)=st_crs(groups)
     cur_data$new_area  <- as.numeric(st_intersects(cur_data,groups))
     subregions <- cur_data[!is.na(cur_data$new_area) & !is.na(cur_data$new_periods),]
     subregions$phase <- paste(subregions$new_periods,subregions$new_area,sep="-")
     sel <- input$data_type_selector
     if(sel == "Faunal") cnt <- "NISP"
     if(sel == "Botanical") cnt <- "TotalCount"

     print("running CA")
     cts <- tapply(subregions[[cnt]],list(subregions$phase,subregions$new_txgroups),sum,na.rm=T)

     if(any(is.na(cts)))
     { 
         showNotification("Some taxon are missing from some periods/area, we will assume they are asbent (ie replace NA by 0) of the group", type = "warning")
         cts[is.na(cts)] <- 0
     }
     cts.ca <<- FactoMineR::CA(cts,graph=F)
     an <- rownames(cts)
     perarea <- do.call("rbind",strsplit(an,"-"))
     colnames(perarea) <- c("Phase start","Polygon ID")
     combined_data  <- cbind.data.frame(perarea,cts,cts.ca$row$coord)
     combinedData(combined_data)
     output$plot1 <- renderPlot({
         countTotal(subregions,cnt)
     })

     output$plot2 <- renderPlot({
         bySpeciesComposition(subregions,cnt)
     })
     output$plot3 <- renderPlot({

         if(!(is.null(dim(cts.ca$row$coord)))){
             plotCAarrows(cts.ca)
         }

     })
     output$plot4 <- renderPlot({
         print(groupInfo$ngroup)
         plot2dim(cts.ca,groupInfo$ngroup)
         saveRDS(file="cts.ca.RDS",cts.ca)
     })
     grp=as.numeric(st_intersects(points_sf,do.call("rbind",drawnPolygons$polygons)))
     grp=ifelse(is.na(grp),"black",areapal[grp+1])
     groupInfo$group <- grp
     leafletProxy("map") |>
       clearGroup("sites") |>
       addPolygons(data = do.call(rbind, drawnPolygons$polygons), fillColor = "transparent", stroke = TRUE) |>
       addCircleMarkers(data = points_sf, radius = 3, color = groupInfo$group, group = "sites")
  })


  observeEvent(input$clear_polygons_btn, {
    drawnPolygons$polygons <- list()  # Clear the stored polygons
    groupInfo$ngroup <- 1
    groupInfo$group <- 1
    
    leafletProxy("map") |>
      clearMarkers()  |>
      addCircleMarkers(data = points_sf, radius = 3, color = "black", group = "sites") 

  })


})
