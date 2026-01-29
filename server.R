library(leaflet)
library(leaflet.extras)
library(shiny)
library(jsonlite)
library(geojsonsf)

split_periods <- function(dataset,selection,start,end,len){
     if(selection == "auto"){
         abs_periods=seq(start,end,-(len))
         new_periods <- cut(dataset$GMM,breaks=abs_periods,label=rev(abs_periods[-length(abs_periods)]))
     }
     else if(selection == "groupings/periods/periods.csv"){
         new_periods <- groupPeriod(dataset$Period,selection)
         new_periods <- factor(new_periods,levels=c("Neolithic","Eneolithic","Bronze Age","Early Iron Age"),ordered=T)
     }
     return(new_periods)
}

shinyServer(function(input, output, session) {
  print(Sys.time())
  check.conn(conn)
  taxon_dataset <- reactiveVal()
  polygons <- reactiveVal()
  combinedData <- reactiveVal()  # Define a reactiveVal to store the combined data

  groupInfo <- reactiveValues(ngroup = length(polygons), group = rep(areapal[1], length(points_sf)))
  taxon_dataset(allBTaxa)  # You should load the initial data here

  observeEvent(input$data_type_selector, {
    print("changing dataset")
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
     dataset$new_periods  <- split_periods( dataset=dataset,selection=input$file_selector_per,start=input$start_value,end=input$end_value,len=input$duration)
    head(dataset)
    taxon_dataset(dataset)
    output$plot1 <- renderPlot({})
    output$plot2 <- renderPlot({})
    output$plot3 <- renderPlot({})
    output$plot4 <- renderPlot({})
  })
    
  output$map <- renderLeaflet({

      polygons <- polygons()
      maps  <- leaflet() |>
        addTiles() |>
        addCircleMarkers(data = points_sf, radius = 3, color = groupInfo$group, group = "sites") 
      if(!is.null(polygons) && length(polygons)>0){
          maps  <- maps |>  
          addPolygons(data = polygons, fillColor = "transparent", stroke = TRUE) 
      }
      maps <- maps |>
        addDrawToolbar(
           targetGroup = 'polygons',
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

      polygons <- polygons()
      if (!is.null(polygons) && length(polygons) > 0) {
        st_write(polygons, file, driver = "gpkg")
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

  observeEvent(input$upload_own_grouping, {
    # Trigger file input
    shinyjs::reset("file_upload") # Reset file input
    updateTabsetPanel(session, "sidebar", "tab_upload") # Ensure file input is visible to the user
  })
  observeEvent(input$file_selector_tx, {
    handleFileSelection()
  })

  observeEvent(input$file_selector_per, {
     print("updating dataset period")
     taxon_new <- taxon_dataset()
     taxon_new$new_periods  <- split_periods( dataset=taxon_new,selection=input$file_selector_per,start=input$start_value,end=input$end_value,len=input$duration)
     taxon_dataset(taxon_new)  
  })

  #maybe merge the one above and the one below?
  observeEvent({
      input$start_value
      input$end_value
      input$duration
  }, {
      taxon_new <- taxon_dataset()
      print("update periods")
     taxon_new$new_periods  <- split_periods( dataset=taxon_new,selection="auto",start=input$start_value,end=input$end_value,len=input$duration)
      taxon_dataset(taxon_new)
      output$plot1 <- renderPlot({})
      output$plot2 <- renderPlot({})
      output$plot3 <- renderPlot({})
      output$plot4 <- renderPlot({})
  })

  observeEvent(input$map_draw_new_feature, {
    feature <- input$map_draw_new_feature
    pol <- geojson_sf(jsonlite::toJSON(feature,auto_unbox=T))
    prev_pol <- polygons()
    polygons <- rbind(prev_pol, pol)
    polygons(polygons)
    groupInfo$ngroup <- groupInfo$ngroup + 1
    groupInfo$group <- ifelse(lengths(st_intersects(points_sf,pol)),areapal[groupInfo$ngroup],groupInfo$group)

    print("---- grouped poly:")
    print(polygons)
    print("----")

    leafletProxy("map") |>
      clearGroup("sites") |>
      addPolygons(data = polygons, fillColor = "transparent", stroke = TRUE) |>
      addCircleMarkers(data = points_sf, radius = 3, color = groupInfo$group, group = "sites")
  })
  
  
  observeEvent(input$run_btn, { 
     polygons <- polygons()
     if(is.null(polygons) || length(polygons)==0)
     {
         showNotification("Please draw or upload at least one area before proceeding.", type = "warning");
         return()
     }
     #saveRDS(drawnPolygons$polygons, file = "groupings/spatial/groups.RDS");
     groups <- length(polygons) #that should be a multipolygon with ach area manually selected
     print("---start analysis")
     cur_data <- taxon_dataset()
     st_crs(cur_data) <- st_crs(polygons)
     cur_data$new_area  <- as.numeric(st_intersects(cur_data,polygons))
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
         nullrow <- apply(cts,1,sum) == 0
         cts <- cts[!nullrow,,drop=FALSE] 
     }
     cts.ca <- FactoMineR::CA(cts,graph=F)
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
         plot2dim(cts.ca,groupInfo$ngroup)
     })
     grp=as.numeric(st_intersects(points_sf,polygons))
     grp=ifelse(is.na(grp),"black",areapal[grp+1])
     groupInfo$group <- grp
     leafletProxy("map") |>
       clearGroup("sites") |>
       addPolygons(data =  polygons, fillColor = "transparent", stroke = TRUE) |>
       addCircleMarkers(data = points_sf, radius = 3, color = groupInfo$group, group = "sites")
  })


  observeEvent(input$clear_polygons_btn, {
    polygons(NULL)
    groupInfo$ngroup <- 1
    groupInfo$group <- 1
    leafletProxy("map") |>
      clearMarkers()  |>
      addCircleMarkers(data = points_sf, radius = 3, color = "black", group = "sites") 

  })

  handleFileSelection <- function() {
    selected_file <- if (!is.null(input$file_upload)) {
      input$file_upload$datapath
    } else {
      input$file_selector_tx
    }
    
    taxon_new <- taxon_dataset()
    taxon_new <- groupTaxons(taxon_new, selected_file)
    taxon_dataset(taxon_new)  # Update the reactive value
    output$taxon_table <- renderText({
      ns <- lengths(tapply(taxon_new$TaxonCode, taxon_new$new_txgroups, unique))
      paste("\t\t<i> groups are: ", paste0(paste(names(ns), ":", ns, "tx"), collapse = ", "), "</i>  ")
    })
  }
  

  observeEvent(input$shapefile, {
    print("--upload shp--")
    req(input$shapefile)
    uploaded_polygons <- st_read(input$shapefile$datapath)
    print("--upload poly--")
    print(uploaded_polygons)
    polygons(uploaded_polygons)
    grp=as.numeric(st_intersects(points_sf,uploaded_polygons))
    grp=ifelse(is.na(grp),"black",areapal[grp+1])
    groupInfo$group <- grp
    groupInfo$ngroup <- length(uploaded_polygons)
  })


})
