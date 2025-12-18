
getABotTaxa <- function(crs=4326,conn=NULL){
    if(is.null(conn))conn <- BIADconnect::init.conn()
    req <- 
        " SELECT 
            a.SampleID,
            a.TaxonCode,
            a.TotalCount,
            p.PhaseId,
            p.Period,
            p.Culture1,
            p.GMM,
            p.GMS,
            s.SiteId,
            s.SiteName,
            s.Country,
            s.Longitude,
            s.Latitude
        FROM 
            ABotSamples a
        JOIN 
            Phases p ON a.PhaseID = p.PhaseID
        JOIN 
            Sites s ON p.SiteID = s.SiteID; "
     alltaxons <- BIADconnect::query.database(conn = conn,sql.command = req)
     alltaxons <- sf::st_as_sf(alltaxons,coords=c("Longitude","Latitude"))
     st_crs(alltaxons) = crs 
     return(alltaxons)
}

