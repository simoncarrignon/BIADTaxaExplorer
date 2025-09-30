
getFaunalTaxa <- function(crs=4326){
    conn <- BIADconnect::init.conn()
    req <- 
        " SELECT 
            f.FaunalSpeciesID,
            f.NISP,
            f.TaxonCode,
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
            FaunalSpecies f
        JOIN 
            Phases p ON f.PhaseID = p.PhaseID
        JOIN 
            Sites s ON p.SiteID = s.SiteID; "
     alltaxons <- BIADconnect::query.database(conn = conn,sql.command = req)
     alltaxons <- sf::st_as_sf(alltaxons,coords=c("Longitude","Latitude"))
     st_crs(alltaxons) = crs 
     return(alltaxons)
}

