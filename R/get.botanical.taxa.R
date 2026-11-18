# Botanical records ------------------------------------------------------------

getABotTaxa <- function(crs = 4326, conn = NULL) {
  if (is.null(conn)) conn <- BIADconnect::init.conn()

  sql <- "
    SELECT
      a.SampleID,
      a.TaxonCode,
      a.TotalCount,
      p.PhaseId,
      p.Period,
      p.Culture1 AS Culture,
      p.GMM,
      p.GMS,
      s.SiteId,
      s.SiteName,
      s.Country,
      s.Longitude,
      s.Latitude
    FROM ABotSamples a
    JOIN Phases p ON a.PhaseID = p.PhaseID
    JOIN Sites s ON p.SiteID = s.SiteID;"

  taxa <- BIADconnect::query.database(conn = conn, sql.command = sql)
  taxa <- sf::st_as_sf(taxa, coords = c("Longitude", "Latitude"))
  sf::st_crs(taxa) <- crs
  taxa
}
