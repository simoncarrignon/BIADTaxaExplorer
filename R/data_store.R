# Cached access to the BIAD database ------------------------------------------

biad_cache_file <- function() {
  file.path("cache", "biad_data.rds")
}

biad_cache_ttl <- function() {
  hours <- suppressWarnings(as.numeric(Sys.getenv("BIAD_CACHE_HOURS", "12")))
  if (!is.finite(hours) || hours <= 0) hours <- 12
  hours * 3600
}

biad_cache_age <- function(path = biad_cache_file()) {
  if (!file.exists(path)) return(Inf)
  as.numeric(difftime(Sys.time(), file.info(path)$mtime, units = "secs"))
}

biad_cache_fresh <- function(path = biad_cache_file()) {
  is.finite(biad_cache_age(path)) && biad_cache_age(path) < biad_cache_ttl()
}

format_cache_age <- function(seconds) {
  if (!is.finite(seconds)) return("no cache")
  if (seconds < 90) return("just now")
  if (seconds < 3600) return(sprintf("%d min ago", round(seconds / 60)))
  if (seconds < 86400) return(sprintf("%.1f h ago", seconds / 3600))
  sprintf("%.1f days ago", seconds / 86400)
}

clear_biad_cache <- function() {
  path <- biad_cache_file()
  if (file.exists(path)) unlink(path)
  invisible(NULL)
}

load_biad_data <- function(refresh = FALSE, quiet = FALSE) {
  path <- biad_cache_file()

  if (!refresh && biad_cache_fresh(path)) {
    cached <- tryCatch(readRDS(path), error = function(e) NULL)
    if (!is.null(cached) && all(c("faunal", "botanical") %in% names(cached))) {
      cached$source <- "cache"
      cached$cache_age <- biad_cache_age(path)
      return(cached)
    }
  }

  if (!quiet) message("Querying the BIAD database ...")
  conn <- BIADconnect::init.conn()
  on.exit(try(BIADconnect::disconnect(), silent = TRUE), add = TRUE)

  faunal <- getFaunalTaxa(conn = conn)
  botanical <- getABotTaxa(conn = conn)

  bundle <- list(
    faunal = faunal,
    botanical = botanical,
    faunal_sites = unique_site_points(faunal),
    botanical_sites = unique_site_points(botanical),
    loaded_at = Sys.time(),
    source = "database",
    cache_age = 0
  )

  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  tryCatch(
    saveRDS(bundle, path),
    error = function(e) warning("Could not write the data cache: ", conditionMessage(e))
  )

  bundle
}

biad_data_summary <- function(bundle) {
  if (is.null(bundle)) {
    return(list(
      faunal_records = 0L,
      botanical_records = 0L,
      faunal_sites = 0L,
      botanical_sites = 0L,
      loaded_at = NULL,
      source = "none",
      cache_age = Inf
    ))
  }

  list(
    faunal_records = nrow(bundle$faunal),
    botanical_records = nrow(bundle$botanical),
    faunal_sites = nrow(bundle$faunal_sites),
    botanical_sites = nrow(bundle$botanical_sites),
    loaded_at = bundle$loaded_at,
    source = bundle$source,
    cache_age = bundle$cache_age
  )
}
