# Lightweight in-memory caches -------------------------------------------------

.app_cache <- new.env(parent = emptyenv())

cache_get <- function(key) {
  if (exists(key, envir = .app_cache, inherits = FALSE)) {
    get(key, envir = .app_cache, inherits = FALSE)
  } else {
    NULL
  }
}

cache_set <- function(key, value) {
  assign(key, value, envir = .app_cache)
  invisible(value)
}

cache_remove <- function(key) {
  if (exists(key, envir = .app_cache, inherits = FALSE)) {
    rm(list = key, envir = .app_cache)
  }
  invisible(NULL)
}

cache_clear <- function(pattern = NULL) {
  keys <- ls(envir = .app_cache, all.names = TRUE)
  if (!is.null(pattern)) {
    keys <- grep(pattern, keys, value = TRUE)
  }
  if (length(keys)) {
    rm(list = keys, envir = .app_cache)
  }
  invisible(NULL)
}

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x)) y else x
}
