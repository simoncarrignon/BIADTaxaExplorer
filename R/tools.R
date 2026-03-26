split_plus <- function(list_taxa) {
  if (is.na(list_taxa) || !nzchar(list_taxa)) {
    return(character(0))
  }

  trimws(strsplit(list_taxa, "\\+")[[1]])
}

read_taxon_grouping <- function(filename) {
  grouping_data <- read.csv(
    filename,
    skip = 1,
    header = FALSE
  )
  grouping_data <- grouping_data[!is.na(grouping_data[, 1]), , drop = FALSE]

  group_labels <- trimws(gsub("=", "", grouping_data[, 2]))
  taxa_lists <- lapply(grouping_data[, 3], split_plus)
  names(taxa_lists) <- group_labels

  taxa_lists
}

build_group_mapping <- function(grouping_definition) {
  mapping <- c()

  for (group_name in names(grouping_definition)) {
    taxa_codes <- grouping_definition[[group_name]]
    if (length(taxa_codes) == 0) {
      next
    }

    mapped_codes <- rep(group_name, length(taxa_codes))
    names(mapped_codes) <- taxa_codes
    mapping <- c(mapping, mapped_codes)
  }

  mapping
}

groupTaxons <- function(
  dataset,
  filename = "groupings/faunal_taxa/group_IV.csv",
  logging = TRUE
) {
  grouping_definition <- read_taxon_grouping(filename)
  mapping <- build_group_mapping(grouping_definition)

  dataset$new_txgroups <- unname(mapping[dataset$TaxonCode])

  if (logging) {
    matched <- !is.na(dataset$new_txgroups)
    message(
      sprintf(
        "Grouping %s unique taxa into %s categories; %s taxa remain unmatched.",
        length(unique(dataset$TaxonCode)),
        length(grouping_definition),
        length(unique(dataset$TaxonCode[!matched]))
      )
    )
  }

  dataset
}

groupPeriod <- function(original_period, filename = "groupings/periods/periods.csv") {
  periods <- read.csv(filename, stringsAsFactors = FALSE)
  grouping <- strsplit(periods$PeriodID, "\\+")
  names(grouping) <- periods$Period
  mapping <- build_group_mapping(grouping)

  message(
    sprintf(
      "Grouping %s periods into %s larger categories.",
      length(unique(original_period)),
      length(grouping)
    )
  )

  unname(mapping[original_period])
}

format_grouping_label <- function(path) {
  label <- tools::file_path_sans_ext(basename(path))
  label <- gsub("^group_", "Group ", label)
  gsub("_", " ", label)
}

available_groupings <- function(data_type = "Faunal") {
  directory <- if (identical(data_type, "Botanical")) {
    "groupings/botanical_taxa"
  } else {
    "groupings/faunal_taxa"
  }

  files <- sort(
    list.files(
      directory,
      pattern = "^group_.*\\.csv$",
      full.names = TRUE
    )
  )

  stats::setNames(files, vapply(files, format_grouping_label, character(1)))
}

group_summary_stats <- function(dataset) {
  empty_groups <- data.frame(
    group_name = character(0),
    taxon_count = integer(0),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  empty_groups$taxa <- I(vector("list", 0))

  if (!"new_txgroups" %in% names(dataset)) {
    return(list(groups = empty_groups, unmatched = 0))
  }

  matched_rows <- !is.na(dataset$new_txgroups)
  grouped_taxa <- split(
    as.character(dataset$TaxonCode[matched_rows]),
    as.character(dataset$new_txgroups[matched_rows])
  )

  if (length(grouped_taxa) == 0) {
    unmatched_taxa <- length(unique(dataset$TaxonCode[!matched_rows]))
    return(list(groups = empty_groups, unmatched = unmatched_taxa))
  }

  taxa_lists <- lapply(grouped_taxa, function(taxa_codes) {
    sort(unique(as.character(taxa_codes)))
  })
  groups <- data.frame(
    group_name = names(taxa_lists),
    taxon_count = lengths(taxa_lists),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  groups$taxa <- I(unname(taxa_lists))
  groups <- groups[order(-groups$taxon_count, groups$group_name), , drop = FALSE]
  rownames(groups) <- NULL

  unmatched_taxa <- length(unique(dataset$TaxonCode[!matched_rows]))

  list(groups = groups, unmatched = unmatched_taxa)
}

unique_site_points <- function(dataset) {
  dataset[!duplicated(dataset$SiteId), , drop = FALSE]
}

count_column_for_type <- function(data_type) {
  switch(
    data_type,
    "Faunal" = "NISP",
    "Botanical" = "TotalCount",
    stop("Unknown data type selected.")
  )
}

first_intersection_id <- function(intersections) {
  vapply(
    intersections,
    function(indexes) {
      if (length(indexes) == 0) {
        return(NA_integer_)
      }

      as.integer(indexes[[1]])
    },
    integer(1)
  )
}
