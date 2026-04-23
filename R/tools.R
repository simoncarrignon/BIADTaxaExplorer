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

groupTaxons_from_definition <- function(
  dataset,
  grouping_definition,
  logging = TRUE,
  selected_groups = NULL,
  group_prefix = NULL,
  allow_multiple_matches = FALSE
) {
  if (!is.null(selected_groups)) {
    grouping_definition <- grouping_definition[names(grouping_definition) %in% selected_groups]
  }

  mapping <- build_group_mapping(grouping_definition)
  mapping_by_taxon <- split(unname(mapping), names(mapping))
  matched_groups <- lapply(as.character(dataset$TaxonCode), function(taxon_code) {
    unique(unname(mapping_by_taxon[[taxon_code]]))
  })

  if (allow_multiple_matches) {
    match_counts <- lengths(matched_groups)
    matched_rows <- match_counts > 0

    if (!any(matched_rows)) {
      dataset$new_txgroups <- rep(NA_character_, nrow(dataset))
    } else {
      expanded_row_index <- rep(which(matched_rows), match_counts[matched_rows])
      dataset <- dataset[expanded_row_index, , drop = FALSE]
      dataset$new_txgroups <- unlist(matched_groups[matched_rows], use.names = FALSE)
    }
  } else {
    dataset$new_txgroups <- vapply(
      matched_groups,
      function(group_names) {
        if (!length(group_names)) {
          return(NA_character_)
        }

        group_names[[1]]
      },
      character(1)
    )
  }

  if (!is.null(group_prefix) && nzchar(group_prefix)) {
    matched <- !is.na(dataset$new_txgroups)
    dataset$new_txgroups[matched] <- paste(group_prefix, dataset$new_txgroups[matched], sep = ": ")
  }

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

groupTaxons <- function(
  dataset,
  filename = "groupings/faunal_taxa/group_IV.csv",
  logging = TRUE,
  selected_groups = NULL,
  group_prefix = NULL
) {
  grouping_definition <- read_taxon_grouping(filename)
  groupTaxons_from_definition(
    dataset = dataset,
    grouping_definition = grouping_definition,
    logging = logging,
    selected_groups = selected_groups,
    group_prefix = group_prefix
  )
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

group_catalog <- function(data_type = "Botanical") {
  grouping_files <- available_groupings(data_type)
  catalog_rows <- lapply(seq_along(grouping_files), function(index) {
    scheme_path <- unname(grouping_files)[[index]]
    scheme_label <- names(grouping_files)[[index]]
    grouping_definition <- read_taxon_grouping(scheme_path)

    if (!length(grouping_definition)) {
      return(NULL)
    }

    data.frame(
      option_id = sprintf("%s::%s", scheme_path, names(grouping_definition)),
      scheme_path = scheme_path,
      scheme_label = scheme_label,
      group_name = names(grouping_definition),
      output_label = make.unique(sprintf("%s / %s", scheme_label, names(grouping_definition))),
      stringsAsFactors = FALSE,
      row.names = NULL
    ) |> transform(taxa = I(unname(grouping_definition)))
  })

  catalog_rows <- Filter(Negate(is.null), catalog_rows)
  if (!length(catalog_rows)) {
    return(data.frame(
      option_id = character(0),
      scheme_path = character(0),
      scheme_label = character(0),
      group_name = character(0),
      output_label = character(0),
      stringsAsFactors = FALSE,
      row.names = NULL
    ))
  }

  do.call(rbind, catalog_rows)
}

group_catalog_choices <- function(data_type = "Botanical") {
  catalog <- group_catalog(data_type)
  if (!nrow(catalog)) {
    return(list())
  }

  split(
    stats::setNames(catalog$option_id, catalog$group_name),
    catalog$scheme_label
  )
}

group_definition_overlap_summary <- function(grouping_definition) {
  taxon_lookup <- data.frame(
    taxon_code = unlist(grouping_definition, use.names = FALSE),
    group_name = rep(names(grouping_definition), lengths(grouping_definition)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  taxon_lookup <- unique(taxon_lookup)

  if (!nrow(taxon_lookup)) {
    return(list())
  }

  overlap_counts <- table(taxon_lookup$taxon_code)
  overlapping_taxa <- names(overlap_counts[overlap_counts > 1L])

  stats::setNames(
    lapply(overlapping_taxa, function(taxon_code) {
      unique(as.character(taxon_lookup$group_name[taxon_lookup$taxon_code == taxon_code]))
    }),
    overlapping_taxa
  )
}

group_definition_from_catalog <- function(selection_ids, data_type = "Botanical") {
  catalog <- group_catalog(data_type)
  selected_rows <- catalog[match(selection_ids, catalog$option_id), , drop = FALSE]
  selected_rows <- selected_rows[!is.na(selected_rows$option_id), , drop = FALSE]

  if (!nrow(selected_rows)) {
    return(list())
  }

  grouping_definition <- selected_rows$taxa
  names(grouping_definition) <- selected_rows$output_label
  grouping_definition
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

standardize_taxa_dataset <- function(dataset, count_column, source_label) {
  culture_column <- intersect(c("Culture", "Culture1"), names(dataset))
  culture_values <- if (length(culture_column)) {
    as.character(dataset[[culture_column[[1]]]])
  } else {
    rep(NA_character_, nrow(dataset))
  }
  source_record_id <- if ("FaunalSpeciesID" %in% names(dataset)) {
    paste("Faunal", as.character(dataset$FaunalSpeciesID), sep = ":")
  } else if ("SampleID" %in% names(dataset)) {
    paste("Botanical", as.character(dataset$SampleID), sep = ":")
  } else {
    paste(source_label, seq_len(nrow(dataset)), sep = ":")
  }
  group_values <- if ("new_txgroups" %in% names(dataset)) {
    as.character(dataset$new_txgroups)
  } else {
    rep(NA_character_, nrow(dataset))
  }

  sf::st_sf(
    data.frame(
      TaxonCode = as.character(dataset$TaxonCode),
      PhaseId = as.character(dataset$PhaseId),
      Period = as.character(dataset$Period),
      Culture = culture_values,
      GMM = dataset$GMM,
      GMS = dataset$GMS,
      SiteId = as.character(dataset$SiteId),
      SiteName = as.character(dataset$SiteName),
      Country = as.character(dataset$Country),
      Dataset = source_label,
      SourceRecordId = source_record_id,
      Count = dataset[[count_column]],
      new_txgroups = group_values,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    geometry = sf::st_geometry(dataset),
    crs = sf::st_crs(dataset)
  )
}

combine_taxa_datasets <- function(faunal_dataset, botanical_dataset) {
  rbind(
    standardize_taxa_dataset(faunal_dataset, "NISP", "Faunal"),
    standardize_taxa_dataset(botanical_dataset, "TotalCount", "Botanical")
  )
}

combine_grouped_taxa_datasets <- function(
  faunal_dataset,
  botanical_dataset,
  faunal_grouping_file,
  botanical_group_ids
) {
  if (is.null(faunal_grouping_file) || !nzchar(faunal_grouping_file)) {
    stop("Select a faunal grouping for the Combined dataset.")
  }
  if (is.null(botanical_group_ids) || !length(botanical_group_ids)) {
    stop("Select at least one botanical subgroup for the Combined dataset.")
  }

  botanical_grouping_definition <- group_definition_from_catalog(botanical_group_ids, "Botanical")
  overlapping_taxa <- group_definition_overlap_summary(botanical_grouping_definition)

  faunal_grouped <- groupTaxons(
    dataset = faunal_dataset,
    filename = faunal_grouping_file,
    logging = FALSE,
    group_prefix = "Faunal"
  )
  botanical_grouped <- groupTaxons_from_definition(
    dataset = botanical_dataset,
    grouping_definition = botanical_grouping_definition,
    logging = FALSE,
    group_prefix = "Botanical",
    allow_multiple_matches = TRUE
  )

  combined_dataset <- combine_taxa_datasets(faunal_grouped, botanical_grouped)
  analysis_notes <- character(0)

  if (length(overlapping_taxa)) {
    overlap_messages <- vapply(
      head(names(overlapping_taxa), 5),
      function(taxon_code) {
        sprintf("%s (%s)", taxon_code, paste(overlapping_taxa[[taxon_code]], collapse = ", "))
      },
      character(1)
    )
    extra_message <- if (length(overlapping_taxa) > 5) {
      sprintf(" and %s more", length(overlapping_taxa) - 5)
    } else {
      ""
    }

    analysis_notes <- c(
      analysis_notes,
      sprintf(
        "Selected botanical subgroups overlap on %s taxa. Matching botanical records were counted in each selected subgroup, so some counts are duplicated. Examples: %s%s.",
        length(overlapping_taxa),
        paste(overlap_messages, collapse = "; "),
        extra_message
      )
    )
  }

  attr(combined_dataset, "analysis_notes") <- analysis_notes
  combined_dataset
}

combine_site_points <- function(faunal_sites, botanical_sites) {
  standardize_sites <- function(dataset) {
    if (is.null(dataset) || !nrow(dataset)) {
      return(NULL)
    }

    sf::st_sf(
      data.frame(
        SiteId = as.character(dataset$SiteId),
        SiteName = as.character(dataset$SiteName),
        Country = as.character(dataset$Country),
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),
      geometry = sf::st_geometry(dataset),
      crs = sf::st_crs(dataset)
    )
  }

  site_frames <- Filter(
    Negate(is.null),
    list(standardize_sites(faunal_sites), standardize_sites(botanical_sites))
  )

  if (!length(site_frames)) {
    return(NULL)
  }

  combined_sites <- do.call(rbind, site_frames)
  unique_site_points(combined_sites)
}

count_column_for_type <- function(data_type) {
  switch(
    data_type,
    "Faunal" = "NISP",
    "Botanical" = "TotalCount",
    "Combined" = "Count",
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
