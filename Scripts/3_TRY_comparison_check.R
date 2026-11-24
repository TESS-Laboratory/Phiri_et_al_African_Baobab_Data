# ------------------------------------------------------------------------------
# Introduction
# ------------------------------------------------------------------------------

# Compares every Adansonia observation in a TRY release with the baobab
# records of the public BAM data file written by 2_QAQC_BAM_Data.qmd. A TRY
# observation within 50 m of a BAM baobab that also agrees in height or
# contributor is classed as a probable duplicate. The Adansonia rows of the
# release are saved next to it, so later runs do not read the large file
# again.


# 1. Pipeline folder and paths

# rstudioapi is optional: it finds the folder of this script when run in
# RStudio.
if (requireNamespace("rstudioapi", quietly = TRUE)) library(rstudioapi)

script_folders <- c(
  tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) NA_character_),
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    tryCatch(dirname(rstudioapi::getSourceEditorContext()$path), error = function(e) NA_character_)
  } else NA_character_,
  getwd())
script_folders <- script_folders[!is.na(script_folders) & file.exists(file.path(script_folders, "pipeline_paths.csv"))]
if (length(script_folders) == 0) {
  stop("pipeline_paths.csv was not found. Run 1_Harmonise_tree_data.R first, then open this script ",
       "from the same folder and click Source.", call. = FALSE)
}
setwd(script_folders[1])
pipeline_paths <- utils::read.csv("pipeline_paths.csv", colClasses = "character", na.strings = character())

# pipeline_value() returns one path from pipeline_paths.csv, where every path
# of the pipeline is set by 1_Harmonise_tree_data.R.
pipeline_value <- function(name) {
  value <- pipeline_paths$value[pipeline_paths$name == name]
  if (length(value) == 0) stop(name, " is missing from pipeline_paths.csv. Run 1_Harmonise_tree_data.R again.",
                               call. = FALSE)
  value[1]
}

# Input files
public_data_file <- pipeline_value("public_data_file")   # public BAM data, written by 2_QAQC_BAM_Data.qmd
dictionary_file  <- pipeline_value("dictionary_file")    # BAM_variable_dictionary.csv
try_data_file    <- pipeline_value("try_data_file")      # TRY release (.txt file or its folder)

# Output folders
try_results_folder <- pipeline_value("try_results_folder")   # tables of this script


# Output files
try_observations_output <- file.path(try_results_folder, "try_adansonia_observations.csv")
duplicate_check_output  <- file.path(try_results_folder, "try_bam_duplicate_check.csv")
trait_coverage_output   <- file.path(try_results_folder, "try_trait_coverage.csv")
summary_output          <- file.path(try_results_folder, "try_cross_check_summary.csv")
# Map of the comparison, for checking only; it is not a figure of the
# manuscript.
figure_output_stem      <- file.path(try_results_folder, "try_bam_locations_check")  # .pdf and .png

message("Files used by 3_TRY_comparison_check.R:",
        "\n  public data file:   ", public_data_file,
        "\n  dictionary:         ", dictionary_file,
        "\n  TRY release:        ", try_data_file,
        "\n  results folder:     ", try_results_folder,
        "\n  figure:             ", figure_output_stem, ".pdf/.png")

# Columns are renamed to their dictionary roles, so variables can be renamed
# in the dictionary.
dictionary_roles <- utils::read.csv(dictionary_file, colClasses = "character")[, c("variable_name", "role")]
rename_to_roles <- function(table) {
  hit <- names(table) %in% dictionary_roles$variable_name
  names(table)[hit] <- dictionary_roles$role[match(names(table)[hit], dictionary_roles$variable_name)]
  table
}


# 2. Cross-check settings 

genus_pattern  <- "^Adansonia"
target_species <- "Adansonia digitata"

# Version of the TRY release, as stated by TRY when the data were released
# (for example "6.0"); it is reported with the results and cited in the
# manuscript.
try_release_version <- ""

# TRY identifiers of latitude, longitude, altitude and vegetative plant height.
latitude_data_id  <- 59
longitude_data_id <- 60
altitude_data_id  <- 61
height_trait_id   <- 3106

# Original names of map-grid coordinates, which cannot be converted without
# their zone.
map_grid_names <- c("utm x", "utm y", "utm_x", "utm_y", "easting", "northing")

# TRY ErrorRisk above this value marks a probable outlier.
error_risk_limit <- 4

# Distances for the same tree and the same site; 50 m allows for older GPS
# receivers and rounded coordinates.

same_tree_metres <- 50
same_site_kilometres <- 1

# Heights within this proportion, or 1 m, are taken to agree.
height_tolerance <- 0.10

# BAM variables and the words that identify the matching TRY traits.
trait_matching <- tibble::tribble(
  ~bam_variable,                 ~try_trait_words,
  "height",                      "plant height",
  "circumference",               "stem diameter|circumference|girth|dbh",
  "crown_1",                     "crown|canopy",
  "leaves",                      "leaf phenology|leafing|leaf flushing|deciduous",
  "flowers",                     "flower|flowering",
  "fruit",                       "fruit|fruiting",
  "multiple_stems",              "stem number|multi.?stem|branching"
)


# 3. Packages and helper functions 

required_packages <- c("data.table", "readr", "dplyr", "tidyr", "stringr", "tibble", "ggplot2",
                       "sf", "rnaturalearth", "rnaturalearthdata", "ragg", "systemfonts")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                              logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Please install these packages first: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}
# data.table is loaded first, so that dplyr's first(), last() and between()
# take precedence over the data.table functions of the same names.
suppressPackageStartupMessages({
  library(data.table)         # fast reading of the TRY release
  library(readr)              # reading and writing CSV files
  library(dplyr)              # table manipulation
  library(tidyr)              # reshaping tables
  library(stringr)            # text handling
  library(tibble)             # tables
  library(ggplot2)            # map of the comparison
  library(sf)                 # distances between TRY and BAM trees
  library(rnaturalearth)      # country boundaries
  library(rnaturalearthdata)  # data used by rnaturalearth
  library(ragg)               # PNG output
  library(systemfonts)        # finding Arial or Helvetica
})
dir.create(try_results_folder, recursive = TRUE, showWarnings = FALSE)

to_dictionary_style <- function(names) {
  names %>% tolower() %>% str_replace_all("[^a-z0-9]+", "_") %>% str_remove_all("^_+|_+$")
}

collapse_unique <- function(values) {
  values <- unique(values[!is.na(values) & values != ""])
  if (length(values) == 0) NA_character_ else paste(values, collapse = "; ")
}

given_path <- function(value) is.character(value) && length(value) == 1 && !is.na(value) && value != ""

# TRY coordinates without a standardised value are read from the original
# text with the same rules as step 1.
compass_letter_pattern <- "(?<![A-Z])[NSEW](?![A-Z])"

# standardise_coordinate_text() reduces a coordinate to digits, spaces, signs,
# full stops and compass letters.
standardise_coordinate_text <- function(text) {
  text <- as.character(text)
  unmarked <- !is.na(text) & Encoding(text) == "unknown" & validUTF8(text)
  Encoding(text)[unmarked] <- "UTF-8"
  text <- str_to_upper(enc2utf8(text))
  text <- str_replace_all(text, "[\\u00A0\\u2007\\u202F]", " ")
  text <- str_replace_all(text, "[\\u2212\\u2013\\u2014\\u2010\\u2012]", "-")
  text <- str_replace_all(text, "\\bNORTH(ERN)?\\b", " N ")
  text <- str_replace_all(text, "\\bSOUTH(ERN)?\\b", " S ")
  text <- str_replace_all(text, "\\bEAST(ERN)?\\b",  " E ")
  text <- str_replace_all(text, "\\bWEST(ERN)?\\b",  " W ")
  text <- str_replace_all(text, "\\b(LATITUDE|LONGITUDE|LAT|LONG|LON|LNG|COORDINATES?|GPS|POSITION)\\b\\s*[:=]?", " ")
  text <- str_replace_all(text, "\\b(DEGREES?|DEGS?|MINUTES?|MINS?|SECONDS?|SECS?)\\b", " ")
  text <- str_replace_all(text, "[\\u00B0\\u00BA\\u02DA\\u2218'\\u2032\\u2019\\u2018\\u00B4`\"\\u2033\\u201D\\u201C:*]", " ")
  # A comma is a decimal mark only when the text has no full stop.
  if (!str_detect(text, fixed("."))) text <- str_replace_all(text, "(\\d),(\\d)", "\\1.\\2")
  text <- str_replace_all(text, "[,;]", " ")
  text <- str_replace_all(text, "(?<=\\s|^)\\+(?=\\d)", "")
  text <- str_squish(text)
  # Letters D and M used as degree and minute markers become spaces.
  if (str_detect(text, "\\d\\s*D\\s*\\d|\\d\\s*M\\s*(\\d|$|[NSEW])")) {
    text <- str_replace_all(text, "(?<=\\d)\\s*D\\s*(?=\\d)", " ")
    text <- str_replace_all(text, "(?<=\\d)\\s*M(?=\\s*(\\d|$|[NSEW]))", " ")
    text <- str_replace_all(text, "(?<=\\d)\\s*S(?=\\s+[NSEW](\\s|$))", " ")
    text <- str_squish(text)
  }
  text
}

# read_single_coordinate() reads one latitude or longitude and notes how it was
# written.
read_single_coordinate <- function(text, which_one) {
  answer <- function(value, style, note = NA_character_, compass = NA_character_) {
    list(value = value, style = style, note = note, compass = compass)
  }
  text <- standardise_coordinate_text(text)
  if (text == "") return(answer(NA_real_, NA_character_))
  if (grepl("^-?\\d+(\\.\\d+)?E[+-]?\\d+$", text, perl = TRUE)) {
    return(answer(as.numeric(text), "decimal degrees (scientific notation)"))
  }
  compass <- regmatches(text, gregexpr(compass_letter_pattern, text, perl = TRUE))[[1]]
  text <- gsub("\\s+", " ", trimws(gsub(compass_letter_pattern, " ", text, perl = TRUE)), perl = TRUE)
  if (grepl("[A-Z]", text, perl = TRUE)) {
    return(answer(NA_real_, "not recognised", "contains letters (perhaps a place name or a map grid without its zone)"))
  }
  if (length(unique(compass)) > 1) {
    return(answer(NA_real_, "not recognised", "two compass letters that do not belong together"))
  }
  compass <- if (length(compass) > 0) compass[1] else NA_character_
  minus_sign <- grepl("^-", text, perl = TRUE)
  if (minus_sign && compass %in% c("N", "E")) {
    return(answer(NA_real_, "not recognised", "a minus sign contradicts the compass letter N or E", compass))
  }
  sign <- if (minus_sign || compass %in% c("S", "W")) -1 else 1
  with_letter <- if (!is.na(compass)) " with compass letter" else ""
  text <- sub("^-\\s*", "", text, perl = TRUE)
  numbers_text <- regmatches(text, gregexpr("\\d+(\\.\\d+)?", text, perl = TRUE))[[1]]
  numbers <- as.numeric(numbers_text)
  if (grepl("-", text, perl = TRUE)) {
    return(answer(NA_real_, "not recognised", "a minus sign inside the value", compass))
  }
  if (length(numbers) == 1) {
    value <- numbers[1]
    limit <- if (which_one == "latitude") 90 else 180
    whole_digits <- nchar(sub("\\..*$", "", numbers_text[1], perl = TRUE))
    degree_digits <- if (which_one == "latitude") 2 else 3
    # NMEA output such as 2233.835S, read only when a compass letter is present.
    if (!is.na(compass) && value > limit && whole_digits %in% c(degree_digits + 2, degree_digits + 1) &&
        grepl("\\.", numbers_text[1], perl = TRUE)) {
      degrees <- floor(value / 100)
      minutes <- value - degrees * 100
      if (minutes < 60) {
        return(answer(sign * (degrees + minutes / 60), "NMEA (degrees and decimal minutes run together)",
                      compass = compass))
      }
    }
    # Packed degrees, minutes and seconds, such as 223350S.
    if (!is.na(compass) && value > limit && whole_digits %in% c(degree_digits + 4, degree_digits + 3)) {
      whole   <- floor(value)
      degrees <- floor(whole / 10000)
      minutes <- floor((whole %% 10000) / 100)
      seconds <- (whole %% 100) + (value - whole)
      if (minutes < 60 && seconds < 60) {
        return(answer(sign * (degrees + minutes / 60 + seconds / 3600),
                      "packed degrees, minutes and seconds", compass = compass))
      }
    }
    return(answer(sign * value, paste0("decimal degrees", with_letter), compass = compass))
  }
  if (length(numbers) %in% 2:3) {
    degrees <- numbers[1]
    minutes <- numbers[2]
    seconds <- if (length(numbers) == 3) numbers[3] else 0
    if (degrees != round(degrees) || minutes >= 60 || seconds >= 60 ||
        (length(numbers) == 3 && minutes != round(minutes))) {
      return(answer(NA_real_, "not recognised", "degrees, minutes or seconds out of range", compass))
    }
    style <- if (length(numbers) == 2) "degrees and decimal minutes" else "degrees, minutes and seconds"
    return(answer(sign * (degrees + minutes / 60 + seconds / 3600), paste0(style, with_letter),
                  compass = compass))
  }
  answer(NA_real_, "not recognised", "number pattern not recognised", compass)
}



# 4. Reading the TRY release ----

# The release is read once, in pieces, and its Adansonia rows are saved next to
# it for later runs.

# Text that every kept line must contain; change it only if genus_pattern
# changes.
line_filter <- "dansonia"

try_path   <- if (given_path(try_data_file)) sub("(?<!:)[/\\\\]+$", "", try_data_file, perl = TRUE) else ""
try_folder <- if (try_path != "" && dir.exists(try_path)) try_path else if (try_path != "") dirname(try_path) else ""

# find_try_release_file() returns the .txt file whose first line holds the TRY
# column names.
find_try_release_file <- function(path) {
  if (path == "" || !file.exists(path)) return(NA_character_)
  if (!dir.exists(path)) return(path)
  candidates <- list.files(path, pattern = "\\.txt$", ignore.case = TRUE, recursive = TRUE, full.names = TRUE)
  is_release <- vapply(candidates, function(file) {
    header <- tryCatch(readLines(file, n = 1, warn = FALSE), error = function(e) "")
    all(vapply(c("AccSpeciesName", "ObservationID", "TraitID", "DataID"), grepl, logical(1),
               x = header, fixed = TRUE))
  }, logical(1))
  releases <- candidates[is_release]
  if (length(releases) == 0) return(NA_character_)
  if (length(releases) > 1) message("Several TRY release files were found; using the largest.")
  releases[which.max(file.size(releases))]
}

can_write <- function(folder) {
  test_file <- file.path(folder, ".bam_write_test")
  created <- tryCatch(file.create(test_file), warning = function(w) FALSE, error = function(e) FALSE)
  if (isTRUE(created)) unlink(test_file)
  isTRUE(created)
}

rows_name <- function(release_name) paste0(tools::file_path_sans_ext(release_name), "_adansonia_rows.tsv")
info_name <- function(release_name) paste0(tools::file_path_sans_ext(release_name), "_adansonia_rows_info.csv")

read_info <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(utils::read.csv(path, colClasses = "character"), error = function(e) NULL)
}

release_file <- find_try_release_file(try_path)

if (!is.na(release_file)) {
  message("TRY release file: ", release_file, " (", round(file.size(release_file) / 1e9, 1), " GB)")
  release_name <- basename(release_file)
  release_size <- format(file.size(release_file), scientific = FALSE)
  
  # The saved rows go next to the release, or to the results folder if that
  # folder is read-only.
  rows_folder <- if (can_write(dirname(release_file))) dirname(release_file) else try_results_folder
  subset_file <- file.path(rows_folder, rows_name(release_name))
  info_file   <- file.path(rows_folder, info_name(release_name))
  
  # Saved rows are reused if they come from a release of the same name and
  # size.
  saved_info <- read_info(info_file)
  rows_are_current <- file.exists(subset_file) && !is.null(saved_info) &&
    identical(saved_info$release[1], release_name) &&
    identical(saved_info$size_bytes[1], release_size) &&
    identical(saved_info$line_filter[1], line_filter)
  
  if (rows_are_current) {
    message("Using the rows saved by an earlier run (", basename(subset_file), "); the release is not read again")
    release_records <- as.numeric(saved_info$records[1])
  } else {
    message("Reading the TRY release once, in pieces, and keeping the rows with '", line_filter,
            "' (this takes several minutes; later runs use the saved rows)")
    # Rows go to a temporary file first, so that an interrupted run never
    # leaves a complete-looking file.
    partial_file <- paste0(subset_file, ".part")
    input  <- file(release_file, open = "rb")
    output <- file(partial_file, open = "wb")
    header <- readLines(input, n = 1, warn = FALSE)
    writeLines(header, output, useBytes = TRUE)
    release_records <- 0
    kept_rows <- 0
    repeat {
      lines <- readLines(input, n = 200000L, warn = FALSE)
      if (length(lines) == 0) break
      release_records <- release_records + length(lines)
      keep <- grepl(line_filter, lines, fixed = TRUE, useBytes = TRUE)
      if (any(keep)) writeLines(lines[keep], output, useBytes = TRUE)
      kept_rows <- kept_rows + sum(keep)
      if (release_records %% 5e6 < 200000) {
        message("  ", format(release_records, big.mark = ","), " lines read; ",
                format(kept_rows, big.mark = ","), " kept")
      }
    }
    close(input); close(output)
    if (file.exists(subset_file)) file.remove(subset_file)
    file.rename(partial_file, subset_file)
    utils::write.csv(data.frame(release = release_name, size_bytes = release_size,
                                records = format(release_records, scientific = FALSE),
                                line_filter = line_filter),
                     info_file, row.names = FALSE)
    message("  finished: ", format(release_records, big.mark = ","), " lines read; ",
            format(kept_rows, big.mark = ","), " rows saved in ", subset_file)
  }
} else {
  # Without a release, the rows saved by an earlier run are used.
  search_folders <- unique(Filter(function(folder) folder != "" && dir.exists(folder),
                                  c(try_folder, try_results_folder)))
  saved_rows <- unlist(lapply(search_folders, list.files, pattern = "_adansonia_rows\\.tsv$", full.names = TRUE))
  if (length(saved_rows) == 0) {
    stop("No TRY release was found at ", try_data_file, ", and no rows saved by an earlier run were found. ",
         "Set try_data_file in 1_Harmonise_tree_data.R (Section 1) to the TRY release (for example ",
         ".../52124.txt) or to the folder that holds it.", call. = FALSE)
  }
  subset_file <- saved_rows[which.max(file.mtime(saved_rows))]
  saved_info  <- read_info(sub("\\.tsv$", "_info.csv", subset_file))
  release_records <- if (!is.null(saved_info)) as.numeric(saved_info$records[1]) else NA_real_
  message("The TRY release was not found; using the rows saved by an earlier run: ", subset_file)
}

# Every column is read as text so that no value is changed.
try_data <- data.table::fread(subset_file, sep = "\t", quote = "", encoding = "Latin-1",
                              colClasses = "character", na.strings = c("", "NA"), fill = TRUE,
                              showProgress = FALSE)
try_data <- tibble::as_tibble(try_data)

# TRY files end with an empty column.
try_data <- try_data[, colSums(!is.na(try_data)) > 0 | !grepl("^V\\d+$", names(try_data))]

needed_columns <- c("LastName", "FirstName", "DatasetID", "Dataset", "SpeciesName", "AccSpeciesName",
                    "ObservationID", "TraitID", "TraitName", "DataID", "DataName", "OriglName",
                    "OrigValueStr", "OrigUnitStr", "StdValue", "UnitName", "ErrorRisk", "Reference")
missing_columns <- setdiff(needed_columns, names(try_data))
if (length(missing_columns) > 0) {
  stop("These TRY columns are missing: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}

try_data <- try_data %>%
  mutate(across(c(TraitID, DataID), as.integer),
         across(c(StdValue, ErrorRisk), ~ suppressWarnings(as.numeric(.x))))

try_overview <- tibble(records = release_records)
message("TRY release: ", format(release_records, big.mark = ","), " records; Adansonia rows read: ",
        format(nrow(try_data), big.mark = ","))


# 5. One row per TRY Adansonia observation ----

adansonia <- filter(try_data, str_detect(AccSpeciesName, regex(genus_pattern, ignore_case = TRUE)))
if (nrow(adansonia) == 0) stop("The TRY release holds no records of the genus Adansonia.", call. = FALSE)

# An observation may combine records from several datasets; all are listed.
observations <- adansonia %>%
  group_by(ObservationID) %>%
  summarise(species = collapse_unique(AccSpeciesName),
            species_as_recorded = collapse_unique(SpeciesName),
            dataset_id = collapse_unique(DatasetID),
            dataset = collapse_unique(Dataset),
            contributor = collapse_unique(str_squish(paste(FirstName, LastName))),
            contributor_surname = collapse_unique(LastName),
            reference = collapse_unique(Reference),
            .groups = "drop")

# coordinate_value() uses the standardised value, or converts the original text
# where TRY gives none.
coordinate_value <- function(data_id, which_one) {
  adansonia %>%
    filter(is.na(TraitID), DataID == data_id, !tolower(OriglName) %in% map_grid_names) %>%
    mutate(value = ifelse(!is.na(StdValue), StdValue,
                          vapply(OrigValueStr, function(text) {
                            if (is.na(text)) NA_real_ else read_single_coordinate(text, which_one)$value
                          }, numeric(1)))) %>%
    filter(!is.na(value)) %>%
    group_by(ObservationID) %>%
    summarise(value = first(value), .groups = "drop")
}
latitudes  <- rename(coordinate_value(latitude_data_id, "latitude"), latitude = value)
longitudes <- rename(coordinate_value(longitude_data_id, "longitude"), longitude = value)
altitudes  <- adansonia %>%
  filter(is.na(TraitID), DataID == altitude_data_id, !is.na(StdValue)) %>%
  group_by(ObservationID) %>% summarise(altitude_m = first(StdValue), .groups = "drop")
map_grid_only <- adansonia %>%
  filter(is.na(TraitID), tolower(OriglName) %in% map_grid_names) %>%
  distinct(ObservationID) %>% mutate(has_map_grid_coordinates = TRUE)

# Each numeric trait becomes one column, averaged where an observation has
# several values.
trait_values <- adansonia %>%
  filter(!is.na(TraitID), !is.na(StdValue)) %>%
  mutate(column = paste0("try_", TraitID, "_", to_dictionary_style(TraitName),
                         ifelse(is.na(UnitName), "", paste0("_", to_dictionary_style(UnitName))))) %>%
  group_by(ObservationID, column) %>%
  summarise(value = mean(StdValue), .groups = "drop") %>%
  pivot_wider(names_from = column, values_from = value)
heights <- adansonia %>%
  filter(TraitID %in% height_trait_id, !is.na(StdValue)) %>%
  group_by(ObservationID) %>%
  summarise(try_height_m = mean(StdValue),
            try_height_error_risk = if (all(is.na(ErrorRisk))) NA_real_ else max(ErrorRisk, na.rm = TRUE),
            .groups = "drop")

try_observations <- observations %>%
  left_join(latitudes, by = "ObservationID") %>%
  left_join(longitudes, by = "ObservationID") %>%
  left_join(altitudes, by = "ObservationID") %>%
  left_join(map_grid_only, by = "ObservationID") %>%
  left_join(heights, by = "ObservationID") %>%
  left_join(trait_values, by = "ObservationID") %>%
  mutate(has_map_grid_coordinates = coalesce(has_map_grid_coordinates, FALSE),
         valid_position = !is.na(latitude) & !is.na(longitude) & abs(latitude) <= 90 & abs(longitude) <= 180)
message("TRY Adansonia observations: ", nrow(try_observations), " (",
        sum(try_observations$valid_position), " with coordinates)")


# 6. BAM baobab records ----

if (!file.exists(public_data_file)) {
  stop("The public data file was not found: ", public_data_file,
       ". Render 2_QAQC_BAM_Data.qmd first.", call. = FALSE)
}
# The baobab records of the public file. The species as recorded is withheld
# from that file, so baobabs are identified by tree type.
public_baobabs <- read_csv(public_data_file, col_types = cols(.default = col_character()),
                           na = c("", "NA"), progress = FALSE) %>% rename_to_roles() %>%
  filter(tree_type %in% "baobab")
bam <- public_baobabs %>%
  mutate(latitude = as.numeric(latitude), longitude = as.numeric(longitude),
         height = as.numeric(height)) %>%
  filter(!is.na(latitude), !is.na(longitude))
message("BAM baobab records with coordinates in the public file: ", format(nrow(bam), big.mark = ","))

# Words that identify BAM contributors: the observer names published with
# consent. "Anonymous" identifies no one and is left out.
bam_name_words <- unique(tolower(unlist(str_split(na.omit(bam$observer_name), "[^\\p{L}]+"))))
bam_name_words <- setdiff(bam_name_words[nchar(bam_name_words) >= 3], "anonymous")


# 7. Comparing TRY with BAM ----

sf_use_s2(TRUE)
bam_points <- st_as_sf(bam, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
located <- filter(try_observations, valid_position)
try_points <- st_as_sf(located, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

comparison <- tibble(ObservationID = character())
if (nrow(try_points) > 0) {
  nearest <- st_nearest_feature(try_points, bam_points)
  distance_m <- as.numeric(st_distance(try_points, bam_points[nearest, ], by_element = TRUE))
  bam_within_site <- lengths(st_is_within_distance(try_points, bam_points, dist = same_site_kilometres * 1000))
  comparison <- tibble(
    ObservationID = located$ObservationID,
    nearest_bam_distance_m = round(distance_m, 1),
    nearest_bam_record_id = bam$record_id[nearest],
    nearest_bam_height_m = bam$height[nearest],
    bam_records_within_site = bam_within_site
  )
}

countries <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  select(country = name_long) %>% st_make_valid()
country_of <- if (nrow(try_points) > 0) {
  index <- vapply(suppressMessages(st_intersects(try_points, countries)),
                  function(i) if (length(i)) i[1] else NA_integer_, integer(1))
  tibble(ObservationID = located$ObservationID, country = countries$country[index])
} else tibble(ObservationID = character(), country = character())

duplicate_check <- try_observations %>%
  left_join(comparison, by = "ObservationID") %>%
  left_join(country_of, by = "ObservationID") %>%
  mutate(
    heights_agree = !is.na(try_height_m) & !is.na(nearest_bam_height_m) &
      abs(try_height_m - nearest_bam_height_m) <= pmax(1, height_tolerance * nearest_bam_height_m),
    contributor_in_bam = vapply(tolower(coalesce(contributor_surname, "")), function(names) {
      surnames <- str_split(names, ";\\s*")[[1]]
      any(surnames %in% bam_name_words & nchar(surnames) >= 3)
    }, logical(1)),
    same_tree = !is.na(nearest_bam_distance_m) & nearest_bam_distance_m <= same_tree_metres,
    same_site = !is.na(nearest_bam_distance_m) & nearest_bam_distance_m <= same_site_kilometres * 1000,
    assessment = case_when(
      !valid_position & has_map_grid_coordinates ~ "map-grid coordinates only: check by hand",
      !valid_position                              ~ "no coordinates: check by hand",
      same_tree & (heights_agree | contributor_in_bam) ~ "probable duplicate",
      same_tree                                    ~ "possible duplicate",
      same_site & contributor_in_bam               ~ "possible duplicate",
      same_site                                    ~ "same site as BAM records, different tree",
      TRUE                                         ~ "not in BAM"
    ),
    reason = case_when(
      assessment == "probable duplicate" ~ paste0("within ", same_tree_metres, " m of a BAM tree",
                                                  ifelse(heights_agree, "; heights agree", ""),
                                                  ifelse(contributor_in_bam, "; contributor also in BAM", "")),
      assessment == "possible duplicate" & same_tree ~ paste0("within ", same_tree_metres, " m of a BAM tree"),
      assessment == "possible duplicate" ~ paste0("within ", same_site_kilometres, " km of BAM records; contributor also in BAM"),
      TRUE ~ NA_character_)
  ) %>%
  select(ObservationID, assessment, reason, species, country, latitude, longitude, altitude_m,
         try_height_m, nearest_bam_distance_m, nearest_bam_record_id,
         nearest_bam_height_m, bam_records_within_site, heights_agree,
         contributor_in_bam, contributor, dataset_id, dataset, reference) %>%
  arrange(factor(assessment, levels = c("probable duplicate", "possible duplicate",
                                        "same site as BAM records, different tree", "not in BAM",
                                        "map-grid coordinates only: check by hand",
                                        "no coordinates: check by hand")), ObservationID)


# 8. Trait coverage for A. digitata ----

target <- filter(adansonia, AccSpeciesName == target_species, !is.na(TraitID))
located_ids <- try_observations$ObservationID[try_observations$valid_position]
try_traits <- target %>%
  group_by(TraitID, TraitName) %>%
  summarise(try_records = n(), try_numeric_values = sum(!is.na(StdValue)),
            try_observations = n_distinct(ObservationID),
            try_observations_with_coordinates = n_distinct(ObservationID[ObservationID %in% located_ids]),
            try_values_high_error_risk = sum(!is.na(ErrorRisk) & ErrorRisk > error_risk_limit),
            try_units = collapse_unique(UnitName), .groups = "drop")

# All baobab records of the public file, with or without coordinates.
all_bam <- public_baobabs
trait_coverage <- try_traits %>%
  mutate(bam_variable = vapply(tolower(TraitName), function(name) {
    hit <- trait_matching$bam_variable[str_detect(name, trait_matching$try_trait_words)]
    if (length(hit) == 0) NA_character_ else hit[1]
  }, character(1))) %>%
  full_join(select(trait_matching, bam_variable), by = "bam_variable") %>%
  mutate(bam_baobab_records_with_value = vapply(bam_variable, function(variable) {
    if (is.na(variable) || !variable %in% names(all_bam)) NA_integer_ else sum(!is.na(all_bam[[variable]]))
  }, integer(1))) %>%
  arrange(is.na(TraitID), desc(try_observations))


# 9. Saving the results 

summary_table <- tibble(
  measure = c("TRY records in the release", "TRY Adansonia observations",
              paste("TRY", target_species, "observations"),
              "TRY Adansonia observations with coordinates",
              "Probable duplicates of BAM records", "Possible duplicates of BAM records",
              "Same site as BAM records, different tree", "Not in BAM",
              "Without usable coordinates (check by hand)", "BAM baobab records with coordinates"),
  value = c(try_overview$records, nrow(try_observations),
            sum(str_detect(try_observations$species, fixed(target_species))),
            sum(try_observations$valid_position),
            sum(duplicate_check$assessment == "probable duplicate"),
            sum(duplicate_check$assessment == "possible duplicate"),
            sum(duplicate_check$assessment == "same site as BAM records, different tree"),
            sum(duplicate_check$assessment == "not in BAM"),
            sum(str_detect(duplicate_check$assessment, "check by hand")),
            nrow(bam))
)

# Figures for A. digitata alone, as quoted in the Technical Validation.
target_check <- filter(duplicate_check, str_detect(species, fixed(target_species)))
target_trait_names <- sort(unique(na.omit(target$TraitName)))
summary_table <- bind_rows(
  mutate(summary_table, value = as.character(value)),
  tibble(
    measure = paste(target_species, c("observations in TRY", "observations with coordinates",
                                      "traits in TRY", "trait names in TRY",
                                      paste("observations within", same_tree_metres, "m of a BAM baobab"),
                                      "probable duplicates of BAM records")),
    value = c(nrow(target_check),
              sum(!str_detect(target_check$assessment, "check by hand")),
              length(target_trait_names),
              paste(target_trait_names, collapse = "; "),
              sum(target_check$nearest_bam_distance_m <= same_tree_metres, na.rm = TRUE),
              sum(target_check$assessment == "probable duplicate")) %>% as.character()
  )
)

# The MD5 checksum of the public file read here lets 5_Summary_report_BAM.qmd
# confirm that these results match the current public file.

summary_table <- bind_rows(summary_table, tibble(
  measure = c("TRY release version", "Public data file MD5 checksum"),
  value = c(if (given_path(try_release_version)) try_release_version else "not recorded",
            unname(tools::md5sum(public_data_file)))))

write_excel_csv(try_observations, try_observations_output, na = "NA")
write_excel_csv(duplicate_check, duplicate_check_output, na = "NA")
write_excel_csv(trait_coverage, trait_coverage_output, na = "NA")
write_excel_csv(summary_table, summary_output, na = "NA")


# 10. Map of BAM and TRY locations (for checking; not in the manuscript)

text_size_pt <- 7
small_text_size_pt <- 6
line_width <- function(points) points / (ggplot2::.pt * 72 / 96)  # points to ggplot2 linewidth
figure_font <- intersect(c("Arial", "Helvetica", "Liberation Sans", "Arimo"),
                         unique(systemfonts::system_fonts()$family))[1]
if (is.na(figure_font)) figure_font <- "sans"
format_longitude <- function(x) paste0(abs(x), "\u00b0", ifelse(x > 0, "E", ifelse(x < 0, "W", "")))
format_latitude  <- function(x) paste0(abs(x), "\u00b0", ifelse(x > 0, "N", ifelse(x < 0, "S", "")))

sf_use_s2(FALSE)
world_countries <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
africa <- filter(world_countries, continent == "Africa")

# The map shows mainland Africa and Madagascar; small outlying islands do not
# set the frame.
africa_parts <- st_cast(st_geometry(africa), "POLYGON")
africa_frame <- st_bbox(africa_parts[as.numeric(st_area(st_set_crs(africa_parts, NA))) >= 1])
map_extent <- c(xmin = africa_frame[["xmin"]] - 2, xmax = africa_frame[["xmax"]] + 2,
                ymin = africa_frame[["ymin"]] - 2, ymax = africa_frame[["ymax"]] + 2)
in_frame <- function(x, y) x >= map_extent[["xmin"]] & x <= map_extent[["xmax"]] &
  y >= map_extent[["ymin"]] & y <= map_extent[["ymax"]]

# colours remain distinct under common colour-vision deficiencies.

assessment_colours <- c("probable duplicate" = "#D55E00", "possible duplicate" = "#E69F00",
                        "same site as BAM records, different tree" = "#009E73", "not in BAM" = "#0072B2")
map_points <- filter(duplicate_check, !is.na(latitude), !is.na(longitude),
                     assessment %in% names(assessment_colours), in_frame(longitude, latitude)) %>%
  mutate(assessment = factor(assessment, levels = names(assessment_colours))) %>%
  arrange(desc(assessment))  # duplicates are drawn last
outside_frame <- sum(try_observations$valid_position, na.rm = TRUE) - nrow(map_points)
if (outside_frame > 0) message(outside_frame, " TRY observations lie outside Africa and are not mapped.")

bam_in_frame <- filter(bam, in_frame(longitude, latitude))
legend_labels <- paste0(c("Probable duplicate", "Possible duplicate", "Same site as BAM records, different tree",
                          "Not in BAM"), " (n = ",
                        table(map_points$assessment)[names(assessment_colours)], ")")

supplementary_figure_1 <- ggplot() +
  geom_sf(data = world_countries, fill = "white", colour = "grey75", linewidth = line_width(0.25)) +
  geom_sf(data = africa, fill = "grey92", colour = "grey55", linewidth = line_width(0.25)) +
  geom_point(data = bam_in_frame, aes(x = longitude, y = latitude, colour = "BAM baobab records"),
             size = 0.4, alpha = 0.6, shape = 16) +
  geom_point(data = map_points, aes(x = longitude, y = latitude, fill = assessment),
             shape = 21, colour = "white", stroke = line_width(0.25), size = 1.8) +
  scale_colour_manual(values = c("BAM baobab records" = "red"), name = NULL,
                      labels = paste0("BAM baobab records (n = ", format(nrow(bam_in_frame), big.mark = ","), ")")) +
  scale_fill_manual(values = assessment_colours, labels = legend_labels, drop = FALSE,
                    name = "TRY Adansonia observations") +
  scale_x_continuous(breaks = seq(-20, 60, by = 20), labels = format_longitude) +
  scale_y_continuous(breaks = seq(-30, 30, by = 15), labels = format_latitude) +
  coord_sf(xlim = unname(map_extent[c("xmin", "xmax")]), ylim = unname(map_extent[c("ymin", "ymax")]),
           expand = FALSE) +
  guides(colour = guide_legend(order = 1, override.aes = list(size = 1.5, alpha = 1)),
         fill = guide_legend(order = 2, ncol = 1, title.position = "top", override.aes = list(size = 2))) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = text_size_pt, base_family = figure_font) +
  theme(text = element_text(colour = "black"),
        axis.text = element_text(size = small_text_size_pt, colour = "black"),
        axis.ticks = element_line(linewidth = line_width(0.5), colour = "black"),
        axis.line = element_blank(),
        panel.background = element_rect(fill = "#eef3f8", colour = NA),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = line_width(0.5)),
        legend.position = "bottom", legend.box = "vertical", legend.box.just = "left",
        legend.title = element_text(size = text_size_pt), legend.text = element_text(size = small_text_size_pt),
        legend.key.size = unit(3, "mm"), legend.margin = margin(0, 0, 0, 0),
        legend.spacing.y = unit(1, "mm"), plot.margin = margin(2, 3, 2, 2, unit = "mm"))

figure_file <- figure_output_stem
ggsave(paste0(figure_file, ".pdf"), supplementary_figure_1, width = 120, height = 150, units = "mm",
       bg = "white", device = cairo_pdf)
ggsave(paste0(figure_file, ".png"), supplementary_figure_1, width = 120, height = 150, units = "mm",
       dpi = 600, bg = "white", device = ragg::agg_png)


# 11. Report 

message("\nCross-check of BAM with TRY:")
for (i in seq_len(nrow(summary_table))) {
  message(sprintf("  %-50s %s", summary_table$measure[i], format(summary_table$value[i], big.mark = ",")))
}
message("\nResults and the check map saved in ", try_results_folder,
        "\nRead try_bam_duplicate_check.csv before deciding whether any record is a duplicate.",
        if (!given_path(try_release_version)) "\nSet try_release_version (Section 2): the TRY version is cited in the manuscript." else "")