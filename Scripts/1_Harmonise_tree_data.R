# ------------------------------------------------------------------------------ 
# Introduction
# ------------------------------------------------------------------------------

# This script reads every contributor file (CSV, TSV, Excel, GPX, KML) in the 
# raw input folder, matches each column to BAM_variable_dictionary.csv and 
# writes one record per tree observation. Coordinates are converted to WGS84 
# decimal degrees, sizes to dictionary units and observer names to one spelling 
# per person; every change is logged.
#
# The pipeline runs in five numbered steps, one script each, in this order:
#   1_Harmonise_tree_data.R       harmonisation (this script)
#   2_QAQC_BAM_Data.qmd           quality control and the public data file
#   3_TRY_comparison_check.R      comparison with the TRY Plant Trait Database
#   4_Environmental_covariates.R  environmental covariates and Figure 3
#   5_Summary_report_BAM.qmd      dataset summaries and Figure 2
# Every path and setting of the pipeline is set in Section 1 of this script.


# Paths and setup


input_tables_folder <- "data/input_tables"  # correction lists (not public)
raw_input_folder    <- "data/raw_input"     # one sub-folder per data collection method
results_folder      <- "data/outputs"       # results of every step

# The variable dictionary, read from the input tables folder 

dictionary_name <- "BAM_variable_dictionary.csv"

# World Flora Online backbone (classification.csv or its folder); "" skips

wfo_backbone_file <- "data/WFO"

# Tree photographs, searched with all sub-folders

photo_source_folder <- "data/Bam Photos"

# TRY release (the .txt file or its folder) 

try_data_file <- "data/TRY_Dataset"

# Environmental data for step 4: download folder, climate period and optional
# local copies ("" downloads the data).

environment_data_folder           <- "data/environment"
climate_period                    <- "1991:2020"
chirps_local_folder               <- ""
terraclimate_pet_local_source     <- ""
era5land_temperature_local_source <- ""
isda_local_folder                 <- ""
protected_areas_local_source      <- ""
dem_local_source                  <- ""  # elevation model (GeoTIFF); "" uses Copernicus DEM GLO-30

# IUCN Red List range polygons of the African elephants, used in step 2 to
# check reports of elephant damage: the folder holding the GeoPackages (or
# shapefiles) downloaded from https://www.iucnredlist.org for the savanna
# elephant (Loxodonta africana) and the forest elephant (L. cyclotis)

elephant_range_folder <- "data/elephant_geopackage"
iucn_red_list_version <- "3"

# Files in which a tree recorded on several sheets carries the same
# identifier on every sheet, so that its position may be copied from one
# sheet to another

files_with_trees_shared_across_sheets <- c()

# UTM zone for files that do not state one, for example "36S"; NA infers it.

utm_zone_when_not_stated <- NA

# Species given to records from files that never name one; NA leaves it empty.

species_when_not_recorded <- "Adansonia digitata"

# TRUE reads ambiguous dates such as 03/04/2021 as 3 April.

day_first_dates <- TRUE

# Files whose observation dates are taken from the photograph file names.

files_dated_from_photographs <- c("Kremetart_negative_trees.csv")

# TRUE also removes trees recorded as dead in a status column

remove_trees_recorded_as_dead <- TRUE


# packages and saved paths

# rstudioapi is optional: it finds the folder of this script when run in
# RStudio.

if (requireNamespace("rstudioapi", quietly = TRUE)) library(rstudioapi)

script_folders <- c(
  tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) NA_character_),
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    tryCatch(dirname(rstudioapi::getSourceEditorContext()$path), error = function(e) NA_character_)
  } else NA_character_,
  getwd())
searched_folders <- unique(script_folders[!is.na(script_folders) & script_folders != ""])
script_folders <- searched_folders[file.exists(file.path(searched_folders, "1_Harmonise_tree_data.R"))]
if (length(script_folders) == 0) {
  stop("The scripts folder was not found (no 1_Harmonise_tree_data.R in: ",
       paste(searched_folders, collapse = "; "), "). Open this script in RStudio from the scripts ",
       "folder and click Source, or first run setwd() with the path of the scripts folder.", call. = FALSE)
}
setwd(script_folders[1])
scripts_folder <- normalizePath(getwd(), winslash = "/")
project_folder <- dirname(scripts_folder)

required_packages <- c("readr", "readxl", "dplyr", "purrr", "stringr", "stringdist",
                       "lubridate", "tibble", "sf", "rnaturalearth", "rnaturalearthdata")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Please install these packages first: install.packages(c(\"",
       paste(missing_packages, collapse = "\", \""), "\"))", call. = FALSE)
}
suppressPackageStartupMessages({
  library(readr)              # reading and writing CSV files
  library(readxl)             # reading Excel workbooks
  library(dplyr)              # table manipulation
  library(purrr)              # applying functions to lists
  library(stringr)            # text handling
  library(stringdist)         # similarity of observer names
  library(lubridate)          # reading dates
  library(tibble)             # tables
  library(sf)                 # spatial checks of coordinates
  library(rnaturalearth)      # country boundaries
  library(rnaturalearthdata)  # data used by rnaturalearth
})

in_project_folder <- function(path) {
  if (is.na(path) || path == "") return(path)
  path <- gsub("\\\\", "/", path)
  if (grepl("^([A-Za-z]:)?/", path) || startsWith(path, "~")) path else file.path(project_folder, path)
}
for (path_name in c("input_tables_folder", "raw_input_folder", "results_folder", "wfo_backbone_file",
                    "photo_source_folder",
                    "try_data_file", "environment_data_folder", "chirps_local_folder",
                    "terraclimate_pet_local_source", "era5land_temperature_local_source",
                    "isda_local_folder", "protected_areas_local_source", "dem_local_source",
                    "elephant_range_folder")) {
  assign(path_name, in_project_folder(get(path_name)))
}

# One results folder per step, numbered as the scripts; figures of steps 4
# and 5 share one folder.

harmonised_results_folder  <- file.path(results_folder, "1_HARMONISED")
qaqc_results_folder        <- file.path(results_folder, "2_QAQC")
try_results_folder         <- file.path(results_folder, "3_TRY_CROSS_CHECK")
environment_results_folder <- file.path(results_folder, "4_ENVIRONMENT")
figures_folder             <- file.path(results_folder, "5_FIGURES")

# The dictionary is read from the input tables folder, or from the scripts
# folder if it is kept there; the first found is used.
dictionary_places <- c(file.path(input_tables_folder, dictionary_name), file.path(scripts_folder, dictionary_name))
dictionary_file <- dictionary_places[file.exists(dictionary_places)][1]
if (is.na(dictionary_file)) {
  stop("The dictionary ", dictionary_name, " was not found in ", input_tables_folder, " or ", scripts_folder, ".",
       call. = FALSE)
}
message("Dictionary used: ", dictionary_file)
species_corrections_file <- file.path(input_tables_folder, "species_name_corrections.csv")
name_corrections_file    <- file.path(input_tables_folder, "name_corrections.csv")
harmonised_dataset_file  <- file.path(harmonised_results_folder, "baobab_tree_dataset.csv")
public_data_file         <- file.path(qaqc_results_folder, "public_data", "BAM_data_for_public.csv")

if (!file.exists(dictionary_file)) stop("The dictionary was not found: ", dictionary_file, call. = FALSE)
if (!dir.exists(raw_input_folder)) stop("The raw input folder was not found: ", raw_input_folder, call. = FALSE)
dir.create(harmonised_results_folder, recursive = TRUE, showWarnings = FALSE)

# Every path and setting is saved in pipeline_paths.csv, which the other
# scripts read: 2_QAQC_BAM_Data.qmd, 3_TRY_comparison_check.R,
# 4_Environmental_covariates.R and 5_Summary_report_BAM.qmd. Paths are set
# here only.

path_names <- c("raw_input_folder", "results_folder", "harmonised_results_folder", "try_results_folder",
                "qaqc_results_folder", "figures_folder", "environment_results_folder", "dictionary_file",
                "input_tables_folder", "species_corrections_file", "name_corrections_file",
                "harmonised_dataset_file", "public_data_file", "wfo_backbone_file", "photo_source_folder",
                "try_data_file", "environment_data_folder", "climate_period", "chirps_local_folder",
                "terraclimate_pet_local_source", "era5land_temperature_local_source", "isda_local_folder",
                "protected_areas_local_source", "dem_local_source", "elephant_range_folder",
                "iucn_red_list_version")
pipeline_paths <- tibble(name = path_names,
                         value = vapply(path_names, function(n) as.character(get(n)), character(1)))
paths_file <- file.path(scripts_folder, "pipeline_paths.csv")
write_csv(pipeline_paths, paths_file)
message("Paths saved in ", paths_file, ":\n",
        paste0("  ", format(pipeline_paths$name), "  ", pipeline_paths$value, collapse = "\n"))
start_time <- Sys.time()


# create some helper function

# tidy_text() writes text in a comparable form: lower case, no accents and
# single spaces.

tidy_text <- function(x) {
  x <- iconv(as.character(x), from = "UTF-8", to = "ASCII//TRANSLIT", sub = " ")
  x <- str_replace_all(str_to_lower(x), "&", " and ")
  str_squish(str_replace_all(x, "[^a-z0-9]+", " "))
}

is_blank <- function(x) is.na(x) | str_trim(x) == ""

as_number <- function(x) suppressWarnings(as.numeric(str_replace(str_trim(x), ",", ".")))

# Answers that mean yes, no or unknown on the field forms.

yes_words     <- c("yes", "y", "true", "t", "1", "present", "ja", "oui", "si")
no_words      <- c("no", "n", "false", "f", "0", "absent", "none", "nee", "non", "nil")
unknown_words <- c("unknown", "unsure", "not sure", "dont know", "don t know", "na", "n a", "not known")

# Values that cannot be a person's name.

not_a_name_words <- c("dead", "baobab", "tree", "unknown", "none", "test", "testing", "na", "n a",
                      "yes", "no", "missing", "plot", "site", "email")

# Distance in metres between two points (haversine formula).

metres_apart <- function(lat1, lon1, lat2, lon2) {
  to_rad <- pi / 180
  a <- sin((lat2 - lat1) * to_rad / 2)^2 +
    cos(lat1 * to_rad) * cos(lat2 * to_rad) * sin((lon2 - lon1) * to_rad / 2)^2
  a <- pmin(pmax(a, 0), 1)
  6371008.8 * 2 * atan2(sqrt(a), sqrt(1 - a))
}

# Logs for the files checked after a run.

column_log     <- list()
coordinate_log <- list()
value_log      <- list()
sheet_log      <- list()

add_log <- function(log_name, rows) {
  current <- get(log_name, envir = .GlobalEnv)
  current[[length(current) + 1]] <- rows
  assign(log_name, current, envir = .GlobalEnv)
}


# Coordinate conversion 

mgrs_alphabet <- strsplit("ABCDEFGHJKLMNPQRSTUVWXYZ", "")[[1]]   # I and O are not used.
mgrs_band_letters <- strsplit("CDEFGHJKLMNPQRSTUVWX", "")[[1]]   # 8-degree bands from 80 S.
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

# utm_to_degrees() applies the inverse transverse Mercator series on the WGS84
# ellipsoid.

utm_to_degrees <- function(easting, northing, zone, southern) {
  a   <- 6378137
  f   <- 1 / 298.257223563
  k0  <- 0.9996
  e2  <- f * (2 - f)
  ep2 <- e2 / (1 - e2)
  x <- easting - 500000
  y <- northing - ifelse(southern, 10000000, 0)
  mu <- (y / k0) / (a * (1 - e2 / 4 - 3 * e2^2 / 64 - 5 * e2^3 / 256))
  e1 <- (1 - sqrt(1 - e2)) / (1 + sqrt(1 - e2))
  phi1 <- mu + (3 * e1 / 2 - 27 * e1^3 / 32) * sin(2 * mu) +
    (21 * e1^2 / 16 - 55 * e1^4 / 32) * sin(4 * mu) +
    (151 * e1^3 / 96) * sin(6 * mu) + (1097 * e1^4 / 512) * sin(8 * mu)
  n1 <- a / sqrt(1 - e2 * sin(phi1)^2)
  t1 <- tan(phi1)^2
  c1 <- ep2 * cos(phi1)^2
  r1 <- a * (1 - e2) / (1 - e2 * sin(phi1)^2)^1.5
  d  <- x / (n1 * k0)
  latitude <- phi1 - (n1 * tan(phi1) / r1) *
    (d^2 / 2 - (5 + 3 * t1 + 10 * c1 - 4 * c1^2 - 9 * ep2) * d^4 / 24 +
       (61 + 90 * t1 + 298 * c1 + 45 * t1^2 - 252 * ep2 - 3 * c1^2) * d^6 / 720)
  longitude <- (d - (1 + 2 * t1 + c1) * d^3 / 6 +
                  (5 - 2 * c1 + 28 * t1 - 3 * c1^2 + 8 * ep2 + 24 * t1^2) * d^5 / 120) / cos(phi1)
  central_meridian <- (zone - 1) * 6 - 177
  list(latitude = latitude * 180 / pi, longitude = central_meridian + longitude * 180 / pi)
}

# mgrs_to_degrees() returns NULL if the text is not a valid MGRS reference.

mgrs_to_degrees <- function(text) {
  compact_text <- gsub("\\s", "", toupper(text), perl = TRUE)
  parts <- regmatches(compact_text, regexec(
    "^(\\d{1,2})([C-HJ-NP-X])([A-HJ-NP-Z])([A-HJ-NP-V])(\\d*)$", compact_text, perl = TRUE))[[1]]
  if (length(parts) == 0) return(NULL)
  zone <- as.integer(parts[2]); band <- parts[3]
  column_letter <- parts[4]; row_letter <- parts[5]; digits <- parts[6]
  if (zone < 1 || zone > 60 || nchar(digits) %% 2 != 0 || nchar(digits) > 10) return(NULL)
  # Column letters cycle through three sets of eight; row letters every
  # 2,000 km, starting five letters later in even zones.
  column_index <- match(column_letter, mgrs_alphabet) - ((zone - 1) %% 3) * 8
  if (column_index < 1 || column_index > 8) return(NULL)
  row_index <- match(row_letter, mgrs_alphabet[1:20]) - 1
  if (is.na(row_index)) return(NULL)
  if (zone %% 2 == 0) row_index <- (row_index - 5) %% 20
  precision  <- nchar(digits) / 2
  scale      <- 10^(5 - precision)
  easting_m  <- if (precision > 0) as.numeric(substr(digits, 1, precision)) * scale else 0
  northing_m <- if (precision > 0) as.numeric(substr(digits, precision + 1, 2 * precision)) * scale else 0
  easting    <- column_index * 100000 + easting_m
  northing   <- row_index * 100000 + northing_m
  # The latitude band fixes the 2,000 km cycle of the northing.
  band_index <- match(band, mgrs_band_letters)
  lowest     <- -80 + 8 * (band_index - 1)
  highest    <- if (band == "X") 84 else lowest + 8
  southern   <- band_index <= match("M", mgrs_band_letters)
  for (cycle in 0:5) {
    candidate <- northing + cycle * 2000000
    if (candidate > 10000000) break
    position <- utm_to_degrees(easting, candidate, zone, southern)
    if (position$latitude >= lowest - 0.5 && position$latitude <= highest + 0.5) {
      return(list(latitude = position$latitude, longitude = position$longitude,
                  style = paste0("MGRS grid reference (precision ", scale, " m)")))
    }
  }
  NULL
}

read_grid_number <- function(text) {
  text <- toupper(trimws(as.character(text)))
  text <- gsub("[EN]$|^[EN]", "", text, perl = TRUE)
  text <- gsub("(?<=\\d)[ ,](?=\\d{3}(\\D|$))", "", text, perl = TRUE)
  suppressWarnings(as.numeric(trimws(gsub(",", ".", text, perl = TRUE))))
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

# read_whole_position() reads a cell holding both coordinates, or returns NULL.

read_whole_position <- function(text) {
  result <- function(latitude, longitude, style, note = NA_character_) {
    list(latitude = latitude, longitude = longitude, style = style, note = note)
  }
  raw <- toupper(trimws(as.character(text)))
  if (is.na(raw) || raw == "") return(NULL)
  
  # WKT point: POINT (longitude latitude).
  
  wkt <- regmatches(raw, regexec(
    "^POINT\\s*Z?M?\\s*\\(\\s*([-+]?[0-9.]+(?:E[-+]?\\d+)?)\\s+([-+]?[0-9.]+(?:E[-+]?\\d+)?)", raw, perl = TRUE))[[1]]
  if (length(wkt) == 3) return(result(as.numeric(wkt[3]), as.numeric(wkt[2]), "WKT point (longitude first)"))
  
  # MGRS grid reference.
  
  if (grepl("^\\d{1,2}\\s*[C-HJ-NP-X]\\s*[A-HJ-NP-Z][A-HJ-NP-V]\\s*\\d*\\s*\\d*$", raw, perl = TRUE)) {
    position <- mgrs_to_degrees(raw)
    if (!is.null(position)) return(result(position$latitude, position$longitude, position$style))
    return(result(NA_real_, NA_real_, "not recognised", "looks like an MGRS reference but is not valid"))
  }
  
  # UTM with the zone written in the cell.
  
  utm_text <- trimws(gsub("\\s+", " ", gsub("\\b(UTM|ZONE|WGS ?84)\\b|[,;:]", " ", raw, perl = TRUE)))
  utm <- regmatches(utm_text, regexec(
    "^(\\d{1,2})\\s*([C-HJ-NP-X]|NORTH|SOUTH)?\\s+([EN]?\\s*\\d[\\d.]*[EN]?)\\s+([EN]?\\s*\\d[\\d.]*[EN]?)$",
    utm_text, perl = TRUE))[[1]]
  if (length(utm) == 5) {
    zone <- as.integer(utm[2])
    first <- utm[4]; second <- utm[5]
    first_number <- read_grid_number(first); second_number <- read_grid_number(second)
    if (zone >= 1 && zone <= 60 && !is.na(first_number) && !is.na(second_number) &&
        first_number >= 1000 && second_number >= 1000) {
      letter <- utm[3]
      if (letter == "") {
        return(result(NA_real_, NA_real_, "not recognised",
                      "UTM position without a hemisphere or latitude band; add it (for example 36S)"))
      }
      southern <- letter %in% c("S", "SOUTH") ||
        (nchar(letter) == 1 && match(letter, mgrs_band_letters) <= match("M", mgrs_band_letters))
      # Letters E and N identify each value; otherwise the easting comes first.
      if (grepl("N", first, perl = TRUE) || grepl("E", second, perl = TRUE) ||
          (first_number >= 1000000 && second_number < 1000000)) {
        easting <- second_number; northing <- first_number
      } else {
        easting <- first_number; northing <- second_number
      }
      if (easting < 100000 || easting > 900000 || northing > 10000000) {
        return(result(NA_real_, NA_real_, "not recognised", "UTM easting or northing out of range"))
      }
      position <- utm_to_degrees(easting, northing, zone, southern)
      return(result(position$latitude, position$longitude,
                    paste0("UTM in the cell (zone ", zone, if (southern) " south" else " north", ", WGS 84)")))
    }
  }
  
  text <- standardise_coordinate_text(raw)
  
  # Both coordinates with compass letters.
  
  compass <- regmatches(text, gregexpr(compass_letter_pattern, text, perl = TRUE))[[1]]
  if (length(compass) == 2 && any(compass %in% c("N", "S")) && any(compass %in% c("E", "W"))) {
    parts <- regmatches(text, gregexpr("(?:[NSEW]\\s*-?[\\d.\\s]+)|(?:-?[\\d.\\s]+\\s*[NSEW])", text,
                                       perl = TRUE))[[1]]
    latitude_part  <- parts[grepl("[NS]", parts, perl = TRUE)]
    longitude_part <- parts[grepl("[EW]", parts, perl = TRUE)]
    if (length(latitude_part) == 1 && length(longitude_part) == 1) {
      latitude  <- read_single_coordinate(latitude_part, "latitude")
      longitude <- read_single_coordinate(longitude_part, "longitude")
      note <- c(latitude$note, longitude$note)
      note <- if (all(is.na(note))) NA_character_ else paste(na.omit(note), collapse = "; ")
      return(result(latitude$value, longitude$value,
                    paste0("both coordinates in one cell (", latitude$style, ")"), note))
    }
    return(result(NA_real_, NA_real_, "not recognised", "both coordinates in one cell, could not be separated"))
  }
  
  # A decimal pair or a KoboToolbox geopoint; both numbers need decimals.
  pair <- regmatches(text, regexec(
    "^(-?\\d+\\.\\d+)\\s+(-?\\d+\\.\\d+)((\\s+-?\\d+(\\.\\d+)?){0,2})$", text, perl = TRUE))[[1]]
  if (length(pair) >= 3) {
    first <- as.numeric(pair[2]); second <- as.numeric(pair[3])
    style <- if (trimws(pair[4]) != "") "GPS geopoint (latitude, longitude, altitude, accuracy)" else
      "latitude-longitude pair in one cell"
    if (abs(first) <= 90 && abs(second) <= 180) return(result(first, second, style))
    if (abs(second) <= 90 && abs(first) <= 180) {
      return(result(second, first, style, "longitude written first; read the other way round"))
    }
    return(result(NA_real_, NA_real_, "not recognised", "pair of numbers out of range for degrees"))
  }
  NULL
}

# convert_coordinates() converts each record's latitude and longitude text to
# decimal degrees and records how it was written.

convert_coordinates <- function(latitude_text, longitude_text) {
  n <- length(latitude_text)
  latitude <- longitude <- rep(NA_real_, n)
  written_as <- note <- rep(NA_character_, n)
  add_note <- function(i, text) {
    if (!is.na(text)) note[i] <<- if (is.na(note[i])) text else paste(note[i], text, sep = "; ")
  }
  for (i in seq_len(n)) {
    lat_cell <- latitude_text[i]; lon_cell <- longitude_text[i]
    lat_empty <- is.na(lat_cell) || trimws(lat_cell) == ""
    lon_empty <- is.na(lon_cell) || trimws(lon_cell) == ""
    if (lat_empty && lon_empty) next
    
    # A whole position in one cell supplies both coordinates.
    whole <- NULL
    for (cell in unique(c(if (!lat_empty) lat_cell, if (!lon_empty) lon_cell))) {
      whole <- read_whole_position(cell)
      if (!is.null(whole)) break
    }
    if (!is.null(whole)) {
      latitude[i] <- whole$latitude; longitude[i] <- whole$longitude
      written_as[i] <- whole$style; add_note(i, whole$note)
      next
    }
    
    lat_reading <- if (!lat_empty) read_single_coordinate(lat_cell, "latitude") else NULL
    lon_reading <- if (!lon_empty) read_single_coordinate(lon_cell, "longitude") else NULL
    # Compass letters reveal swapped latitude and longitude columns.
    swapped <- !is.null(lat_reading) && !is.null(lon_reading) &&
      lat_reading$compass %in% c("E", "W") && lon_reading$compass %in% c("N", "S")
    if (swapped) {
      lat_reading <- read_single_coordinate(lon_cell, "latitude")
      lon_reading <- read_single_coordinate(lat_cell, "longitude")
      add_note(i, "compass letters show the latitude and longitude columns were the wrong way round; swapped back")
    }
    if (!is.null(lat_reading)) {
      latitude[i] <- lat_reading$value
      if (!is.na(lat_reading$note)) add_note(i, paste("latitude:", lat_reading$note))
      if (lat_reading$compass %in% c("E", "W") && !swapped) {
        add_note(i, "latitude has an E or W letter; left empty")
        latitude[i] <- NA
      }
    }
    if (!is.null(lon_reading)) {
      longitude[i] <- lon_reading$value
      if (!is.na(lon_reading$note)) add_note(i, paste("longitude:", lon_reading$note))
      if (lon_reading$compass %in% c("N", "S") && !swapped) {
        add_note(i, "longitude has an N or S letter; left empty")
        longitude[i] <- NA
      }
    }
    styles <- unique(na.omit(c(lat_reading$style, lon_reading$style)))
    written_as[i] <- if (length(styles) == 0) NA_character_ else paste(styles, collapse = " / ")
  }
  data.frame(latitude = latitude, longitude = longitude, written_as = written_as, note = note,
             stringsAsFactors = FALSE)
}


# 5. Variable dictionary 

dictionary <- read_csv(dictionary_file, col_types = cols(.default = col_character()),
                       na = c("", "NA"), progress = FALSE) %>%
  mutate(variable_order = as.integer(variable_order)) %>%
  arrange(variable_order)

if (anyDuplicated(dictionary$variable_name) > 0) stop("A variable name appears twice in the dictionary.")
if (anyDuplicated(dictionary$role) > 0) stop("A role appears twice in the dictionary.")

dictionary$include_list <- lapply(coalesce(dictionary$match_patterns, ""),
                                  function(p) str_trim(str_split(p, ";")[[1]]) %>% .[. != ""])
dictionary$exclude_text <- coalesce(dictionary$exclude_patterns, "")

# variable_for() returns the variable with a given role, so that variables can
# be renamed in the dictionary.
variable_for <- function(role_name) {
  hit <- dictionary$variable_name[dictionary$role %in% role_name]
  if (length(hit) == 0) NA_character_ else hit[1]
}

kind_of <- function(variable) {
  row <- dictionary[dictionary$variable_name == variable, ]
  type <- tolower(row$data_type)
  unit <- tolower(coalesce(row$unit_or_values, ""))
  if (str_detect(unit, "email")) return("email")
  if (type %in% c("numeric", "integer")) return("number")
  if (type %in% c("datetime", "date")) return("date")
  if (str_detect(unit, "present/absent")) return("presence")
  if (type == "boolean" || str_detect(unit, "^yes/no")) return("yes_no")
  "text"
}
dictionary$kind <- vapply(dictionary$variable_name, kind_of, character(1))
size_roles <- c("circumference", "circumference_2", "dbh", "dbh_2", "height", "crown_1", "crown_2")


# Column names

# observation_in_heading() tells which tree a column belongs to on forms that
# record a baobab and a neighbour.

observation_in_heading <- function(tidy_heading) {
  if (str_detect(tidy_heading, "\\bboth\\b")) return(NA_character_)
  if (str_detect(tidy_heading, "\\bnon ?baobabs?\\b")) return("non-baobab")
  if (str_detect(tidy_heading, "\\bbaobabs?\\b")) return("baobab")
  NA_character_
}

# match_heading() returns the dictionary variable of a column name and the
# priority of the pattern that matched.

heading_cache <- new.env()
match_heading <- function(heading) {
  key <- coalesce(as.character(heading), "<NA>")
  if (!is.null(heading_cache[[key]])) return(heading_cache[[key]])
  result <- match_heading_uncached(heading)
  assign(key, result, envir = heading_cache)
  result
}
match_heading_uncached <- function(heading) {
  tidy <- tidy_text(heading)
  none <- list(variable = NA_character_, priority = NA_integer_)
  if (is.na(tidy) || tidy == "") return(none)
  for_matching <- str_squish(str_replace_all(tidy, "\\bnon ?baobabs?\\b|\\bbaobabs?\\b", " "))
  if (for_matching == "") for_matching <- tidy
  for (i in seq_len(nrow(dictionary))) {
    excluded <- dictionary$exclude_text[i] != "" &&
      str_detect(for_matching, regex(dictionary$exclude_text[i]))
    if (excluded) next
    patterns <- dictionary$include_list[[i]]
    if (length(patterns) == 0) next
    hits <- vapply(patterns, function(p) str_detect(for_matching, regex(p)), logical(1))
    if (any(hits)) return(list(variable = dictionary$variable_name[i], priority = which(hits)[1]))
  }
  none
}

# unit_in_heading() returns the unit written in a column name, preferring units
# in brackets.

unit_in_heading <- function(heading) {
  text <- str_to_lower(coalesce(heading, ""))
  text <- str_remove_all(text, "1[.,]3 ?m(etres|eters)?( above ground)?|breast height|\\(m2\\)|m2")
  bracketed <- str_match_all(text, "[\\(\\[]\\s*(?:in\\s+)?(cm|mm|m|metres?|meters?|centimet(?:re|er)s?)\\s*[\\)\\]]")[[1]]
  found <- if (nrow(bracketed) > 0) bracketed[1, 2] else
    str_match(text, "(?:^|[^a-z])(cm|mm|m)(?:$|[^a-z0-9])")[1, 2]
  case_when(is.na(found) ~ NA_character_,
            str_detect(found, "^cm|^centi") ~ "cm",
            found == "mm" ~ "mm",
            TRUE ~ "m")
}

looks_like_person_name <- function(x) {
  !is.na(x) & !str_detect(x, "[0-9@]") &
    str_detect(x, "^\\s*[A-Za-z\u00C0-\u017F'.-]+(\\s+[A-Za-z\u00C0-\u017F'.-]+){1,3}\\s*$") &
    !tidy_text(x) %in% not_a_name_words
}


# Reading files 

# read_any_file() returns one grid of text per sheet, so that values are not
# changed on reading.

read_any_file <- function(path) {
  extension <- tolower(tools::file_ext(path))
  # GPX and KML layers become tables of point name, description, coordinates,
  # elevation and time.
  if (extension %in% c("gpx", "kml")) {
    layers <- sf::st_layers(path)$name
    grids <- list()
    for (layer in layers) {
      points <- tryCatch(sf::st_read(path, layer = layer, quiet = TRUE), error = function(e) NULL)
      if (is.null(points) || nrow(points) == 0) next
      points <- points[sf::st_geometry_type(points) %in% c("POINT", "MULTIPOINT"), ]
      if (nrow(points) == 0) next
      xy <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(sf::st_transform(points, 4326))))
      fields <- sf::st_drop_geometry(points)
      names(fields) <- tolower(names(fields))
      take <- function(candidates) { hit <- intersect(candidates, names(fields)); if (length(hit)) as.character(fields[[hit[1]]]) else NA_character_ }
      table <- data.frame(`Tree ID` = take(c("name")), Comments = take(c("desc", "description", "cmt")),
                          Latitude = as.character(xy[, 2]), Longitude = as.character(xy[, 1]),
                          Elevation = take(c("ele", "elevation", "altitude")), Date = take(c("time", "timestamp", "date")),
                          check.names = FALSE, stringsAsFactors = FALSE)
      extra <- fields[, setdiff(names(fields), c("name", "desc", "description", "cmt", "ele", "elevation",
                                                 "altitude", "time", "timestamp", "date")), drop = FALSE]
      extra <- extra[, vapply(extra, function(x) any(!is.na(x) & as.character(x) != ""), logical(1)), drop = FALSE]
      if (ncol(extra) > 0) table <- cbind(table, as.data.frame(lapply(extra, as.character), check.names = FALSE))
      grids[[layer]] <- rbind(names(table), as.matrix(table))
    }
    if (length(grids) == 0) stop("no points found in the file")
    return(grids)
  }
  if (extension %in% c("xlsx", "xlsm", "xls")) {
    sheets <- excel_sheets(path)
    grids <- lapply(sheets, function(sheet) {
      table <- suppressMessages(read_excel(path, sheet = sheet, col_names = FALSE,
                                           col_types = "text", .name_repair = "minimal"))
      as.matrix(table)
    })
    names(grids) <- sheets
    return(grids)
  }
  # CSV and TSV files are read as UTF-8, or as Latin-1 if that fails.
  delimiter <- if (extension == "tsv") "\t" else ","
  encoding <- if (is.na(guess_encoding(path)$encoding[1])) "UTF-8" else guess_encoding(path)$encoding[1]
  if (!encoding %in% c("UTF-8", "ASCII")) encoding <- "ISO-8859-1"
  table <- read_delim(path, delim = delimiter, col_names = FALSE, col_types = cols(.default = col_character()),
                      locale = locale(encoding = encoding), progress = FALSE, na = character(),
                      name_repair = "minimal")
  list("(single table)" = as.matrix(table))
}

# tidy_grid() trims spaces, empties blank cells and drops headings added by R.

tidy_grid <- function(grid) {
  if (length(grid) == 0) return(matrix(character(0), 0, 0))
  grid[] <- str_squish(str_replace_all(grid, "\u00a0", " "))
  grid[grid %in% c("", "NA", "N/A", "n/a", "#N/A", "-")] <- NA
  grid[!is.na(grid) & str_detect(grid, "^\\.\\.\\.\\d+$")] <- NA
  grid
}


# Sheet layout

# find_constant_columns() sets aside columns that repeat one value on every
# row, such as an observer's name.

find_constant_columns <- function(grid) {
  n_rows <- nrow(grid)
  result <- list()
  for (j in seq_len(ncol(grid))) {
    values <- grid[, j]
    filled <- which(!is.na(values))
    if (length(filled) < max(5, 0.5 * n_rows)) next
    counts <- sort(table(values[filled]), decreasing = TRUE)
    main_value <- names(counts)[1]
    if (counts[1] / length(filled) < 0.9 || counts[1] < 5) next
    if (!is.na(as_number(main_value))) next
    first_row <- min(which(values == main_value))
    above <- values[seq_len(first_row - 1)]
    above <- above[!is.na(above) & above != main_value]
    heading <- if (length(above) > 0) tail(above, 1) else NA_character_
    # A column headed by its own value is named from its content.
    role <- if (str_detect(main_value, "@")) "observer_email" else
      if (looks_like_person_name(main_value)) "observer_name" else NA_character_
    result[[length(result) + 1]] <- tibble(column = j, value = main_value, heading = heading,
                                           content_role = role)
  }
  bind_rows(tibble(column = integer(), value = character(), heading = character(),
                   content_role = character()), result)
}

score_row <- function(cells) {
  cells <- cells[!is.na(cells)]
  if (length(cells) == 0) return(0)
  sum(vapply(cells, function(cell) !is.na(match_heading(cell)$variable), logical(1)))
}

# has_decimal_numbers() marks rows of data rather than column names.

has_decimal_numbers <- function(cells) {
  numbers <- as_number(cells[!is.na(cells)])
  any(!is.na(numbers) & abs(numbers - round(numbers)) > 1e-9)
}

# find_layout() finds the header rows, group labels, title lines and first data
# row of a sheet.

find_layout <- function(grid, constant_columns) {
  usable <- setdiff(seq_len(ncol(grid)), constant_columns$column)
  rows_to_search <- seq_len(min(30, nrow(grid)))
  scores <- vapply(rows_to_search, function(r) score_row(grid[r, usable]), numeric(1))
  if (length(scores) == 0 || max(scores) < 2) return(NULL)
  header_row <- which.max(scores)
  upper_row <- NA_integer_
  # A second header row, such as "Lat" and "Long" under "GPS".
  below <- header_row + 1
  if (below <= nrow(grid) && score_row(grid[below, usable]) >= 2 &&
      !has_decimal_numbers(grid[below, usable])) {
    upper_row <- header_row
    header_row <- below
  } else if (header_row > 1 && score_row(grid[header_row - 1, usable]) >= 2 &&
             !has_decimal_numbers(grid[header_row - 1, usable])) {
    upper_row <- header_row - 1
  }
  # Upper labels cover the empty cells to their right.
  
  lower <- grid[header_row, ]
  headings <- lower
  if (!is.na(upper_row)) {
    upper <- grid[upper_row, ]
    carried <- upper
    for (j in seq_along(carried)[-1]) {
      if (is.na(carried[j]) && !is.na(carried[j - 1]) && !is.na(lower[j])) carried[j] <- carried[j - 1]
    }
    headings <- ifelse(is.na(lower), carried, ifelse(is.na(carried), lower, paste(carried, lower)))
  }
  for (j in which(is.na(headings))) {
    above <- grid[seq_len(header_row - 1), j]
    above <- above[!is.na(above)]
    if (length(above) > 0) headings[j] <- tail(above, 1)
  }
  # Group labels above the header, such as survey dates, cover the columns to
  # their right.
  
  label_row <- if (is.na(upper_row)) header_row - 1 else upper_row - 1
  group_labels <- rep(NA_character_, ncol(grid))
  if (label_row >= 1) {
    labels <- grid[label_row, ]
    labels[constant_columns$column] <- NA
    if (sum(!is.na(labels[usable])) >= 2 && score_row(labels[usable]) == 0) {
      current <- NA_character_
      for (j in seq_along(labels)) {
        if (!is.na(labels[j])) current <- labels[j]
        group_labels[j] <- current
      }
    }
  }
  first_layout_row <- min(c(header_row, upper_row), na.rm = TRUE)
  title_cells <- as.vector(grid[seq_len(first_layout_row - 1), usable])
  title_cells <- title_cells[!is.na(title_cells)]
  list(header_row = header_row, upper_row = upper_row, headings = headings,
       group_labels = group_labels, title_cells = title_cells, first_data_row = header_row + 1)
}

labelled_values <- function(texts) {
  parts <- str_match(texts, "^\\s*([^:]{2,40}?)\\s*:\\s*(.+?)\\s*$")
  tibble(label = parts[, 2], value = parts[, 3]) %>% filter(!is.na(label), !is.na(value))
}

# names_listed_in_sheet() reads a list of people from a sheet without a table.

names_listed_in_sheet <- function(grid) {
  first_cells <- apply(grid, 1, function(r) { r <- r[!is.na(r)]; if (length(r) == 1) r else NA_character_ })
  start <- which(!is.na(first_cells) &
                   str_detect(tidy_text(first_cells), "collected by|recorded by|observers?$|recorders?$|surveyors?$|field team"))
  if (length(start) == 0) return(character())
  found <- character()
  for (r in seq(start[1] + 1, length.out = 15)) {
    if (r > length(first_cells) || is.na(first_cells[r])) break
    found <- c(found, str_squish(str_remove_all(first_cells[r], "\\(.*?\\)")))
  }
  found[looks_like_person_name(found)]
}


# From sheet to records

# Workbook text, used later to find the country, UTM zone and species of a
# file.
workbook_text <- list()

sheet_to_records <- function(grid, file_label, sheet_name) {
  empty_result <- function(reason) {
    add_log("sheet_log", tibble(file = file_label, sheet = sheet_name, records = 0L, result = reason))
    NULL
  }
  if (nrow(grid) < 2 || ncol(grid) < 2) return(empty_result("empty sheet"))
  constant <- find_constant_columns(grid)
  layout <- find_layout(grid, constant)
  if (is.null(layout)) return(empty_result("no row of recognised column names; no tree records"))
  
  headings <- layout$headings
  n_columns <- ncol(grid)
  
  matches <- lapply(seq_len(n_columns), function(j) {
    is_constant <- j %in% constant$column
    heading <- headings[j]
    if (is_constant) {
      info <- constant[constant$column == j, ]
      heading <- coalesce(info$heading, heading)
      if (!is.na(info$content_role) && (is.na(heading) || heading == info$value)) {
        return(list(variable = variable_for(info$content_role), priority = 1L,
                    heading = coalesce(heading, info$value), how = "recognised from its content"))
      }
    }
    # Survey-form system columns such as "_id" are not tree data.
    if (!is.na(heading) && str_detect(heading, "^_") && !str_detect(tidy_text(heading), " ")) {
      return(list(variable = NA_character_, priority = NA_integer_, heading = heading, how = "form system column"))
    }
    m <- match_heading(heading)
    list(variable = m$variable, priority = m$priority, heading = heading,
         how = if (is.na(m$variable)) "not recognised" else "dictionary pattern")
  })
  column_table <- tibble(
    column = seq_len(n_columns),
    heading = vapply(matches, function(m) coalesce(as.character(m$heading), NA_character_), character(1)),
    variable = vapply(matches, function(m) m$variable, character(1)),
    priority = vapply(matches, function(m) as.integer(m$priority), integer(1)),
    how = vapply(matches, function(m) m$how, character(1)),
    observation = vapply(tidy_text(headings), function(h) if (is.na(h)) NA_character_ else
      observation_in_heading(h), character(1)),
    group = layout$group_labels
  )
  
  data_rows <- seq(layout$first_data_row, length.out = max(0, nrow(grid) - layout$first_data_row + 1))
  usable <- setdiff(seq_len(n_columns), constant$column)
  header_texts <- tidy_text(headings)
  keep_row <- rep(TRUE, length(data_rows))
  section_value <- rep(NA_character_, length(data_rows))
  section_fields <- list()
  current_labels <- list()
  for (k in seq_along(data_rows)) {
    cells <- grid[data_rows[k], usable]
    filled <- which(!is.na(cells))
    if (length(filled) == 0) { keep_row[k] <- FALSE; next }
    if (mean(tidy_text(cells[filled]) %in% header_texts[usable][filled]) >= 0.5 && length(filled) >= 3) {
      keep_row[k] <- FALSE; next
    }
    # Label rows such as "Rep A" apply to the rows below them.
    if (length(filled) == 1 && is.na(as_number(cells[filled]))) {
      text <- cells[filled]
      parts <- labelled_values(text)
      key <- if (nrow(parts) == 1) tidy_text(parts$label) else "section"
      current_labels[[key]] <- text
      keep_row[k] <- FALSE
      next
    }
    if (length(current_labels) > 0) {
      section_value[k] <- paste(unlist(current_labels), collapse = " | ")
      section_fields[[k]] <- current_labels
    }
  }
  if (!any(keep_row)) return(empty_result("column names found but no data rows"))
  
  # Work out how the observations are laid out in the sheet.
  recognised <- column_table %>% filter(!is.na(variable))
  named_observations <- unique(na.omit(recognised$observation))
  repeated_groups <- recognised %>% filter(!is.na(group)) %>% count(variable, group) %>%
    count(variable) %>% filter(n >= 2)
  numbered <- recognised %>%
    mutate(number = str_match(tidy_text(heading), "(?:^| )([12])$")[, 2]) %>%
    filter(!is.na(number))
  numbered_pairs <- numbered %>% distinct(variable, number) %>% count(variable) %>% filter(n == 2)
  observation_of_column <- rep(NA_character_, n_columns)
  how_split <- "one observation per row"
  if (length(named_observations) >= 2) {
    observation_of_column[recognised$column] <- recognised$observation
    how_split <- "baobab and non-baobab tree on one row"
  } else if (nrow(repeated_groups) >= 2) {
    in_group <- recognised$variable %in% repeated_groups$variable & !is.na(recognised$group)
    observation_of_column[recognised$column[in_group]] <- recognised$group[in_group]
    how_split <- "repeated visits side by side"
  } else if (nrow(numbered_pairs) >= 2) {
    pair_rows <- numbered$variable %in% numbered_pairs$variable
    observation_of_column[numbered$column[pair_rows]] <- paste("tree", numbered$number[pair_rows])
    how_split <- "two numbered trees on one row"
  } else if (length(named_observations) == 1) {
    observation_of_column[recognised$column] <- NA_character_
    single_observation <- named_observations
  }
  observations <- unique(na.omit(observation_of_column))
  if (length(observations) == 0) observations <- NA_character_
  
  title_values <- labelled_values(layout$title_cells) %>%
    mutate(variable = vapply(label, function(l) match_heading(l)$variable, character(1))) %>%
    filter(!is.na(variable)) %>% distinct(variable, .keep_all = TRUE)
  
  # Names written anywhere in the sheet, used later for records without a name.
  name_variable <- variable_for("observer_name")
  name_columns <- recognised$column[recognised$variable %in% name_variable]
  names_in_sheet <- unique(c(na.omit(as.vector(grid[data_rows, name_columns])),
                             title_values$value[title_values$variable %in% name_variable]))
  
  rows_used <- data_rows[keep_row]
  records <- list()
  for (observation in observations) {
    own_columns <- if (is.na(observation)) integer() else which(observation_of_column %in% observation)
    shared_columns <- recognised$column[is.na(observation_of_column[recognised$column])]
    columns_here <- recognised %>% filter(column %in% c(own_columns, shared_columns)) %>%
      arrange(variable, priority, column)
    table <- tibble(.rows = length(rows_used))
    unit_notes <- character()
    for (variable in unique(columns_here$variable)) {
      cols <- columns_here$column[columns_here$variable == variable]
      # A second column with the same name passes its values to the partner
      # variable.
      partner <- dictionary$overflow_to[dictionary$variable_name == variable]
      if (length(cols) >= 2 && !is.na(partner) && !partner %in% columns_here$variable) {
        table[[partner]] <- grid[rows_used, cols[2]]
        cols <- cols[-2]
      }
      values <- grid[rows_used, cols, drop = FALSE]
      rule <- coalesce(dictionary$combine_rule[dictionary$variable_name == variable], "first")
      combined <- if (rule == "join") {
        apply(values, 1, function(v) { v <- unique(v[!is.na(v)]); if (length(v) == 0) NA_character_ else paste(v, collapse = "; ") })
      } else if (rule == "any_present") {
        apply(values, 1, function(v) {
          v <- v[!is.na(v)]
          if (length(v) == 0) return(NA_character_)
          meaning <- presence_meaning(v)
          if (any(meaning %in% "yes")) "present" else if (all(meaning %in% "no")) "absent" else v[1]
        })
      } else {
        apply(values, 1, function(v) { v <- v[!is.na(v)]; if (length(v) == 0) NA_character_ else v[1] })
      }
      table[[variable]] <- combined
      if (variable %in% dictionary$variable_name[dictionary$role %in% size_roles]) {
        attr(table[[variable]], "unit") <- unit_in_heading(headings[cols[1]])
      }
    }
    for (k in seq_along(rows_used)) {
      fields <- if (length(section_fields) >= which(keep_row)[k]) section_fields[[which(keep_row)[k]]] else NULL
      if (is.null(fields)) next
      for (key in names(fields)) {
        parts <- labelled_values(fields[[key]])
        if (nrow(parts) == 0) next
        variable <- match_heading(parts$label)$variable
        if (is.na(variable)) next
        if (!variable %in% names(table)) table[[variable]] <- NA_character_
        if (is.na(table[[variable]][k])) table[[variable]][k] <- parts$value
      }
    }
    for (i in seq_len(nrow(title_values))) {
      variable <- title_values$variable[i]
      if (!variable %in% names(table)) table[[variable]] <- NA_character_
      table[[variable]][is.na(table[[variable]])] <- title_values$value[i]
    }
    # A yes/no date column takes its date from the column name.
    date_variable <- variable_for("date")
    if (!is.na(date_variable) && date_variable %in% names(table)) {
      date_col <- columns_here$column[columns_here$variable == date_variable][1]
      answers <- tidy_text(table[[date_variable]])
      if (!is.na(date_col) && mean(answers %in% c(yes_words, no_words), na.rm = TRUE) > 0.8) {
        heading_date <- str_extract(headings[date_col],
                                    "(?i)(jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)[a-z]*\\.? \\d{4}")
        table[[date_variable]] <- ifelse(answers %in% yes_words & !is.na(heading_date), heading_date, NA)
      }
    }
    if (how_split == "repeated visits side by side" && !is.na(date_variable)) {
      table[[date_variable]] <- coalesce(if (date_variable %in% names(table)) table[[date_variable]] else NULL,
                                         rep(observation, length(rows_used)))
    }
    # Keep only rows that hold a tree observation.
    observation_variables <- dictionary$variable_name[dictionary$counts_as_observation == "yes"]
    check_columns <- if (length(own_columns) > 0) {
      intersect(unique(recognised$variable[recognised$column %in% own_columns]), observation_variables)
    } else intersect(names(table), observation_variables)
    holds_observation <- if (length(check_columns) == 0) rep(FALSE, length(rows_used)) else
      rowSums(!is.na(as.data.frame(table[check_columns]))) > 0
    table <- table[holds_observation, , drop = FALSE]
    if (nrow(table) == 0) next
    table$source_file <- file_label
    table$source_sheet <- sheet_name
    table$source_row <- rows_used[holds_observation]
    table$source_section <- section_value[keep_row][holds_observation]
    table$source_column_group <- observation
    table$observation_kind <- if (how_split == "baobab and non-baobab tree on one row") observation else
      if (exists("single_observation", inherits = FALSE)) single_observation else NA_character_
    records[[length(records) + 1]] <- table
  }
  records <- bind_rows(records)
  add_log("column_log", column_table %>% mutate(file = file_label, sheet = sheet_name, header_row = layout$header_row,
                                                .before = 1))
  add_log("sheet_log", tibble(file = file_label, sheet = sheet_name, records = nrow(records),
                              result = paste0("header on row ", layout$header_row,
                                              if (!is.na(layout$upper_row)) paste0(" (with row ", layout$upper_row, ")") else "",
                                              "; ", how_split)))
  list(records = records, title_cells = c(layout$title_cells, unique(na.omit(section_value))),
       names_in_sheet = names_in_sheet)
}

# presence_meaning() reads presence answers, where counts above zero mean yes.
presence_meaning <- function(values) {
  tidy <- tidy_text(values)
  number <- as_number(values)
  case_when(
    is.na(values) ~ NA_character_,
    !is.na(number) & number > 0 ~ "yes",
    !is.na(number) & number == 0 ~ "no",
    tidy %in% yes_words ~ "yes",
    tidy %in% no_words | str_detect(tidy, "^no\\b|\\bno leaves\\b|^none\\b") ~ "no",
    tidy %in% unknown_words ~ NA_character_,
    TRUE ~ "yes"
  )
}


# Reading every file 

file_list <- list.files(raw_input_folder, pattern = "\\.(csv|tsv|xlsx|xlsm|xls|gpx|kml)$",
                        recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
file_list <- file_list[!startsWith(basename(file_list), "~$")]
if (length(file_list) == 0) stop("No CSV, Excel, GPX or KML files were found in ", raw_input_folder, call. = FALSE)

# The data collection method of a file is the name of its folder directly
# under the raw input folder, exactly as written.
relative_path <- substring(normalizePath(file_list, winslash = "/"),
                           nchar(normalizePath(raw_input_folder, winslash = "/")) + 2)
collection_method_of_path <- ifelse(str_detect(relative_path, "/"), str_extract(relative_path, "^[^/]+"), NA)
if (anyNA(collection_method_of_path)) {
  stop("These files are not in a data collection method folder. Move each into one of the folders in ",
       raw_input_folder, ":\n  ", paste(relative_path[is.na(collection_method_of_path)], collapse = "\n  "),
       call. = FALSE)
}

# Identical copies of a file are read once.
fingerprints <- tools::md5sum(file_list)
copy_of <- match(fingerprints, fingerprints)
for (i in which(copy_of != seq_along(file_list))) {
  note <- if (collection_method_of_path[i] != collection_method_of_path[copy_of[i]]) {
    paste0(" (in folder ", collection_method_of_path[i], "; read from folder ",
           collection_method_of_path[copy_of[i]], ")")
  } else ""
  if (note != "") warning("The same file is in two data collection method folders: ",
                          relative_path[i], " and ", relative_path[copy_of[i]], call. = FALSE)
  add_log("sheet_log", tibble(file = basename(file_list[i]), sheet = NA, records = 0L,
                              result = paste0("skipped: identical copy of ", basename(file_list[copy_of[i]]), note)))
}
keep_file <- copy_of == seq_along(file_list)
file_list <- file_list[keep_file]
collection_method_of_path <- collection_method_of_path[keep_file]

# Records are linked to their file by its name, so two different files may not
# share a name.
repeated_names <- unique(basename(file_list)[duplicated(basename(file_list))])
if (length(repeated_names) > 0) {
  stop("Different files share these names; please rename one of each: ",
       paste(repeated_names, collapse = ", "), call. = FALSE)
}
collection_method_of_file <- setNames(collection_method_of_path, basename(file_list))
message("Files to read: ", length(file_list), " (",
        paste(names(table(collection_method_of_path)), table(collection_method_of_path), sep = ": ", collapse = "; "),
        ")")

all_records <- list()
file_names_listed <- list()
sheet_names_written <- list()
for (path in file_list) {
  file_label <- basename(path)
  message("Reading: ", file_label)
  grids <- tryCatch(read_any_file(path), error = function(e) e)
  if (inherits(grids, "error")) {
    add_log("sheet_log", tibble(file = file_label, sheet = NA, records = 0L,
                                result = paste("could not be opened:", conditionMessage(grids))))
    next
  }
  texts <- c(file_label, names(grids))
  listed <- character()
  for (sheet_name in names(grids)) {
    grid <- tidy_grid(grids[[sheet_name]])
    outcome <- tryCatch(sheet_to_records(grid, file_label, sheet_name), error = function(e) {
      add_log("sheet_log", tibble(file = file_label, sheet = sheet_name, records = 0L,
                                  result = paste("error:", conditionMessage(e))))
      NULL
    })
    if (!is.null(outcome) && nrow(outcome$records) > 0) all_records[[length(all_records) + 1]] <- outcome$records
    if (!is.null(outcome)) {
      texts <- c(texts, outcome$title_cells)
      sheet_names_written[[paste(file_label, sheet_name, sep = " | ")]] <- outcome$names_in_sheet
    }
    if (is.null(outcome) && length(grid) > 0) {
      texts <- c(texts, grid[!is.na(grid) & is.na(as_number(grid))])
      listed <- c(listed, names_listed_in_sheet(grid))
    }
  }
  workbook_text[[file_label]] <- paste(unique(texts), collapse = " ")
  file_names_listed[[file_label]] <- unique(listed)
}

records <- bind_rows(all_records)
if (nrow(records) == 0) stop("No tree records were found. See sheet_summary.csv.", call. = FALSE)
for (variable in dictionary$variable_name) if (!variable %in% names(records)) records[[variable]] <- NA_character_
records$record_number <- seq_len(nrow(records))
records$data_collection_method <- unname(collection_method_of_file[records$source_file])
records_read <- nrow(records)
original_values <- records


# Tidying values

log_values <- function(variable, rows, problem) {
  if (length(rows) == 0) return(invisible(NULL))
  add_log("value_log", tibble(variable_name = variable, record_number = records$record_number[rows],
                              source_file = records$source_file[rows], source_sheet = records$source_sheet[rows],
                              source_row = records$source_row[rows],
                              original_value = as.character(original_values[[variable]][rows]),
                              problem = problem))
}

# clean_numbers() keeps one number per cell; values given as limits (< or >)
# are left empty.
clean_numbers <- function(variable) {
  text <- records[[variable]]
  scientific <- suppressWarnings(as.numeric(str_trim(text)))
  scientific[!str_detect(coalesce(text, ""), "(?i)^\\s*-?[0-9.]+e[-+]?[0-9]+\\s*$")] <- NA
  text[!is.na(scientific)] <- format(scientific[!is.na(scientific)], scientific = FALSE, digits = 15, trim = TRUE)
  text <- str_replace_all(text, "(?<=\\d),(?=\\d{1,2}$)", ".")
  text <- str_remove_all(text, ",")
  censored <- !is.na(text) & str_detect(text, "[<>]")
  found <- str_extract_all(coalesce(text, ""), "-?\\d+(?:\\.\\d+)?|-?\\.\\d+")
  n_found <- lengths(found)
  value <- rep(NA_real_, length(text))
  value[n_found == 1 & !censored] <- as.numeric(unlist(found[n_found == 1 & !censored]))
  log_values(variable, which(censored), "written as a limit (< or >); left empty")
  log_values(variable, which(n_found > 1 & !censored), "more than one number in the cell; left empty")
  log_values(variable, which(!is.na(text) & n_found == 0), "no number found; left empty")
  value
}

clean_dates <- function(variable) {
  text <- str_squish(str_replace_all(records[[variable]], "\\s*/\\s*", "/"))
  result <- rep(NA_character_, length(text))
  number <- as_number(text)
  year_only <- !is.na(number) & number == round(number) & number >= 1900 & number <= 2100
  result[year_only] <- as.character(number[year_only])
  serial <- !is.na(number) & !year_only & number > 2100 & number < 80000
  stamp <- as.POSIXct(round(number[serial] * 86400), origin = "1899-12-30", tz = "UTC")
  result[serial] <- ifelse(format(stamp, "%H:%M:%S") == "00:00:00", format(stamp, "%Y-%m-%d"),
                           format(stamp, "%Y-%m-%dT%H:%M:%S"))
  to_read <- is.na(result) & !is.na(text)
  month_year <- str_match(str_to_lower(text), "(jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)[a-z]*\\.?,? (\\d{4})")
  label_only <- to_read & !is.na(month_year[, 1]) & !str_detect(text, "\\d{1,2}[/.-]\\d{1,2}")
  result[label_only] <- format(as.Date(paste(month_year[label_only, 3], month_year[label_only, 2], "01"),
                                       "%Y %b %d"), "%Y-%m")
  to_read <- is.na(result) & !is.na(text)
  orders <- c("Ymd HMS", "Ymd HM", "Ymd", if (day_first_dates) c("dmY HMS", "dmY HM", "dmY", "dmy") else
    c("mdY HMS", "mdY HM", "mdY", "mdy"), "dbY", "bdY", "bY")
  parsed <- suppressWarnings(parse_date_time(text[to_read], orders = orders, tz = "UTC", quiet = TRUE))
  has_time <- !is.na(parsed) & format(parsed, "%H:%M:%S") != "00:00:00"
  result[to_read] <- ifelse(is.na(parsed), NA, ifelse(has_time, format(parsed, "%Y-%m-%dT%H:%M:%S"),
                                                      format(parsed, "%Y-%m-%d")))
  log_values(variable, which(!is.na(text) & is.na(result)), "date not recognised; left empty")
  ambiguous <- !is.na(text) & str_detect(text, "^(\\d{1,2})/(\\d{1,2})/\\d{2,4}") &
    as.integer(str_match(text, "^(\\d{1,2})")[, 2]) <= 12 & as.integer(str_match(text, "^\\d{1,2}/(\\d{1,2})")[, 2]) <= 12
  log_values(variable, which(ambiguous), paste("day and month could be swapped; read as",
                                               if (day_first_dates) "day first" else "month first"))
  result
}

# Month labels without a year take the year of the previous visit.
fill_missing_years <- function(variable) {
  text <- records[[variable]]
  needs_year <- !is.na(text) & str_detect(str_to_lower(text), "^(early|mid|late|end|beginning)? ?(jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)[a-z]*$")
  for (i in which(needs_year)) {
    same_sheet <- which(records$source_file == records$source_file[i] & records$source_sheet == records$source_sheet[i] &
                          !needs_year & !is.na(text) & str_detect(text, "\\d{4}"))
    if (length(same_sheet) == 0) next
    earlier <- records$source_column_group[same_sheet]
    latest <- str_match(str_to_lower(tail(earlier[!is.na(earlier)], 1)), "(jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)[a-z]* (\\d{4})")
    if (is.na(latest[1])) next
    month_here <- match(str_sub(str_extract(str_to_lower(text[i]), "(jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)"), 1, 3),
                        tolower(month.abb))
    year <- as.integer(latest[3]) + (month_here < match(latest[2], tolower(month.abb)))
    text[i] <- paste(text[i], year)
  }
  text
}

clean_answers <- function(variable) {
  kind <- dictionary$kind[dictionary$variable_name == variable]
  values <- records[[variable]]
  tidy <- tidy_text(values)
  meaning <- if (kind == "presence") presence_meaning(values) else
    case_when(str_detect(coalesce(values, ""), "^_+1") ~ "yes",
              tidy %in% yes_words ~ "yes", tidy %in% no_words ~ "no",
              str_detect(tidy, "^more than one|^two|^2|^yes|multi|__1") ~ "yes",
              str_detect(tidy, "^one|^1 stem|^single|^no ") ~ "no")
  result <- if (kind == "presence") recode(meaning, yes = "present", no = "absent") else meaning
  result[tidy %in% unknown_words] <- "unknown"
  log_values(variable, which(!is.na(values) & is.na(result)), "answer not recognised; left empty")
  result
}

date_variable <- variable_for("date")

# Dates taken from photograph file names, for the files listed in
# files_dated_from_photographs (Section 1). The first run of eight digits that
# forms a valid date (YYYYMMDD) in a photograph name is used; every replaced
# date is logged in value_problems.csv with its original value.

photo_date_variables <- na.omit(c(variable_for("photo"), variable_for("photo_2")))
date_from_photograph <- function(references) {
  digits <- str_match(coalesce(references, ""), "(?<!\\d)(\\d{4})(\\d{2})(\\d{2})\\d{6,}")
  dates <- as.Date(paste(digits[, 2], digits[, 3], digits[, 4], sep = "-"), format = "%Y-%m-%d")
  ifelse(is.na(dates), NA_character_, format(dates, "%Y-%m-%d"))
}
if (!is.na(date_variable)) {
  for (file_name in files_dated_from_photographs) {
    rows <- which(records$source_file == file_name)
    if (length(rows) == 0) {
      warning(file_name, " is listed in files_dated_from_photographs but was not read.", call. = FALSE)
      next
    }
    photo_dates <- rep(NA_character_, length(rows))
    for (variable in photo_date_variables) {
      photo_dates <- coalesce(photo_dates, date_from_photograph(records[[variable]][rows]))
    }
    dated <- !is.na(photo_dates)
    records[[date_variable]][rows[dated]] <- photo_dates[dated]
    log_values(date_variable, rows[dated], "date taken from the photograph file name (files_dated_from_photographs)")
    log_values(date_variable, rows[!dated], "no date in the photograph file name; recorded date kept")
    message("Dates taken from photograph file names in ", file_name, ": ", sum(dated), " of ", length(rows))
  }
}

if (!is.na(date_variable)) {
  records[[date_variable]] <- fill_missing_years(date_variable)
  records[[date_variable]] <- clean_dates(date_variable)
}

for (variable in dictionary$variable_name) {
  kind <- dictionary$kind[dictionary$variable_name == variable]
  role <- dictionary$role[dictionary$variable_name == variable]
  if (role %in% c("latitude", "longitude", "position", "date")) next
  if (kind == "number") records[[variable]] <- clean_numbers(variable)
  if (kind %in% c("yes_no", "presence")) records[[variable]] <- clean_answers(variable)
  if (kind == "email") {
    email <- str_to_lower(str_remove_all(records[[variable]], "\\s|^mailto:"))
    bad <- !is.na(email) & !str_detect(email, "^[^@]+@[^@]+\\.[^@.]+$")
    log_values(variable, which(bad), "not an e-mail address; left empty")
    email[bad] <- NA
    records[[variable]] <- email
  }
  if (kind == "text" && !role %in% c("observer_name")) records[[variable]] <- str_squish(records[[variable]])
}


# Sizes: units, diameter and circumference
# The unit of each size column comes from its name or, failing that, from the
# range of its values.

unit_target <- setNames(tolower(dictionary$unit_or_values), dictionary$variable_name)
sheet_key <- paste(records$source_file, records$source_sheet, sep = " | ")
unit_log <- list()

for (variable in dictionary$variable_name[dictionary$role %in% size_roles]) {
  role <- dictionary$role[dictionary$variable_name == variable]
  target <- unit_target[[variable]]
  for (key in unique(sheet_key)) {
    rows <- which(sheet_key == key & !is.na(records[[variable]]))
    if (length(rows) == 0) next
    named <- unique(na.omit(unlist(lapply(all_records, function(t)
      if (variable %in% names(t) && paste(t$source_file[1], t$source_sheet[1], sep = " | ") == key)
        attr(t[[variable]], "unit") else NULL))))
    unit <- if (length(named) > 0) named[1] else NA_character_
    values <- records[[variable]][rows]
    reason <- "unit in column name"
    if (is.na(unit)) {
      high <- quantile(values, 0.95, na.rm = TRUE)
      middle <- median(values, na.rm = TRUE)
      unit <- if (role %in% c("circumference", "circumference_2") && high <= 40) "m" else
        if (role %in% c("dbh", "dbh_2") && high <= 15) "m" else
          if (role %in% c("height", "crown_1", "crown_2") && middle > 100) "cm" else target
      reason <- if (unit == target) "no unit given; values consistent with the dictionary unit" else
        "no unit given; unit worked out from the size of the values"
    }
    if (is.na(target)) target <- unit
    factor <- if (unit == target) 1 else if (unit == "m" && target == "cm") 100 else
      if (unit == "cm" && target == "m") 0.01 else if (unit == "mm" && target == "cm") 0.1 else
        if (unit == "mm" && target == "m") 0.001 else 1
    records[[variable]][rows] <- values * factor
    unit_log[[length(unit_log) + 1]] <- tibble(file_and_sheet = key, variable_name = variable,
                                               unit_recorded = unit, multiplied_by = factor, reason = reason)
  }
}

# Diameters recorded in metres are converted to centimetres, using the
# circumference or height of the tree where available.
height_var_for_dbh <- variable_for("height")
for (pair in list(c(variable_for("dbh"), variable_for("circumference")),
                  c(variable_for("dbh_2"), variable_for("circumference_2")))) {
  dbh_v <- pair[1]; circ_v <- pair[2]
  if (is.na(dbh_v)) next
  # Columns already converted as a whole are skipped.
  already_converted <- unique(unlist(lapply(unit_log, function(u)
    u$file_and_sheet[u$variable_name == dbh_v & u$multiplied_by != 1])))
  for (key in setdiff(unique(sheet_key), already_converted)) {
    rows <- which(sheet_key == key & !is.na(records[[dbh_v]]))
    if (length(rows) == 0) next
    values <- records[[dbh_v]][rows]
    has_decimals <- abs(values - round(values)) > 1e-9
    to_metres <- rep(FALSE, length(rows)); reason <- rep(NA_character_, length(rows))
    if (all(values < 5) && mean(has_decimals) >= 0.5) {
      to_metres[] <- TRUE
      reason[] <- "all DBH values of this column are below 5 with decimals: recorded in metres; converted to cm"
    } else {
      circ <- if (!is.na(circ_v)) records[[circ_v]][rows] else rep(NA_real_, length(rows))
      height <- if (!is.na(height_var_for_dbh)) records[[height_var_for_dbh]][rows] else rep(NA_real_, length(rows))
      matches_circumference <- !is.na(circ) & values < 20 &
        abs(values * 100 * pi - circ) <= 0.15 * circ
      too_thin_as_cm <- !is.na(height) & height >= 3 & values / height < 0.5 &
        (values * 100 / height) >= 3 & (values * 100 / height) <= 80
      column_middle <- if (any(values >= 5)) median(values[values >= 5]) else NA_real_
      small_without_check <- is.na(circ) & is.na(height) & values <= 3 & has_decimals &
        isTRUE(column_middle >= 30)
      to_metres <- matches_circumference | too_thin_as_cm | small_without_check
      reason <- case_when(
        matches_circumference ~ "DBH in metres (DBH x 100 x pi matches the circumference); converted to cm",
        too_thin_as_cm ~ "DBH in metres (in cm it would be thinner than the tree is tall); converted to cm",
        small_without_check ~ "DBH in metres (value of 3 or less with decimals among values in cm); converted to cm")
    }
    if (any(to_metres)) {
      changed <- rows[to_metres]
      for (why in unique(reason[to_metres])) log_values(dbh_v, changed[reason[to_metres] == why], why)
      records[[dbh_v]][changed] <- records[[dbh_v]][changed] * 100
      unit_log[[length(unit_log) + 1]] <- tibble(file_and_sheet = key, variable_name = dbh_v,
                                                 unit_recorded = paste0("m (", length(changed), " of ", length(rows), " values)"),
                                                 multiplied_by = 100, reason = "DBH values recorded in metres; see value_problems.csv")
    }
  }
}

# DBH is kept as recorded. A column whose DBH equals circumference / pi was
# probably calculated by the contributor; it is logged for review, not removed.
circ_var <- variable_for("circumference"); dbh_var <- variable_for("dbh")
if (!is.na(circ_var) && !is.na(dbh_var)) {
  for (key in unique(sheet_key)) {
    rows <- which(sheet_key == key & !is.na(records[[circ_var]]) & !is.na(records[[dbh_var]]))
    if (length(rows) < 3) next
    ratio <- records[[dbh_var]][rows] * pi / records[[circ_var]][rows]
    if (mean(abs(ratio - 1) < 0.03) >= 0.8) {
      log_values(dbh_var, rows, "DBH equals circumference / pi (probably calculated by the contributor); kept as recorded")
    }
  }
}

# Values outside physical limits are left empty.
for (variable in dictionary$variable_name[!is.na(dictionary$plausible_max)]) {
  low <- as.numeric(dictionary$plausible_min[dictionary$variable_name == variable])
  high <- as.numeric(dictionary$plausible_max[dictionary$variable_name == variable])
  values <- records[[variable]]
  outside <- which(!is.na(values) & (values < low | values > high))
  log_values(variable, outside, paste0("outside the physical limits (", low, " to ", high, "); left empty"))
  records[[variable]][outside] <- NA
  if (dictionary$role[dictionary$variable_name == variable] %in% size_roles) {
    zero <- which(!is.na(records[[variable]]) & records[[variable]] == 0)
    log_values(variable, zero, "size of 0 treated as not measured; left empty")
    records[[variable]][zero] <- NA
  }
}


# 13. Coordinates ----

lat_var <- variable_for("latitude"); lon_var <- variable_for("longitude")
position_var <- variable_for("position"); country_var <- variable_for("country")

sf_use_s2(TRUE)
world <- ne_countries(scale = "medium", returnclass = "sf")
names(world) <- tolower(names(world))
world <- st_make_valid(world[, c("name_long", "name", "admin", "continent")])
africa_land <- st_union(world[world$continent == "Africa", ])

# Countries named in each file, used to check positions.
country_names <- unique(c(world$name_long, world$name, world$admin))
countries_named <- lapply(names(workbook_text), function(f) {
  text <- str_to_lower(paste(workbook_text[[f]],
                             paste(unique(na.omit(records[[country_var]][records$source_file == f])), collapse = " ")))
  text <- str_replace_all(text, "_", " ")
  found <- country_names[vapply(country_names, function(n)
    str_detect(text, fixed(str_to_lower(n))) &&
      str_detect(text, paste0("\\b", str_escape(str_to_lower(n)), "\\b")), logical(1))]
  unique(world$name_long[world$name_long %in% found | world$name %in% found | world$admin %in% found])
})
names(countries_named) <- names(workbook_text)

# share_inside() gives the share of points on land in Africa or in the named
# countries.
points_of <- function(lat, lon) st_as_sf(tibble(x = lon, y = lat), coords = c("x", "y"), crs = 4326)
area_cache <- list()
share_inside <- function(lat, lon, countries = character()) {
  ok <- !is.na(lat) & !is.na(lon) & abs(lat) <= 90 & abs(lon) <= 180
  if (!any(ok)) return(0)
  area_key <- if (length(countries) > 0) paste(sort(countries), collapse = "|") else "Africa"
  if (is.null(area_cache[[area_key]])) {
    area_cache[[area_key]] <<- if (length(countries) > 0) st_union(world[world$name_long %in% countries, ]) else africa_land
  }
  area <- area_cache[[area_key]]
  pts <- points_of(lat[ok], lon[ok])
  inside <- lengths(st_is_within_distance(pts, area, dist = 5000)) > 0
  mean(inside)
}

latitude_text  <- records[[lat_var]]
longitude_text <- records[[lon_var]]
position_text  <- if (!is.na(position_var)) records[[position_var]] else rep(NA_character_, nrow(records))
lat_number <- as_number(latitude_text); lon_number <- as_number(longitude_text)

# Numbers above 1,000 are UTM eastings and northings.
is_grid <- !is.na(lat_number) & !is.na(lon_number) & pmax(abs(lat_number), abs(lon_number)) > 1000
latitude <- rep(NA_real_, nrow(records)); longitude <- latitude
written_as <- rep(NA_character_, nrow(records)); coordinate_note <- written_as

not_grid <- which(!is_grid)
converted <- convert_coordinates(latitude_text[not_grid], longitude_text[not_grid])
latitude[not_grid] <- converted$latitude; longitude[not_grid] <- converted$longitude
written_as[not_grid] <- converted$written_as; coordinate_note[not_grid] <- converted$note
from_position <- which(is.na(latitude) & !is.na(position_text))
if (length(from_position) > 0) {
  whole <- convert_coordinates(position_text[from_position], rep(NA_character_, length(from_position)))
  latitude[from_position] <- whole$latitude; longitude[from_position] <- whole$longitude
  written_as[from_position] <- whole$written_as; coordinate_note[from_position] <- whole$note
}
# Record the decimal places of each position as written, before conversion.
decimals_written <- function(number_text) {
  # Spreadsheet storage noise beyond the 12th significant digit is ignored.
  long <- !is.na(number_text) & nchar(coalesce(str_extract(number_text, "(?<=[.,])\\d+$"), "")) > 9
  number_text[long] <- sub("\\.?0+$", "", formatC(signif(as.numeric(sub(",", ".", number_text[long])), 12),
                                                  format = "f", digits = 12, drop0trailing = FALSE))
  d <- nchar(str_extract(number_text, "(?<=[.,])\\d+$"))
  ifelse(is.na(number_text), NA_integer_, ifelse(is.na(d), 0L, as.integer(d)))
}
equivalent_decimals <- function(text, style, whole_position = FALSE) {
  # Numbers in scientific notation are written out before counting.
  scientific <- !is.na(text) & str_detect(text, "(?i)^\\s*-?[0-9]*[.,]?[0-9]+e[-+]?[0-9]+\\s*$")
  text[scientific] <- sub("\\.?0+$", "", formatC(signif(as.numeric(sub(",", ".", text[scientific])), 12),
                                                 format = "f", digits = 12))
  numbers <- str_extract_all(coalesce(text, ""), "\\d+(?:[.,]\\d+)?")
  index <- ifelse(str_detect(coalesce(style, ""), "seconds"), 3L,
                  ifelse(str_detect(coalesce(style, ""), "minutes|NMEA"), 2L, 1L))
  pick <- mapply(function(n, i, packed) {
    if (length(n) == 0) return(NA_character_)
    if (!whole_position || packed) return(tail(n, 1))
    if (length(n) >= i) n[i] else tail(n, 1)
  }, numbers, index, str_detect(coalesce(style, ""), "packed|NMEA"))
  d <- decimals_written(pick)
  resolution <- ifelse(str_detect(coalesce(style, ""), "seconds"), 10^-d / 3600,
                       ifelse(str_detect(coalesce(style, ""), "minutes|NMEA"), 10^-d / 60, 10^-d))
  result <- floor(-log10(resolution) + 1e-9)
  result[is.na(style) | str_detect(coalesce(style, ""), "not recognised")] <- NA
  as.integer(result)
}
coordinate_decimals <- rep(NA_integer_, nrow(records))
separate <- not_grid[!not_grid %in% from_position]
coordinate_decimals[separate] <- pmin(equivalent_decimals(latitude_text[separate], written_as[separate]),
                                      equivalent_decimals(longitude_text[separate], written_as[separate]), na.rm = TRUE)
single_cell <- separate[is.na(longitude_text[separate]) & !is.na(latitude_text[separate])]
coordinate_decimals[single_cell] <- equivalent_decimals(latitude_text[single_cell], written_as[single_cell], TRUE)
coordinate_decimals[from_position] <- equivalent_decimals(position_text[from_position], written_as[from_position], TRUE)
coordinate_decimals[is_grid] <- 5L
mgrs <- str_match(coalesce(written_as, ""), "MGRS grid reference \\(precision (\\d+) m\\)")[, 2]
coordinate_decimals[!is.na(mgrs)] <- as.integer(floor(-log10(as.numeric(mgrs[!is.na(mgrs)]) / 111320) + 1e-9))

impossible <- !is.na(latitude) & (abs(latitude) > 90 | abs(longitude) > 180)
coordinate_note[impossible] <- "values out of range for degrees; left empty"
latitude[impossible] <- NA; longitude[impossible] <- NA

# Swapped or unsigned coordinates are corrected when the corrected reading puts
# the trees in the named country.
for (key in unique(sheet_key)) {
  rows <- which(sheet_key == key & !is.na(latitude) & !is_grid)
  if (length(rows) == 0) next
  named <- countries_named[[records$source_file[rows[1]]]]
  as_recorded <- share_inside(latitude[rows], longitude[rows], named)
  if (as_recorded >= 0.8) next
  options <- expand.grid(swap = c(FALSE, TRUE), lat_sign = c(1, -1), lon_sign = c(1, -1))
  options$share <- apply(options, 1, function(o) {
    y <- if (o[["swap"]]) longitude[rows] else latitude[rows]
    x <- if (o[["swap"]]) latitude[rows] else longitude[rows]
    share_inside(y * o[["lat_sign"]], x * o[["lon_sign"]], named)
  })
  best <- options[which.max(options$share), ]
  if (best$share >= 0.8 && best$share > as_recorded) {
    y <- if (best$swap) longitude[rows] else latitude[rows]
    x <- if (best$swap) latitude[rows] else longitude[rows]
    change <- paste(c(if (best$swap) "latitude and longitude columns swapped",
                      if (best$lat_sign < 0) "latitude sign reversed",
                      if (best$lon_sign < 0) "longitude sign reversed"), collapse = "; ")
    add_log("coordinate_log", tibble(file_and_sheet = key, records = length(rows),
                                     action = change,
                                     reason = paste0(round(100 * as_recorded), "% of trees were in ",
                                                     if (length(named) > 0) paste(named, collapse = "/") else "Africa",
                                                     " as recorded, ", round(100 * best$share), "% after the change")))
    latitude[rows] <- y * best$lat_sign; longitude[rows] <- x * best$lon_sign
    coordinate_note[rows] <- coalesce(paste0(coordinate_note[rows], "; "), "") %>% paste0(change)
  }
}
coordinate_note <- str_remove(coordinate_note, "^NA; ")

# Grid zones come from the file, the setting above, or the zone that puts the
# trees in the named country.
reference <- which(!is.na(latitude))
for (key in unique(sheet_key[is_grid])) {
  rows <- which(sheet_key == key & is_grid)
  a <- abs(lat_number[rows]); b <- abs(lon_number[rows])
  northing <- ifelse(a > b, a, b); easting <- ifelse(a > b, b, a)
  southern <- median(northing) > 5e6
  file_label <- records$source_file[rows[1]]
  text <- workbook_text[[file_label]]
  zone_text <- str_match(text, "(?i)(?:utm\\s*(?:zone)?\\s*|zone\\s*)(\\d{1,2})\\s*([C-X])?\\b")
  zone <- NA_integer_; zone_reason <- NA_character_
  if (!is.na(zone_text[1])) {
    zone <- as.integer(zone_text[2]); zone_reason <- paste0("zone written in the file (", zone_text[1], ")")
  } else if (!is.na(utm_zone_when_not_stated)) {
    zone <- as.integer(str_extract(utm_zone_when_not_stated, "\\d+"))
    southern <- !str_detect(toupper(utm_zone_when_not_stated), "N$")
    zone_reason <- "zone from the setting utm_zone_when_not_stated"
  } else {
    named <- countries_named[[file_label]]
    candidates <- tibble(zone = 26:42)
    centre <- lapply(candidates$zone, function(z) utm_to_degrees(median(easting), median(northing), z, southern))
    candidates$lat <- vapply(centre, function(p) p$latitude, numeric(1))
    candidates$lon <- vapply(centre, function(p) p$longitude, numeric(1))
    candidates$on_land <- vapply(seq_len(nrow(candidates)), function(i)
      share_inside(candidates$lat[i], candidates$lon[i], named) == 1, logical(1))
    candidates$km_to_other_records <- vapply(seq_len(nrow(candidates)), function(i) {
      if (length(reference) == 0) return(NA_real_)
      min(metres_apart(candidates$lat[i], candidates$lon[i], latitude[reference], longitude[reference])) / 1000
    }, numeric(1))
    usable_zones <- candidates %>% filter(on_land) %>% arrange(km_to_other_records)
    if (nrow(usable_zones) > 0) {
      zone <- usable_zones$zone[1]
      zone_reason <- paste0("zone not stated; zone ", zone, " places the trees on land",
                            if (length(named) > 0) paste0(" in ", paste(named, collapse = "/")) else "",
                            ", ", round(usable_zones$km_to_other_records[1]), " km from the nearest other record",
                            if (nrow(usable_zones) > 1) paste0(" (next best: zone ", usable_zones$zone[2], ", ",
                                                               round(usable_zones$km_to_other_records[2]), " km)") else "")
    }
  }
  if (is.na(zone)) {
    coordinate_note[rows] <- "grid coordinates without a zone that could be worked out; left empty"
    add_log("coordinate_log", tibble(file_and_sheet = key, records = length(rows), action = "left empty",
                                     reason = "grid coordinates; zone could not be found"))
    next
  }
  position <- utm_to_degrees(easting, northing, zone, southern)
  latitude[rows] <- position$latitude; longitude[rows] <- position$longitude
  written_as[rows] <- paste0("UTM zone ", zone, if (southern) "S" else "N")
  coordinate_note[rows] <- paste0("converted from UTM ", zone, if (southern) "S" else "N",
                                  ifelse(northing == a, "", "; easting and northing columns were the wrong way round"))
  add_log("coordinate_log", tibble(file_and_sheet = key, records = length(rows),
                                   action = paste0("converted from UTM zone ", zone, if (southern) "S" else "N"),
                                   reason = zone_reason))
}

# Half a coordinate pair is removed.
half_pair <- xor(is.na(latitude), is.na(longitude))
coordinate_note[half_pair] <- paste0(coalesce(coordinate_note[half_pair], "only one coordinate could be read"),
                                     "; the pair was left empty")
latitude[half_pair] <- NA; longitude[half_pair] <- NA

# A tree far from the rest of a compact plot is treated as a typing error.
for (key in unique(sheet_key)) {
  rows <- which(sheet_key == key & !is.na(latitude))
  if (length(rows) < 10) next
  km <- metres_apart(latitude[rows], longitude[rows], median(latitude[rows]), median(longitude[rows])) / 1000
  if (mean(km <= 20) < 0.9) next
  for (k in which(km > 50)) {
    near <- sum(metres_apart(latitude[rows[k]], longitude[rows[k]], latitude[rows[-k]], longitude[rows[-k]]) <= 2000)
    if (isTRUE(near >= 2)) next
    coordinate_note[rows[k]] <- paste0(round(km[k]), " km from the other trees of this plot (probably a typing error); left empty")
    add_log("coordinate_log", tibble(file_and_sheet = key, records = 1L, action = "left empty",
                                     reason = coordinate_note[rows[k]]))
    latitude[rows[k]] <- NA; longitude[rows[k]] <- NA
  }
}

records[[lat_var]] <- round(latitude, 6)
records[[lon_var]] <- round(longitude, 6)
records$coordinate_format <- written_as
records$coordinate_note <- coordinate_note
records$coordinate_decimals <- ifelse(is.na(latitude), NA_integer_, coordinate_decimals)


# Records of one tree
# Records of one tree share their position, elevation, species and setting.
# A tree is identified by its file, sheet, section and identifier, because a
# file holding several plots often numbers the trees of each plot from 1.
# Only in the files listed in files_with_trees_shared_across_sheets (Section
# 1) is a tree matched across the sheets of its file.

tree_var <- variable_for("tree_id")
shared_details <- list(c(lat_var, lon_var), variable_for("elevation"), variable_for("species"),
                       variable_for("site_setting"))
across_sheets <- records$source_file %in% files_with_trees_shared_across_sheets
tree_key <- ifelse(is.na(records[[tree_var]]), NA,
                   ifelse(across_sheets,
                          paste(records$source_file, tidy_text(records[[tree_var]]), sep = " | "),
                          paste(records$source_file, records$source_sheet, coalesce(records$source_section, ""),
                                tidy_text(records[[tree_var]]), sep = " | ")))
copy_count <- 0
for (group in shared_details) {
  group <- group[!is.na(group)]
  complete <- rowSums(is.na(as.data.frame(records[group]))) == 0
  empty <- rowSums(!is.na(as.data.frame(records[group]))) == 0
  # A detail is copied only where every complete record of the tree agrees,
  # so an identifier given to two different trees passes nothing on.
  donor_values <- do.call(paste, c(lapply(records[group], as.character), sep = " | "))
  values_per_key <- tapply(donor_values[complete & !is.na(tree_key)], tree_key[complete & !is.na(tree_key)],
                           function(v) length(unique(v)))
  conflicting_keys <- names(values_per_key)[values_per_key > 1]
  donor <- match(tree_key, ifelse(complete, tree_key, NA))
  take <- which(!is.na(tree_key) & empty & !is.na(donor) & !tree_key %in% conflicting_keys)
  for (v in group) records[[v]][take] <- records[[v]][donor[take]]
  if (length(take) > 0 && lat_var %in% group) {
    records$coordinate_note[take] <- "position copied from another record of the same tree"
    records$coordinate_decimals[take] <- records$coordinate_decimals[donor[take]]
  }
  copy_count <- copy_count + length(take)
}


# Species and tree type
# The species comes from a species column, the tree type, the file or the
# default set above.

species_var <- variable_for("species"); type_var <- variable_for("tree_type")
records$species_source <- ifelse(!is.na(records[[species_var]]), "recorded", NA_character_)
mentions_baobab <- vapply(records$source_file, function(f)
  str_detect(str_to_lower(workbook_text[[f]]), "baoba?b|baobob|adansonia"), logical(1))
sheet_has_species <- as.logical(ave(as.numeric(!is.na(records[[species_var]])), sheet_key, FUN = max))
is_baobab_observation <- records$observation_kind %in% "baobab"
is_other_observation <- records$observation_kind %in% "non-baobab"

# A tree whose recorded type is non-baobab is never given a default species.
recorded_non_baobab <- str_detect(coalesce(tidy_text(records[[type_var]]), ""), "\\bnon")

fill <- is.na(records[[species_var]]) & is_baobab_observation
records[[species_var]][fill] <- "Adansonia digitata"; records$species_source[fill] <- "recorded as the baobab of a paired form"
fill <- is.na(records[[species_var]]) & !is_other_observation & !recorded_non_baobab & !sheet_has_species &
  mentions_baobab
records[[species_var]][fill] <- "Adansonia digitata"; records$species_source[fill] <- "inferred: the file is about baobabs"
if (!is.na(species_when_not_recorded)) {
  fill <- is.na(records[[species_var]]) & !is_other_observation & !recorded_non_baobab & !sheet_has_species
  records[[species_var]][fill] <- species_when_not_recorded
  records$species_source[fill] <- "assumed: no species recorded in the file (setting species_when_not_recorded)"
}
recorded_type <- tidy_text(records[[type_var]])
records[[type_var]] <- case_when(
  is_baobab_observation ~ "baobab",
  is_other_observation ~ "non-baobab",
  str_detect(recorded_type, "baobab") & !str_detect(recorded_type, "non") ~ "baobab",
  str_detect(recorded_type, "non") ~ "non-baobab",
  str_detect(str_to_lower(records[[species_var]]), "adansonia|baobab|kremetart|mowana") ~ "baobab",
  !is.na(records[[species_var]]) ~ "non-baobab",
  TRUE ~ NA_character_)


# Observer names
# Names come from a name column, then an e-mail address, then the sheet or
# file.

name_var <- variable_for("observer_name"); email_var <- variable_for("observer_email")
name_log <- list()
raw_names <- records[[name_var]]
records$observer_name_source <- ifelse(is_blank(raw_names), NA_character_, "recorded")

# standardise_name() removes digits, notes, titles and trailing codes from a
# name.
standardise_name <- function(x) {
  x <- str_remove_all(x, "\\(.*?\\)|\\[.*?\\]|[0-9.,]+$|\\d")
  x <- str_remove_all(x, "(?i)\\b(dr|prof|mr|mrs|ms|phd|msc)\\b\\.?")
  x <- str_squish(str_remove(str_squish(x), "\\s+[A-Z]{2,3}$"))
  words <- str_split(x, " ")
  vapply(words, function(w) {
    w <- w[w != ""]
    if (length(w) == 0) return(NA_character_)
    fix <- w == toupper(w) & nchar(w) > 2 | w == tolower(w)
    w[fix] <- str_to_title(w[fix])
    particle <- seq_along(w) > 1 & seq_along(w) < length(w) & tolower(w) %in% c("de", "da", "do", "dos", "van", "von", "der", "du", "la", "le")
    w[particle] <- tolower(w[particle])
    paste(w, collapse = " ")
  }, character(1))
}

split_people <- function(x) str_trim(unlist(str_split(x, "\\s*(;|,|&|/|\\band\\b|\\n)\\s*")))

# Test submissions are removed.
is_test <- !is_blank(raw_names) & str_detect(str_to_lower(raw_names), "\\btest(ing)?\\b")
if (any(is_test)) {
  add_log("value_log", tibble(variable_name = name_var, record_number = records$record_number[is_test],
                              source_file = records$source_file[is_test], source_sheet = records$source_sheet[is_test],
                              source_row = records$source_row[is_test], original_value = raw_names[is_test],
                              problem = "test submission; record removed"))
}

not_name <- !is_blank(raw_names) & (str_detect(raw_names, "@|\\d{3,}") | tidy_text(raw_names) %in% not_a_name_words |
                                      str_detect(tidy_text(raw_names), "\\b(dead|baobab|tree)\\b"))
log_values(name_var, which(not_name & !is_test), "not a person's name; removed")
raw_names[not_name] <- NA

people <- lapply(raw_names, function(x) if (is_blank(x)) character() else standardise_name(split_people(x)))
people <- lapply(people, function(p) p[!is.na(p) & p != ""])
people[is_test] <- list(character())

# Two names written together are split when both are known elsewhere.
known_once <- unique(unlist(people))
split_joined <- function(p) {
  unlist(lapply(p, function(n) {
    w <- str_split(n, " ")[[1]]
    if (length(w) < 4) return(n)
    for (cut in 2:(length(w) - 2)) {
      first <- paste(w[1:cut], collapse = " "); second <- paste(w[(cut + 1):length(w)], collapse = " ")
      if (first %in% known_once && second %in% known_once) return(c(first, second))
    }
    n
  }))
}
people <- lapply(people, split_joined)

# Optional table of name corrections.
corrections <- if (file.exists(name_corrections_file)) read_csv(name_corrections_file, col_types = cols(.default = "c"),
                                                                progress = FALSE) else tibble(name_as_recorded = character(), corrected_name = character())
corrections <- corrections %>% filter(!is.na(name_as_recorded))
apply_corrections <- function(p) {
  hit <- match(tidy_text(p), tidy_text(standardise_name(corrections$name_as_recorded)))
  ifelse(is.na(hit), p, corrections$corrected_name[hit])
}
people <- lapply(people, apply_corrections)

# Spelling variants of one person's name are given the most frequent spelling.
name_counts <- sort(table(unlist(people)), decreasing = TRUE)
all_names <- names(name_counts)
name_words <- function(n) { w <- str_split(coalesce(tidy_text(n), ""), " ")[[1]]; w[w != ""] }
# first_and_last() returns the first and last words of a name, ignoring
# initials.
first_and_last <- function(n) {
  w <- name_words(n); w <- w[nchar(w) > 1]
  if (length(w) == 0) return(c(NA_character_, NA_character_))
  c(w[1], tail(w, 1))
}
canonical <- setNames(all_names, all_names)
for (n in rev(all_names)) {
  fl <- first_and_last(n)
  tokens_n <- name_words(n)
  if (length(tokens_n) == 0) next
  for (m in all_names) {
    if (m == n || name_counts[[m]] < name_counts[[n]] || canonical[[m]] != m) next
    gl <- first_and_last(m)
    tokens_m <- name_words(m)
    if (length(tokens_m) == 0) next
    words_ok <- !anyNA(fl) && !anyNA(gl)
    same_person <-
      (words_ok && length(tokens_n) >= 2 && length(tokens_m) >= 2 && fl[2] == gl[2] &&
         1 - stringdist(fl[1], gl[1], method = "jw", p = 0.1) >= 0.80 && substr(fl[1], 1, 1) == substr(gl[1], 1, 1)) ||
      (length(tokens_n) >= 2 && length(tokens_n) < length(tokens_m) && all(tokens_n == tokens_m[seq_along(tokens_n)])) ||
      (words_ok && length(tokens_n) >= 2 && identical(fl, gl)) ||
      (words_ok && length(tokens_n) >= 2 && length(tokens_m) >= 2 && fl[1] == gl[1] &&
         1 - stringdist(fl[2], gl[2], method = "jw", p = 0.1) >= 0.90) ||
      (length(tokens_n) == 2 && nchar(tokens_n[2]) == 1 && tokens_n[1] == tokens_m[1] &&
         startsWith(tail(tokens_m, 1), tokens_n[2])) ||
      (length(tokens_n) == 2 && nchar(tokens_n[1]) == 1 && length(tokens_m) >= 2 &&
         tokens_n[2] == tail(tokens_m, 1) && startsWith(tokens_m[1], tokens_n[1]))
    if (isTRUE(same_person)) { canonical[[n]] <- m; break }
  }
}
people <- lapply(people, function(p) unname(canonical[p]))
full_names_by_file <- tapply(seq_along(people), records$source_file, function(i) unique(unlist(people[i])))
all_full <- unique(unlist(people)); all_full <- all_full[str_count(all_full, " ") >= 1]
for (i in seq_along(people)) {
  p <- people[[i]]
  single <- which(str_count(p, " ") == 0)
  for (s in single) {
    in_file <- full_names_by_file[[records$source_file[i]]]
    similar_first <- function(candidates) {
      first <- tidy_text(word(candidates, 1))
      exact <- candidates[which(first == tidy_text(p[s]))]
      if (length(exact) > 0) return(exact)
      candidates[which(1 - stringdist(first, tidy_text(p[s]), method = "jw", p = 0.1) >= 0.9)]
    }
    in_file <- similar_first(in_file[str_count(in_file, " ") >= 1])
    everywhere <- similar_first(all_full)
    if (length(in_file) == 1) p[s] <- in_file else if (length(everywhere) == 1) p[s] <- everywhere
  }
  people[[i]] <- unique(p)
}
records[[name_var]] <- vapply(people, function(p) if (length(p) == 0) NA_character_ else paste(p, collapse = "; "), character(1))
records[[name_var]][records[[name_var]] == ""] <- NA

# Fill missing names from the same e-mail address elsewhere.
if (!is.na(email_var)) {
  pairs <- records %>% filter(!is.na(.data[[email_var]]), !is.na(.data[[name_var]])) %>%
    count(email = .data[[email_var]], person = .data[[name_var]]) %>% group_by(email) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>% ungroup()
  fill <- which(is.na(records[[name_var]]) & records[[email_var]] %in% pairs$email)
  records[[name_var]][fill] <- pairs$person[match(records[[email_var]][fill], pairs$email)]
  records$observer_name_source[fill] <- "from the same e-mail address elsewhere in the dataset"
}
for (key in names(sheet_names_written)) {
  written <- unique(unname(coalesce(canonical[standardise_name(sheet_names_written[[key]])],
                                    standardise_name(sheet_names_written[[key]]))))
  written <- written[!is.na(written) & !str_detect(str_to_lower(written), "\\btest")]
  if (length(written) != 1) next
  fill <- which(is.na(records[[name_var]]) & sheet_key == key)
  records[[name_var]][fill] <- written
  records$observer_name_source[fill] <- "the only person named in the same sheet"
}
for (level in c("sheet", "file")) {
  group <- if (level == "sheet") sheet_key else records$source_file
  only_name <- tapply(records[[name_var]], group, function(n) { n <- unique(na.omit(n)); if (length(n) == 1) n else NA })
  fill <- which(is.na(records[[name_var]]) & !is.na(only_name[group]))
  records[[name_var]][fill] <- only_name[group][fill]
  records$observer_name_source[fill] <- paste("the only person named in the same", level)
}
for (f in names(file_names_listed)) {
  listed <- standardise_name(file_names_listed[[f]])
  if (length(listed) == 0) next
  listed <- unique(unname(canonical[listed] %>% coalesce(listed)))
  fill <- which(is.na(records[[name_var]]) & records$source_file == f)
  records[[name_var]][fill] <- paste(listed, collapse = "; ")
  records$observer_name_source[fill] <- "listed on a metadata sheet of the file"
}
name_review <- tibble(name_as_recorded = raw_names, name_written = records[[name_var]],
                      source = records$observer_name_source, source_file = records$source_file) %>%
  count(source_file, name_as_recorded, name_written, source, name = "records")

records <- records[!is_test, ]
removed_test_submissions <- sum(is_test)

# Dead trees ----
# Records that mention "Dead baobab" in any field are removed.
original_rows <- original_values[match(records$record_number, original_values$record_number), ]
text_variables <- intersect(dictionary$variable_name[dictionary$kind %in% c("text", "presence", "yes_no")],
                            names(original_rows))
all_text <- do.call(paste, c(lapply(original_rows[text_variables], function(x) coalesce(as.character(x), "")),
                             sep = " | "))
says_dead_baobab <- str_detect(str_to_lower(all_text), "dead\\s*-?\\s*baobabs?")
status_var <- variable_for("tree_status")
recorded_dead <- if (remove_trees_recorded_as_dead && !is.na(status_var)) {
  tidy_text(original_rows[[status_var]]) %in% c("dead", "d", "died", "dead tree")
} else rep(FALSE, nrow(records))
remove_dead <- says_dead_baobab | recorded_dead
removed_records <- tibble(
  reason = ifelse(says_dead_baobab, "a field says \"Dead baobab\"", "status recorded as dead"),
  source_file = records$source_file, source_sheet = records$source_sheet, source_row = records$source_row,
  tree = records[[variable_for("tree_id")]], text_found = all_text)[remove_dead, ]
records <- records[!remove_dead, ]
message("Records removed as dead trees: ", sum(remove_dead))


# Duplicates
# The same tree on several sheets is merged, and identical records are kept
# once.

output_variables <- dictionary$variable_name[dictionary$role != "position"]
compare_variables <- setdiff(output_variables, c(variable_for("feature_id"), variable_for("tree_id")))
filled_count <- rowSums(!is.na(as.data.frame(records[output_variables])))
position_key <- ifelse(is.na(records[[lat_var]]), NA,
                       paste(records$source_file, round(records[[lat_var]], 5), round(records[[lon_var]], 5)))

# agree() is TRUE when two records never disagree on a variable (numbers within
# 1%).
agree <- function(i, j) {
  for (v in compare_variables) {
    x <- records[[v]][i]; y <- records[[v]][j]
    if (is.na(x) || is.na(y)) next
    if (is.numeric(x)) { if (abs(x - y) > 0.01 * max(abs(x), abs(y))) return(FALSE) } else if (x != y) return(FALSE)
  }
  TRUE
}
duplicate_log <- list()
remove <- rep(FALSE, nrow(records))
for (key in unique(na.omit(position_key[duplicated(position_key)]))) {
  rows <- which(position_key == key)
  if (length(unique(records$source_sheet[rows])) < 2) next
  rows <- rows[order(-filled_count[rows])]
  kept <- integer()
  for (i in rows) {
    target <- NA_integer_
    for (k in kept) {
      if (records$source_sheet[k] != records$source_sheet[i] && agree(k, i)) { target <- k; break }
    }
    if (is.na(target)) { kept <- c(kept, i); next }
    for (v in output_variables) if (is.na(records[[v]][target]) && !is.na(records[[v]][i])) records[[v]][target] <- records[[v]][i]
    remove[i] <- TRUE
    duplicate_log[[length(duplicate_log) + 1]] <- tibble(
      check = "same tree on another sheet of the file (merged into the fuller record)",
      source_file = records$source_file[i], source_sheet = records$source_sheet[i], source_row = records$source_row[i],
      kept = paste(records$source_sheet[target], "row", records$source_row[target]))
  }
}
records <- records[!remove, ]
merged_records <- sum(remove)
identical_rows <- duplicated(records[output_variables])
if (any(identical_rows)) {
  duplicate_log[[length(duplicate_log) + 1]] <- tibble(check = "identical to another record (removed)",
                                                       source_file = records$source_file[identical_rows],
                                                       source_sheet = records$source_sheet[identical_rows],
                                                       source_row = records$source_row[identical_rows], kept = NA)
  records <- records[!identical_rows, ]
}


# Saving the results

provenance <- c("species_source", "observer_name_source", "coordinate_format", "coordinate_note", "coordinate_decimals",
                "data_collection_method", "source_file", "source_sheet", "source_row", "source_section",
                "source_column_group")
for (v in output_variables) records[[v]] <- as.character(records[[v]])
for (v in c(lat_var, lon_var)) records[[v]] <- ifelse(is.na(records[[v]]), NA, sprintf("%.6f", as.numeric(records[[v]])))
baobab_tree_dataset <- records %>% arrange(source_file, source_sheet, source_row) %>%
  select(all_of(c(output_variables, na.omit(position_var), provenance))) %>%
  rename(source_species = species_source, source_observer_name = observer_name_source,
         source_coordinate_format = coordinate_format, source_coordinate_note = coordinate_note,
         source_coordinate_decimals = coordinate_decimals,
         source_data_collection_method = data_collection_method)
write_excel_csv(baobab_tree_dataset, harmonised_dataset_file, na = "")

# Record counts at each step, read by step 2 for the data quality table.
harmonisation_counts <- tibble(
  quantity = c("Tree records read from the contributor files", "Removed: test submissions",
               "Removed: dead trees", "Merged: same tree on another sheet of the file",
               "Removed: identical records", "Records in the harmonised dataset"),
  records = c(records_read, removed_test_submissions, sum(remove_dead), merged_records,
              sum(identical_rows), nrow(baobab_tree_dataset)))
write_excel_csv(harmonisation_counts, file.path(harmonised_results_folder, "harmonisation_counts.csv"))

write_excel_csv(bind_rows(column_log), file.path(harmonised_results_folder, "column_matches.csv"), na = "")
write_excel_csv(bind_rows(sheet_log), file.path(harmonised_results_folder, "sheet_summary.csv"), na = "")
write_excel_csv(bind_rows(coordinate_log), file.path(harmonised_results_folder, "coordinate_repairs.csv"), na = "")
write_excel_csv(bind_rows(value_log), file.path(harmonised_results_folder, "value_problems.csv"), na = "")
write_excel_csv(bind_rows(unit_log), file.path(harmonised_results_folder, "size_units.csv"), na = "")
write_excel_csv(name_review, file.path(harmonised_results_folder, "observer_names.csv"), na = "")
write_excel_csv(bind_rows(duplicate_log), file.path(harmonised_results_folder, "duplicates.csv"), na = "")
write_excel_csv(removed_records, file.path(harmonised_results_folder, "removed_dead_trees.csv"), na = "")
coordinate_formats <- baobab_tree_dataset %>% count(source_file, source_coordinate_format, name = "records")
write_excel_csv(coordinate_formats, file.path(harmonised_results_folder, "coordinate_formats.csv"), na = "")
unrecognised <- bind_rows(column_log) %>% filter(is.na(variable), !is.na(heading)) %>%
  count(heading, name = "sheets", sort = TRUE)
write_excel_csv(unrecognised, file.path(harmonised_results_folder, "columns_not_recognised.csv"), na = "")

message("============================================================")
message("Finished in ", round(as.numeric(difftime(Sys.time(), start_time, units = "secs"))), " s")
message("Records saved:        ", nrow(baobab_tree_dataset))
message("With coordinates:     ", sum(!is.na(baobab_tree_dataset[[lat_var]])))
message("With observer name:   ", sum(!is.na(baobab_tree_dataset[[name_var]])))
message("Results folder:       ", harmonised_results_folder)
message("Paths saved for steps 2 to 5 in: ", normalizePath(paths_file, winslash = "/"))
message("============================================================")
