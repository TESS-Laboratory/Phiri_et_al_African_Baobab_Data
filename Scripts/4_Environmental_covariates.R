#-------------------------------------------------------------------------------
# Introduction
#-------------------------------------------------------------------------------


# Extracts, at every unique baobab location in the public data file, mean
# annual precipitation (CHIRPS), mean annual temperature (ERA5-Land),
# potential evapotranspiration and the aridity index (TerraClimate), topsoil
# texture (iSDAsoil), slope (Copernicus DEM GLO-30) and protected-area status
# (WDPA, cleaned with wdpar). Values are cached, so an interrupted run resumes.


# 1. Pipeline folder and setup

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

# pipeline_value() returns one path or setting from pipeline_paths.csv, where
# every path of the pipeline is set by 1_Harmonise_tree_data.R.
pipeline_value <- function(name) {
  value <- pipeline_paths$value[pipeline_paths$name == name]
  if (length(value) == 0) stop(name, " is missing from pipeline_paths.csv. Run 1_Harmonise_tree_data.R again.",
                               call. = FALSE)
  value[1]
}

# Input files
public_data_file <- pipeline_value("public_data_file")   # public BAM data, written by 2_QAQC_BAM_Data.qmd
dictionary_file  <- pipeline_value("dictionary_file")    # BAM_variable_dictionary.csv

# Environmental data: download folder, climate period and optional local
# copies ("" means the data are downloaded).
environment_data_folder           <- pipeline_value("environment_data_folder")
climate_period                    <- pipeline_value("climate_period")
chirps_local_folder               <- pipeline_value("chirps_local_folder")
terraclimate_pet_local_source     <- pipeline_value("terraclimate_pet_local_source")
era5land_temperature_local_source <- pipeline_value("era5land_temperature_local_source")
isda_local_folder                 <- pipeline_value("isda_local_folder")
protected_areas_local_source      <- pipeline_value("protected_areas_local_source")
dem_local_source                  <- if ("dem_local_source" %in% pipeline_paths$name)
  pipeline_value("dem_local_source") else ""  # absent from older pipeline_paths.csv files

# Output folders
environment_results_folder <- pipeline_value("environment_results_folder")  # tables and cached values
figures_folder             <- pipeline_value("figures_folder")              # Figure 3

# Output files
covariates_output     <- file.path(environment_results_folder, "baobab_environmental_covariates.csv")
covariate_summary_output <- file.path(environment_results_folder, "covariate_summary.csv")
sources_output        <- file.path(environment_results_folder, "data_sources_used.csv")
figure_output_stem    <- file.path(figures_folder, "figure_3_environment")  # .pdf and .png

message("Files used by 4_Environmental_covariates.R:",
        "\n  public data file:        ", public_data_file,
        "\n  dictionary:              ", dictionary_file,
        "\n  environmental data:      ", environment_data_folder,
        "\n  climate period:          ", climate_period,
        "\n  results folder:          ", environment_results_folder,
        "\n  figure:                  ", figure_output_stem, ".pdf/.png")

climate_period <- eval(parse(text = climate_period))  # "1991:2020" becomes 1991 to 2020

# Variables are used by dictionary role, so they can be renamed in the dictionary.

dictionary_roles <- utils::read.csv(dictionary_file, colClasses = "character")[, c("variable_name", "role")]
rename_to_roles <- function(table) {
  hit <- names(table) %in% dictionary_roles$variable_name
  names(table)[hit] <- dictionary_roles$role[match(names(table)[hit], dictionary_roles$variable_name)]
  table
}

# CHIRPS versions tried, in this order. The manuscript cites version 3.0, so
# the script stops if it is not available rather than falling back to 2.0
# without notice; add "2.0" only if the manuscript is changed to match.
chirps_versions <- c("3.0")
keep_downloaded_files <- TRUE        # FALSE deletes raster files once extracted
era5land_margin_degrees <- 0.5       # margin around the locations for ERA5-Land requests
cds_user <- "ecmwfr"                 # name under which the CDS API key is stored

# iSDAsoil topsoil layers (0-20 cm); "metadata" reads the back-transformation
# from the iSDAsoil metadata files.
isda_layers <- c(sand = "sand_content", silt = "silt_content", clay = "clay_content")
isda_base_url <- "https://isdasoil.s3.amazonaws.com/soil_data"
isda_back_transformation <- "metadata"

# UNEP aridity classes (lower limits of the aridity index).
aridity_classes <- data.frame(
  class = c("Hyper-arid", "Arid", "Semi-arid", "Dry sub-humid", "Humid"),
  lower = c(0, 0.05, 0.2, 0.5, 0.65)
)

erase_protected_area_overlaps <- FALSE

# Countries with fewer baobabs are not shown separately in Figure 3e; 1
# shows every country with a baobab, as in the published figure.
minimum_baobabs_per_country <- 1

# Height of Figure 3 in millimetres (Nature Portfolio maximum: 247 mm). The
# figure has three rows of two panels.
figure_3_height_mm <- 215

# In the texture triangle of Figure 3f, a class label also gives the share of
# baobabs in the class when the share is at least this many per cent. The
# shares of every class are saved in texture_class_shares.csv.
texture_share_label_minimum <- 5


# 2. Packages 

required_packages <- c("readr", "dplyr", "tidyr", "stringr", "ggplot2", "patchwork", "scales",
                       "sf", "terra", "rnaturalearth", "rnaturalearthdata", "curl", "jsonlite",
                       "ragg", "systemfonts")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                              logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Please install these packages first: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}
# curl, jsonlite and scales are loaded before readr, so that readr's
# parse_date() and col_factor() take precedence over functions of the same
# names.
suppressPackageStartupMessages({
  library(curl)               # downloads
  library(jsonlite)           # iSDAsoil metadata
  library(scales)             # axis labels
  library(ragg)               # PNG output
  library(systemfonts)        # finding Arial or Helvetica
  library(rnaturalearth)      # country boundaries
  library(rnaturalearthdata)  # data used by rnaturalearth
  library(readr)              # reading and writing CSV files
  library(dplyr)              # table manipulation
  library(tidyr)              # reshaping tables
  library(stringr)            # text handling
  library(ggplot2)            # Figure 3
  library(patchwork)          # combining the panels of Figure 3
  library(sf)                 # vector data
  library(terra)              # raster data
})

# need_package() stops with a clear message if a package needed for one data
# source (ecmwfr for ERA5-Land, wdpar for protected areas) is not installed;
# the package is then loaded with library() where it is used.

need_package <- function(package, purpose) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("The package ", package, " is needed to ", purpose, ". Please install it: install.packages(\"",
         package, "\")", call. = FALSE)
  }
}


# 3. Helper functions 

given <- function(value) is.character(value) && length(value) == 1 && !is.na(value) && value != ""

download_file <- function(url, destination, attempts = 3) {
  if (file.exists(destination) && file.size(destination) > 0) return(TRUE)  
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  partial <- paste0(destination, ".part")
  for (attempt in seq_len(attempts)) {
    
    handle <- curl::new_handle(connecttimeout = 60, low_speed_limit = 100, low_speed_time = 300)
    ok <- tryCatch({ curl::curl_download(url, partial, quiet = TRUE, handle = handle); TRUE },
                   error = function(problem) FALSE)
    if (ok && file.exists(partial) && file.size(partial) > 0) {
      file.rename(partial, destination)
      return(TRUE)
    }
    Sys.sleep(5 * attempt)
  }
  if (file.exists(partial)) file.remove(partial)
  FALSE
}

readable_raster <- function(path) {
  file.exists(path) && file.size(path) > 0 &&
    !inherits(try(terra::rast(path), silent = TRUE), "try-error")
}

list_links <- function(url) {
  handle <- curl::new_handle(timeout = 120, connecttimeout = 30, useragent = "BAM-project R script")
  reply <- tryCatch(curl::curl_fetch_memory(url, handle = handle), error = function(problem) problem)
  if (inherits(reply, "error")) {
    message("    could not connect to ", url, ": ", conditionMessage(reply))
    return(NULL)
  }
  if (reply$status_code != 200) {
    message("    the server answered ", reply$status_code, " for ", url)
    return(NULL)
  }
  page <- rawToChar(reply$content)
  links <- regmatches(page, gregexpr('href="[^"?#]+"', page))[[1]]
  links <- gsub('^href="|"$', "", links)
  folders <- grepl("/$", links)
  links <- sub("/$", "", links)
  links <- basename(links)
  links[folders] <- paste0(links[folders], "/")
  unique(links)
}

ungzip <- function(path, folder, remove_original) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  output <- file.path(folder, sub("\\.gz$", "", basename(path)))
  input <- gzfile(path, "rb")
  writeBin(readBin(input, "raw", n = 2e9), output)
  close(input)
  if (remove_original) file.remove(path)
  output
}

points_as_vector <- function(longitude, latitude, raster) {
  points <- terra::vect(data.frame(x = longitude, y = latitude), geom = c("x", "y"), crs = "EPSG:4326")
  if (!terra::same.crs(points, raster)) points <- terra::project(points, terra::crs(raster))
  points
}

values_at_locations <- function(raster, longitude, latitude) {
  points <- points_as_vector(longitude, latitude, raster)
  values <- as.matrix(terra::extract(raster, points, ID = FALSE))
  missing <- which(rowSums(!is.na(values)) == 0)
  if (length(missing) > 0) {
    xy <- terra::crds(points)[missing, , drop = FALSE]
    centre <- terra::cellFromXY(raster, xy)
    for (k in seq_along(missing)) {
      if (is.na(centre[k])) next
      neighbours <- terra::adjacent(raster, centre[k], directions = 8)
      neighbours <- neighbours[!is.na(neighbours)]
      if (length(neighbours) == 0) next
      neighbour_values <- as.matrix(terra::extract(raster, neighbours))
      has_data <- rowSums(!is.na(neighbour_values)) > 0
      if (!any(has_data)) next
      neighbour_xy <- terra::xyFromCell(raster, neighbours[has_data])
      nearest <- which.min((neighbour_xy[, 1] - xy[k, 1])^2 + (neighbour_xy[, 2] - xy[k, 2])^2)
      values[missing[k], ] <- neighbour_values[has_data, , drop = FALSE][nearest, ]
    }
  }
  attr(values, "neighbour_used") <- sum(rowSums(!is.na(values[missing, , drop = FALSE])) > 0)
  values
}

cds_key_stored <- function() {
  if (nzchar(Sys.getenv("ecmwfr_PAT"))) return(TRUE)
  !inherits(try(ecmwfr::wf_get_key(user = cds_user), silent = TRUE), "try-error")
}

ensure_cds_key <- function() {
  if (cds_key_stored()) return(invisible(TRUE))
  instructions <- paste0(
    "No Climate Data Store (CDS) API key is stored on this computer.\n",
    "  1. Log in at https://cds.climate.copernicus.eu and open your profile page;\n",
    "  2. copy the 'Personal Access Token' (the API key);\n",
    "  3. in the R console run: ecmwfr::wf_set_key(key = \"paste-your-key-here\")\n",
    "  4. accept the ERA5-Land licence under 'Download data' on the dataset page;\n",
    "  then run this script again.")
  if (!interactive()) stop(instructions, call. = FALSE)
  message(instructions)
  key <- if (requireNamespace("getPass", quietly = TRUE)) {
    library(getPass)  # optional: hides the key as it is typed
    getPass::getPass(msg = "Or paste your CDS API key now (leave empty to stop): ")
  } else readline("Or paste your CDS API key now (leave empty to stop): ")
  key <- trimws(if (is.null(key)) "" else key)
  if (key == "") stop("No API key given; ERA5-Land cannot be downloaded.", call. = FALSE)
  stored <- try(ecmwfr::wf_set_key(key = key, user = cds_user), silent = TRUE)
  if (inherits(stored, "try-error") || !cds_key_stored()) {
    
    Sys.setenv(ecmwfr_PAT = key)
    message("  the key could not be saved in the credential store; it is used for this session only")
  } else {
    message("  API key stored; it will be found automatically next time")
  }
  invisible(TRUE)
}

layer_dates <- function(raster, path) {
  dates <- suppressWarnings(as.Date(terra::time(raster)))
  if (length(dates) == terra::nlyr(raster) && !any(is.na(dates))) return(dates)
  seconds <- suppressWarnings(as.numeric(str_match(names(raster), "valid_time=(-?\\d+)")[, 2]))
  if (!any(is.na(seconds))) return(as.Date(as.POSIXct(seconds, origin = "1970-01-01", tz = "UTC")))
  stop("The dates of ", basename(path), " could not be read.", call. = FALSE)
}

days_in_month <- function(year, month) {
  first <- as.Date(sprintf("%04d-%02d-01", year, month))
  as.numeric(as.Date(sprintf("%04d-%02d-01", year + (month == 12), month %% 12 + 1)) - first)
}

cache_path <- function(name) file.path(environment_results_folder, paste0("cache_", name, ".csv"))
read_cache <- function(name, locations) {
  path <- cache_path(name)
  if (!file.exists(path)) return(NULL)
  cached <- read_csv(path, col_types = cols(latitude = col_character(), longitude = col_character(),
                                            .default = col_guess()), progress = FALSE)
  if (!all(paste(locations$latitude_text, locations$longitude_text) %in%
           paste(cached$latitude, cached$longitude))) return(NULL)
  message("  using values saved by an earlier run (", basename(path), ")")
  cached
}
write_cache <- function(name, table) write_csv(table, cache_path(name), na = "NA")

# read_saved_values() returns a whole cache file, so that sources can be
# completed piece by piece.
read_saved_values <- function(name) {
  path <- cache_path(name)
  if (!file.exists(path)) return(NULL)
  read_csv(path, col_types = cols(latitude = col_character(), longitude = col_character(),
                                  .default = col_guess()), progress = FALSE)
}

# as_transformation() accepts only simple arithmetic and exp, expm1, log and
# log1p.
as_transformation <- function(expression_text) {
  expression_text <- trimws(gsub("np\\.", "", expression_text))
  if (expression_text %in% c("", "none", "None", "null")) expression_text <- "x"
  if (!grepl("^[x0-9 .+*/()-]*$", gsub("expm1|exp|log1p|log", "", expression_text))) {
    stop("Unexpected iSDAsoil back-transformation: ", expression_text, call. = FALSE)
  }
  transform <- function(x) eval(parse(text = expression_text), list(x = x))
  attr(transform, "text") <- expression_text
  transform
}

metadata_transformation <- function(metadata_url) {
  metadata <- tryCatch(jsonlite::fromJSON(metadata_url, simplifyVector = FALSE), error = function(e) NULL)
  if (is.null(metadata)) return(NULL)
  found <- character()
  search <- function(item, name = "") {
    if (is.list(item)) {
      for (key in names(item)) search(item[[key]], key)
      if (is.null(names(item))) for (element in item) search(element, name)
    } else if (is.character(item) && grepl("back", name, ignore.case = TRUE)) {
      found <<- c(found, item)
    }
  }
  search(metadata)
  as_transformation(if (length(found) == 0) "x" else found[1])
}

usda_texture_class <- function(sand, silt, clay) {
  total <- sand + silt + clay
  sand <- 100 * sand / total; silt <- 100 * silt / total; clay <- 100 * clay / total
  case_when(
    is.na(total)                                                     ~ NA_character_,
    silt + 1.5 * clay < 15                                           ~ "Sand",
    silt + 1.5 * clay >= 15 & silt + 2 * clay < 30                   ~ "Loamy sand",
    (clay >= 7 & clay < 20 & sand > 52 & silt + 2 * clay >= 30) |
      (clay < 7 & silt < 50 & silt + 2 * clay >= 30)                 ~ "Sandy loam",
    clay >= 7 & clay < 27 & silt >= 28 & silt < 50 & sand <= 52      ~ "Loam",
    (silt >= 50 & clay >= 12 & clay < 27) | (silt >= 50 & silt < 80 & clay < 12) ~ "Silt loam",
    silt >= 80 & clay < 12                                           ~ "Silt",
    clay >= 20 & clay < 35 & silt < 28 & sand > 45                   ~ "Sandy clay loam",
    clay >= 27 & clay < 40 & sand > 20 & sand <= 45                  ~ "Clay loam",
    clay >= 27 & clay < 40 & sand <= 20                              ~ "Silty clay loam",
    clay >= 35 & sand > 45                                           ~ "Sandy clay",
    clay >= 40 & silt >= 40                                          ~ "Silty clay",
    clay >= 40 & sand <= 45 & silt < 40                              ~ "Clay",
    TRUE                                                             ~ NA_character_
  )
}

sources_used <- tibble(covariate = character(), source = character(), version_or_period = character(),
                       accessed = character(), note = character())
record_source <- function(covariate, source, version, note = NA_character_) {
  sources_used <<- bind_rows(sources_used, tibble(covariate = covariate, source = source,
                                                  version_or_period = version,
                                                  accessed = format(Sys.Date()), note = note))
}


# 4. Baobab locations 

if (!file.exists(public_data_file)) {
  stop("The public data file was not found: ", public_data_file, ". Render 2_QAQC_BAM_Data.qmd first.",
       call. = FALSE)
}
for (folder in c(environment_results_folder, figures_folder, environment_data_folder)) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
}


public_data <- read_csv(public_data_file, col_types = cols(.default = col_character()),
                        na = c("", "NA"), progress = FALSE) %>% rename_to_roles()
locations <- public_data %>%
  filter(tree_type %in% "baobab", !is.na(latitude), !is.na(longitude)) %>%
  group_by(latitude_text = latitude, longitude_text = longitude) %>%
  summarise(records = n(), country = first(country), .groups = "drop") %>%
  mutate(location_id = row_number(),
         latitude = as.numeric(latitude_text), longitude = as.numeric(longitude_text))
if (nrow(locations) == 0) stop("The public data file has no baobab records with coordinates.", call. = FALSE)
message("Unique baobab locations: ", nrow(locations))
years <- climate_period
period_text <- paste0(min(years), "-", max(years))


# 5. Precipitation: CHIRPS monthly totals 

message("\nPrecipitation (CHIRPS):")

location_keys <- paste(locations$latitude_text, locations$longitude_text)

months_needed <- expand.grid(month = 1:12, year = years)[, c("year", "month")]

saved_chirps <- read_saved_values("chirps_monthly")
if (!is.null(saved_chirps)) {
  saved_chirps <- saved_chirps %>% filter(year %in% years)
  complete_months <- saved_chirps %>%
    mutate(key = paste(latitude, longitude)) %>%
    filter(key %in% location_keys) %>%
    distinct(year, month, key) %>%
    count(year, month, name = "locations_saved") %>%
    filter(locations_saved == length(location_keys))
} else {
  complete_months <- tibble(year = integer(), month = integer())
}
months_to_do <- anti_join(months_needed, complete_months, by = c("year", "month"))
message("  months of ", period_text, " already saved for all ", length(location_keys), " locations: ",
        nrow(months_needed) - nrow(months_to_do), " of ", nrow(months_needed))

# Values saved from another CHIRPS version are never mixed with the version
# cited in the manuscript.
saved_version_file <- sub("\\.csv$", "_source.txt", cache_path("chirps_monthly"))
if (!is.null(saved_chirps) && nrow(saved_chirps) > 0 && !given(chirps_local_folder)) {
  saved_version <- if (file.exists(saved_version_file)) readLines(saved_version_file, warn = FALSE)[1] else
    "an unrecorded version"
  if (!any(str_detect(saved_version, fixed(paste("CHIRPS", chirps_versions))))) {
    stop("The saved CHIRPS values come from ", saved_version, ", not from CHIRPS ",
         paste(chirps_versions, collapse = " or "), ". Delete ", basename(cache_path("chirps_monthly")),
         " and ", basename(saved_version_file), " in ", environment_results_folder,
         " and run this script again.", call. = FALSE)
  }
}

if (nrow(months_to_do) == 0) {
  message("  CHIRPS precipitation is complete; no download is needed")
  chirps_values <- saved_chirps
  version_file <- sub("\\.csv$", "_source.txt", cache_path("chirps_monthly"))
  chirps_version_used <- if (file.exists(version_file)) readLines(version_file, warn = FALSE)[1] else
    "CHIRPS (version not recorded by the earlier run)"
  message("  source of the saved values: ", chirps_version_used)
  record_source("precipitation", paste0(chirps_version_used, " (values saved by an earlier run)"), period_text)
} else {
  message("  months still to process: ", nrow(months_to_do))
  if (given(chirps_local_folder)) {
    files <- list.files(chirps_local_folder, pattern = "\\.tif(\\.gz)?$", recursive = TRUE, full.names = TRUE)
    chirps_files <- tibble(source = files, name = basename(files))
    chirps_version_used <- paste("CHIRPS (local copies in", chirps_local_folder, ")")
  } else {
    chirps_files <- tibble()
    server_reached <- FALSE
    for (version in chirps_versions) {
      # The monthly folder differs between CHIRPS v3.0 and v2.0.
      start <- if (version == "3.0") "https://data.chc.ucsb.edu/products/CHIRPS/v3.0/monthly/" else
        "https://data.chc.ucsb.edu/products/CHIRPS-2.0/africa_monthly/"
      folder <- start
      for (step in 1:3) {
        links <- list_links(folder)
        if (is.null(links)) { links <- character(); break }
        server_reached <- TRUE
        if (any(grepl("\\.\\d{4}\\.\\d{2}\\.tif(\\.gz)?$", links))) break
        next_folder <- c(grep("^africa.*/$", links, value = TRUE), grep("^tifs/$", links, value = TRUE))
        if (length(next_folder) == 0) { links <- character(); break }
        folder <- paste0(folder, next_folder[1])
      }
      tif_links <- grep("\\.\\d{4}\\.\\d{2}\\.tif(\\.gz)?$", links, value = TRUE)
      if (length(tif_links) > 0) {
        chirps_files <- tibble(source = paste0(folder, basename(tif_links)), name = basename(tif_links))
        chirps_version_used <- paste("CHIRPS", version, "Africa monthly")
        break
      }
      message("  CHIRPS ", version, " monthly files were not found",
              if (version != tail(chirps_versions, 1)) "; trying the next version" else "")
    }
    if (nrow(chirps_files) == 0 && !server_reached) {
      stop("R could not reach data.chc.ucsb.edu (see the messages above). Check the internet connection ",
           "or proxy, or download the Africa monthly GeoTIFFs by hand and set chirps_local_folder ",
           "in 1_Harmonise_tree_data.R.", call. = FALSE)
    }
    if (nrow(chirps_files) == 0) {
      stop("The CHIRPS server was reached, but no monthly GeoTIFFs were listed. The folder layout may ",
           "have changed; download the files by hand and set chirps_local_folder in 1_Harmonise_tree_data.R.",
           call. = FALSE)
    }
  }
  message("  source: ", chirps_version_used)
  
  chirps_files <- chirps_files %>%
    mutate(year = as.integer(str_match(name, "\\.(\\d{4})\\.(\\d{2})\\.tif")[, 2]),
           month = as.integer(str_match(name, "\\.(\\d{4})\\.(\\d{2})\\.tif")[, 3])) %>%
    semi_join(months_to_do, by = c("year", "month")) %>%
    distinct(year, month, .keep_all = TRUE) %>%
    arrange(year, month)
  if (nrow(chirps_files) < nrow(months_to_do)) {
    warning("CHIRPS: only ", nrow(chirps_files), " of the ", nrow(months_to_do), " months still needed were ",
            "found; years with missing months are left out of mean annual precipitation.", call. = FALSE)
  }
  
  chirps_folder <- file.path(environment_data_folder, "CHIRPS")
  downloaded <- !given(chirps_local_folder)
  new_values <- list()
  # save_progress() adds the months processed so far to the cache.
  save_progress <- function() {
    if (length(new_values) == 0) return(invisible(NULL))
    finished <- bind_rows(new_values)
    kept <- if (is.null(saved_chirps)) NULL else
      anti_join(saved_chirps, distinct(finished, year, month), by = c("year", "month"))
    saved_chirps <<- bind_rows(kept, finished) %>% arrange(year, month)
    write_cache("chirps_monthly", saved_chirps)
    writeLines(chirps_version_used, sub("\\.csv$", "_source.txt", cache_path("chirps_monthly")))
    new_values <<- list()
  }
  
  for (i in seq_len(nrow(chirps_files))) {
    path <- chirps_files$source[i]
    if (grepl("^https?://", path)) {
      local_path <- file.path(chirps_folder, chirps_files$name[i])
      unpacked_path <- sub("\\.gz$", "", local_path)
      if (readable_raster(unpacked_path)) {
        path <- unpacked_path
      } else {
        # A GeoTIFF that cannot be opened is removed and downloaded again.
        if (!grepl("\\.gz$", local_path) && file.exists(local_path)) file.remove(local_path)
        if (!download_file(path, local_path)) stop("Could not download ", path, call. = FALSE)
        path <- local_path
      }
    }
    if (grepl("\\.gz$", path)) {
      path <- ungzip(path, chirps_folder, remove_original = downloaded)
    }
    values <- values_at_locations(terra::rast(path), locations$longitude, locations$latitude)[, 1]
    values[values < 0] <- NA  # -9999 marks missing data.
    new_values[[length(new_values) + 1]] <- tibble(latitude = locations$latitude_text,
                                                   longitude = locations$longitude_text,
                                                   year = chirps_files$year[i],
                                                   month = chirps_files$month[i],
                                                   precipitation_mm = values)
    # Local copies are never deleted.
    if (!keep_downloaded_files && downloaded) file.remove(path)
    if (chirps_files$month[i] == 12 || i == nrow(chirps_files)) {
      save_progress()
      message("  ", chirps_files$year[i], " done and saved")
    }
  }
  chirps_values <- saved_chirps
  record_source("precipitation", chirps_version_used, period_text,
                "monthly totals; mean annual precipitation = mean of the annual totals of complete years")
}

chirps_values <- chirps_values %>% filter(paste(latitude, longitude) %in% location_keys)

# Mean annual precipitation is the mean of the complete annual totals.
map_by_location <- chirps_values %>%
  filter(year %in% years) %>%
  group_by(latitude, longitude, year) %>%
  summarise(months = sum(!is.na(precipitation_mm)), total = sum(precipitation_mm), .groups = "drop") %>%
  filter(months == 12) %>%
  group_by(latitude, longitude) %>%
  summarise(precipitation_mm = mean(total), .groups = "drop")


# 6. Temperature: ERA5-Land monthly 2 m air temperature

message("\nTemperature (ERA5-Land):")
era5_values <- read_cache("era5land_temperature", locations)

if (is.null(era5_values)) {
  if (given(era5land_temperature_local_source)) {
    era5_files <- if (dir.exists(era5land_temperature_local_source)) {
      list.files(era5land_temperature_local_source, pattern = "\\.nc$", full.names = TRUE)
    } else era5land_temperature_local_source
  } else {
    need_package("ecmwfr", "download ERA5-Land from the Copernicus Climate Data Store")
    library(ecmwfr)
    era5_folder <- file.path(environment_data_folder, "ERA5-Land")
    dir.create(era5_folder, recursive = TRUE, showWarnings = FALSE)
    # One request per decade keeps each file a manageable size.
    year_groups <- split(years, ceiling(seq_along(years) / 10))
    targets <- vapply(year_groups, function(group) {
      sprintf("era5land_2m_temperature_%d_%d.nc", min(group), max(group))
    }, character(1))
    # A file that cannot be opened is removed and requested again.
    already_there <- vapply(file.path(era5_folder, targets), readable_raster, logical(1))
    for (broken in file.path(era5_folder, targets)[!already_there & file.exists(file.path(era5_folder, targets))]) {
      message("  removing an incomplete file: ", basename(broken))
      file.remove(broken)
    }
    message("  ERA5-Land files already downloaded: ", sum(already_there), " of ", length(targets))
    
    if (!all(already_there)) {
      ensure_cds_key()
      area <- c(ceiling((max(locations$latitude) + era5land_margin_degrees) * 10) / 10,
                floor((min(locations$longitude) - era5land_margin_degrees) * 10) / 10,
                floor((min(locations$latitude) - era5land_margin_degrees) * 10) / 10,
                ceiling((max(locations$longitude) + era5land_margin_degrees) * 10) / 10)
      for (k in which(!already_there)) {
        group <- year_groups[[k]]
        target <- targets[k]
        message("  requesting ", target, " from the Climate Data Store (this can take a while)")
        request <- list(
          dataset_short_name = "reanalysis-era5-land-monthly-means",
          product_type       = "monthly_averaged_reanalysis",
          variable           = "2m_temperature",
          year               = as.character(group),
          month              = sprintf("%02d", 1:12),
          time               = "00:00",
          data_format        = "netcdf",
          download_format    = "unarchived",
          area               = area,
          target             = target
        )
        tryCatch(
          ecmwfr::wf_request(request = request, user = cds_user, path = era5_folder, transfer = TRUE,
                             time_out = 6 * 3600),
          error = function(problem) {
            stop("The ERA5-Land request for ", target, " failed:\n  ", conditionMessage(problem), "\n",
                 "Check that (1) your API key is correct (run ecmwfr::wf_set_key() again if not), and ",
                 "(2) you have accepted the ERA5-Land licence under 'Download data' on ",
                 "https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land-monthly-means. ",
                 "Files already downloaded are kept; run the script again to continue.", call. = FALSE)
          })
        if (!readable_raster(file.path(era5_folder, target))) {
          stop("The Climate Data Store returned ", target, " but it cannot be opened as a raster. ",
               "Delete it and run the script again.", call. = FALSE)
        }
        message("  ", target, " downloaded")
      }
    }
    era5_files <- file.path(era5_folder, targets)
  }
  
  era5_values <- bind_rows(lapply(era5_files, function(path) {
    raster <- terra::rast(path)
    dates <- layer_dates(raster, path)
    values <- values_at_locations(raster, locations$longitude, locations$latitude)
    tibble(latitude = rep(locations$latitude_text, ncol(values)),
           longitude = rep(locations$longitude_text, ncol(values)),
           year = rep(as.integer(format(dates, "%Y")), each = nrow(locations)),
           month = rep(as.integer(format(dates, "%m")), each = nrow(locations)),
           temperature_k = as.vector(values))
  }))
  write_cache("era5land_temperature", era5_values)
  record_source("temperature", "ERA5-Land monthly averaged reanalysis, 2 m temperature", period_text,
                "mean annual temperature = mean of the monthly means weighted by the number of days")
} else {
  record_source("temperature", "ERA5-Land (values saved by an earlier run)", period_text)
}

# Monthly means are weighted by days in the month and converted to degrees C.
mat_by_location <- era5_values %>%
  filter(year %in% years) %>%
  mutate(days = days_in_month(year, month)) %>%
  group_by(latitude, longitude) %>%
  summarise(temperature_c = weighted.mean(temperature_k, days, na.rm = TRUE) - 273.15,
            months = sum(!is.na(temperature_k)), .groups = "drop")
if (any(mat_by_location$months < 12 * length(years))) {
  warning("ERA5-Land: some locations lack months of ", period_text, ".", call. = FALSE)
}


# 7. Potential evapotranspiration: TerraClimate 

message("\nPotential evapotranspiration (TerraClimate):")

terraclimate_servers <- c("https://www.reacchpna.org/thredds",
                          "http://thredds.northwestknowledge.net:8080/thredds")
pet_folder <- file.path(environment_data_folder, "TerraClimate")
dir.create(pet_folder, recursive = TRUE, showWarnings = FALSE)

# The 30-year climatology is used where it matches the climate period.
climatology_files <- c("1961-1990" = "TerraClimate19611990_pet.nc",
                       "1981-2010" = "TerraClimate19812010_pet.nc",
                       "1991-2020" = "TerraClimate19912020_pet.nc")
# Yearly values saved by an earlier run are used first.
saved_yearly_pet <- read_saved_values("terraclimate_pet")
yearly_pet_complete <- !is.null(saved_yearly_pet) && all(years %in% (saved_yearly_pet %>%
                                                                       filter(paste(latitude, longitude) %in% location_keys, !is.na(pet_mm)) %>%
                                                                       distinct(year, latitude, longitude) %>% count(year) %>% filter(n == length(location_keys)) %>% pull(year)))
use_climatology <- period_text %in% names(climatology_files) && !yearly_pet_complete &&
  identical(as.integer(years), min(years):max(years)) && !given(terraclimate_pet_local_source)
pet_by_location <- NULL

if (use_climatology) {
  climatology_name  <- climatology_files[[period_text]]
  saved_climatology <- read_cache("terraclimate_pet_climatology", locations)
  if (!is.null(saved_climatology)) {
    pet_by_location <- saved_climatology
    record_source("potential evapotranspiration",
                  paste("TerraClimate", period_text, "monthly climatology (values saved by an earlier run)"),
                  period_text)
  } else {
    climatology_path <- file.path(pet_folder, climatology_name)
    if (!readable_raster(climatology_path)) {
      if (file.exists(climatology_path)) file.remove(climatology_path)  # An incomplete earlier download.
      for (server in terraclimate_servers) {
        url <- sprintf("%s/fileServer/TERRACLIMATE_ALL/summaries/%s", server, climatology_name)
        message("  downloading ", climatology_name, " (about 120 MB, once only) from ", server)
        if (download_file(url, climatology_path) && readable_raster(climatology_path)) break
        if (file.exists(climatology_path)) file.remove(climatology_path)
        message("  not available from ", server)
      }
    }
    if (readable_raster(climatology_path)) {
      message("  extracting PET at ", nrow(locations), " locations")
      raster <- terra::rast(climatology_path)
      variables <- tryCatch(names(terra::sds(climatology_path)), error = function(e) character())
      if (length(variables) > 1 && "pet" %in% variables) raster <- terra::rast(climatology_path, subds = "pet")
      if (terra::nlyr(raster) != 12) stop(climatology_name, " does not hold 12 monthly layers.", call. = FALSE)
      monthly <- values_at_locations(raster, locations$longitude, locations$latitude)
      # Annual PET is the sum of the 12 monthly means.
      pet_by_location <- tibble(latitude = locations$latitude_text, longitude = locations$longitude_text,
                                pet_mm = rowSums(monthly))
      write_cache("terraclimate_pet_climatology", pet_by_location)
      record_source("potential evapotranspiration",
                    paste("TerraClimate", period_text, "monthly climatology,", climatology_name), period_text,
                    "annual PET = sum of the 12 monthly means of the climate period")
      message("  TerraClimate PET done and saved")
    } else {
      message("  the climatology file could not be downloaded; the yearly files are used instead")
    }
  }
}

# Otherwise one TerraClimate file is used per year.
if (is.null(pet_by_location)) {
  saved_pet <- read_saved_values("terraclimate_pet")
  if (!is.null(saved_pet)) {
    saved_pet <- saved_pet %>% filter(year %in% years)
    complete_years <- saved_pet %>%
      filter(paste(latitude, longitude) %in% location_keys, !is.na(pet_mm)) %>%
      distinct(year, latitude, longitude) %>%
      count(year, name = "locations_saved") %>%
      filter(locations_saved == length(location_keys)) %>%
      pull(year)
  } else {
    complete_years <- integer()
  }
  years_to_do <- setdiff(years, complete_years)
  message("  years of ", period_text, " already saved for all ", length(location_keys), " locations: ",
          length(years) - length(years_to_do), " of ", length(years))
  
  extract_pet <- function(path, year) {
    raster <- terra::rast(path)
    variables <- tryCatch(names(terra::sds(path)), error = function(e) character())
    if (length(variables) > 1 && "pet" %in% variables) raster <- terra::rast(path, subds = "pet")
    values <- values_at_locations(raster, locations$longitude, locations$latitude)
    if (ncol(values) != 12) stop(basename(path), " does not hold 12 monthly layers.", call. = FALSE)
    tibble(latitude = locations$latitude_text, longitude = locations$longitude_text, year = year,
           pet_mm = rowSums(values))
  }
  
  if (length(years_to_do) == 0) {
    message("  TerraClimate PET is complete; no download is needed")
    pet_values <- saved_pet
    record_source("potential evapotranspiration", "TerraClimate (values saved by an earlier run)", period_text)
  } else {
    message("  years still to process: ", length(years_to_do))
    box <- c(north = ceiling(max(locations$latitude)) + 1, south = floor(min(locations$latitude)) - 1,
             west = floor(min(locations$longitude)) - 1, east = ceiling(max(locations$longitude)) + 1)
    pet_sources <- character()
    
    for (year in years_to_do) {
      if (given(terraclimate_pet_local_source)) {
        files <- if (dir.exists(terraclimate_pet_local_source)) {
          list.files(terraclimate_pet_local_source, pattern = "\\.nc$", full.names = TRUE)
        } else terraclimate_pet_local_source
        path <- files[grepl(paste0("_", year, "(_|\\.)"), basename(files))][1]
        if (is.na(path)) stop("No TerraClimate PET file for ", year, " in ", terraclimate_pet_local_source,
                              call. = FALSE)
        source_used <- "TerraClimate yearly files (local copies)"
      } else {
        subset_path <- file.path(pet_folder, sprintf("TerraClimate_pet_%d_baobab_area.nc", year))
        global_path <- file.path(pet_folder, sprintf("TerraClimate_pet_%d.nc", year))
        if (readable_raster(subset_path)) {
          path <- subset_path; source_used <- "TerraClimate PET, NetCDF Subset Service (file already on disk)"
        } else if (readable_raster(global_path)) {
          path <- global_path; source_used <- "TerraClimate PET yearly files (file already on disk)"
        } else {
          path <- NA_character_
          message("  downloading TerraClimate PET for ", year)
          # Both forms of the THREDDS subset address are tried.
          subset_query <- sprintf(paste0("agg_terraclimate_pet_1958_CurrentYear_GLOBE.nc?var=pet",
                                         "&north=%s&west=%s&east=%s&south=%s&horizStride=1",
                                         "&time_start=%d-01-01T00:00:00Z&time_end=%d-12-31T23:59:59Z",
                                         "&accept=netcdf"),
                                  box[["north"]], box[["west"]], box[["east"]], box[["south"]], year, year)
          subset_urls <- as.vector(outer(terraclimate_servers, c("/ncss/grid/", "/ncss/"), paste0))
          for (url in paste0(subset_urls, subset_query)) {
            if (file.exists(subset_path)) file.remove(subset_path)
            if (download_file(url, subset_path, attempts = 2) && readable_raster(subset_path) &&
                terra::nlyr(terra::rast(subset_path)) >= 12) {
              path <- subset_path
              source_used <- paste("TerraClimate PET, NetCDF Subset Service,", sub("/ncss.*$", "", url))
              break
            }
          }
          # Fallback: the global file of the year.
          if (is.na(path)) {
            for (server in terraclimate_servers) {
              if (file.exists(global_path)) file.remove(global_path)
              url <- sprintf("%s/fileServer/TERRACLIMATE_ALL/data/TerraClimate_pet_%d.nc", server, year)
              if (download_file(url, global_path) && readable_raster(global_path)) {
                path <- global_path; source_used <- paste("TerraClimate PET yearly files,", server)
                break
              }
            }
          }
          if (is.na(path)) {
            stop("TerraClimate PET for ", year, " could not be downloaded from either server:\n  ",
                 paste(terraclimate_servers, collapse = "\n  "), "\n",
                 "Check your internet connection (some networks block these servers). You can also ",
                 "download the yearly files TerraClimate_pet_YYYY.nc yourself from ",
                 "https://www.climatologylab.org/terraclimate.html, put them in one folder and give that ",
                 "folder as terraclimate_pet_local_source in 1_Harmonise_tree_data.R (Section 1). Years already saved are kept.",
                 call. = FALSE)
          }
        }
      }
      year_values <- extract_pet(path, year)
      # Each year is cached at once, so an interrupted run resumes.
      saved_pet <- bind_rows(if (!is.null(saved_pet)) filter(saved_pet, year != !!year), year_values) %>%
        arrange(year)
      write_cache("terraclimate_pet", saved_pet)
      if (!keep_downloaded_files && !given(terraclimate_pet_local_source)) file.remove(path)
      pet_sources <- union(pet_sources, source_used)
      message("  ", year, " done and saved")
    }
    pet_values <- saved_pet
    record_source("potential evapotranspiration", paste(pet_sources, collapse = "; "), period_text,
                  "annual totals of monthly PET")
  }
  pet_values <- pet_values %>% filter(paste(latitude, longitude) %in% location_keys)
  
  pet_by_location <- pet_values %>%
    group_by(latitude, longitude) %>%
    summarise(pet_mm = mean(pet_mm, na.rm = TRUE), .groups = "drop")
}
implausible_pet <- sum(pet_by_location$pet_mm < 300 | pet_by_location$pet_mm > 4000, na.rm = TRUE)
if (implausible_pet > 0) {
  warning(implausible_pet, " locations have an annual PET outside 300-4,000 mm; check the scale factor ",
          "of the TerraClimate files.", call. = FALSE)
}


# 8. Soil texture: iSDAsoil sand, silt and clay (0-20 cm) 

message("\nSoil texture (iSDAsoil):")
soil_values <- read_cache("isdasoil_texture", locations)

if (is.null(soil_values)) {
  terra::setGDALconfig("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
  terra::setGDALconfig("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", ".tif")
  raw_soil <- list(); transforms <- list()
  for (fraction in names(isda_layers)) {
    layer <- isda_layers[[fraction]]
    source <- if (given(isda_local_folder)) file.path(isda_local_folder, paste0(layer, ".tif")) else
      sprintf("/vsicurl/%s/%s/%s.tif", isda_base_url, layer, layer)
    raster <- tryCatch(terra::rast(source), error = function(e) {
      stop("Could not open the iSDAsoil layer ", source, ". Check the address on the iSDAsoil ",
           "STAC browser (https://isdasoil.s3.amazonaws.com/index.html).", call. = FALSE)
    })
    band <- which(grepl("mean_0_20", names(raster), ignore.case = TRUE))[1]
    if (is.na(band)) band <- 1
    raw_soil[[fraction]] <- values_at_locations(raster[[band]], locations$longitude, locations$latitude)[, 1]
    transforms[[fraction]] <- if (isda_back_transformation == "metadata") {
      metadata_transformation(sprintf("%s/%s/%s.json", isda_base_url, layer, layer))
    } else as_transformation(isda_back_transformation)
    if (is.null(transforms[[fraction]])) {
      stop("The iSDAsoil metadata file for ", layer, " could not be read, so the back-transformation ",
           "is unknown. Connect to the internet, or set isda_back_transformation in Section 1 after ",
           "checking the iSDAsoil documentation.", call. = FALSE)
    }
    message("  ", layer, ": band ", band, "; back-transformation: ", attr(transforms[[fraction]], "text"))
  }
  chosen <- lapply(names(raw_soil), function(f) transforms[[f]](raw_soil[[f]]))
  names(chosen) <- names(raw_soil)
  # A fraction above 100% means a wrong back-transformation, so the script
  # stops.
  above_100 <- sum(unlist(chosen) > 100.5, na.rm = TRUE)
  if (above_100 > 0) {
    stop("iSDAsoil: ", above_100, " sand, silt or clay values exceed 100% after the back-transformation ",
         attr(transforms$sand, "text"), ". The back-transformation does not fit these layers; check the ",
         "iSDAsoil documentation and set isda_back_transformation in Section 1.", call. = FALSE)
  }
  sum_difference <- median(abs(chosen$sand + chosen$silt + chosen$clay - 100), na.rm = TRUE)
  soil_note <- paste0("back-transformation ", attr(transforms$sand, "text"),
                      "; median |sand + silt + clay - 100| = ", round(sum_difference, 1), "%")
  message("  ", soil_note)
  if (sum_difference > 5) {
    warning("iSDAsoil sand, silt and clay do not add up to about 100% (median difference ",
            round(sum_difference, 1), "%). Check the back-transformation.", call. = FALSE)
  }
  soil_values <- tibble(latitude = locations$latitude_text, longitude = locations$longitude_text,
                        sand_percent = chosen$sand, silt_percent = chosen$silt, clay_percent = chosen$clay)
  write_cache("isdasoil_texture", soil_values)
  record_source("soil texture", "iSDAsoil sand, silt and clay content, predicted mean 0-20 cm", "v1", soil_note)
} else {
  record_source("soil texture", "iSDAsoil (values saved by an earlier run)", "v1")
}
soil_values <- soil_values %>%
  mutate(texture_class = usda_texture_class(sand_percent, silt_percent, clay_percent))


# 9. Protected areas: WDPA cleaned with wdpar 

message("\nProtected areas (WDPA):")
sf_use_s2(FALSE)  # Planar tests are robust to imperfect polygons.
protected_values <- read_cache("protected_areas", locations)

if (is.null(protected_values)) {
  # bind_country_tables() gives every table the same columns before binding.
  bind_country_tables <- function(tables) {
    target_crs <- sf::st_crs(tables[[1]])
    tables <- lapply(tables, function(table) {
      sf::st_geometry(table) <- "geometry"
      if (sf::st_crs(table) != target_crs) table <- sf::st_transform(table, target_crs)
      table
    })
    all_columns <- unique(unlist(lapply(tables, function(table) setdiff(names(table), "geometry"))))
    tables <- lapply(tables, function(table) {
      for (column in setdiff(all_columns, names(table))) table[[column]] <- NA
      table[, c(all_columns, "geometry")]
    })
    do.call(rbind, tables)
  }
  
  if (given(protected_areas_local_source)) {
    protected_areas <- sf::st_read(protected_areas_local_source, quiet = TRUE)
    country_tables <- list(local_copy = protected_areas)
    wdpa_version <- paste("local copy:", basename(protected_areas_local_source))
  } else {
    need_package("wdpar", "download the World Database on Protected Areas")
    library(wdpar)
    countries <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
    countries$iso3 <- ifelse(countries$iso_a3 == "-99", countries$adm0_a3, countries$iso_a3)
    location_points <- st_as_sf(locations, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
    sf_use_s2(FALSE)
    iso3 <- unique(na.omit(countries$iso3[unlist(lapply(
      suppressMessages(st_intersects(location_points, countries)), function(i) if (length(i)) i[1] else NA))]))
    wdpa_folder <- file.path(environment_data_folder, "WDPA")
    dir.create(wdpa_folder, recursive = TRUE, showWarnings = FALSE)
    
    # local_wdpa_file() returns the latest WDPA file of a country already on
    # disk.
    local_wdpa_file <- function(code) {
      files <- list.files(wdpa_folder, pattern = paste0("^WDPA_.*_", code, "-shapefile\\.zip$"),
                          full.names = TRUE)
      if (length(files) == 0) return(NA_character_)
      files[which.max(file.mtime(files))]
    }
    
    # Each country is read from disk, or downloaded only if absent.
    country_tables <- lapply(iso3, function(code) {
      saved_file <- local_wdpa_file(code)
      if (!is.na(saved_file)) {
        message("  reading WDPA for ", code, " (already downloaded: ", basename(saved_file), ")")
        wdpar::wdpa_read(saved_file)
      } else {
        message("  downloading WDPA for ", code)
        wdpar::wdpa_fetch(code, wait = TRUE, download_dir = wdpa_folder, check_version = TRUE)
      }
    })
    
    names(country_tables) <- iso3
    message("  protected areas read: ", sum(vapply(country_tables, nrow, integer(1))), " areas in ",
            length(iso3), " countries")
    
    releases <- unique(na.omit(str_match(basename(vapply(iso3, local_wdpa_file, character(1))),
                                         "^WDPA_([A-Za-z]{3}\\d{4})_")[, 2]))
    wdpa_version <- if (length(releases) > 0) paste("WDPA release", paste(releases, collapse = ", ")) else
      "WDPA (current release on the access date)"
  }
  
  cleaning_crs    <- "ESRI:54017"
  retained_status <- c("Designated", "Inscribed", "Established")
  
  # repair_one_by_one() repairs each geometry at full precision, so one failure
  # cannot stop the others.
  repair_one_by_one <- function(table) {
    table <- sf::st_set_precision(table, 0)
    geometries <- sf::st_geometry(table)
    repaired <- vector("list", length(geometries))
    for (i in seq_along(geometries)) {
      one <- geometries[i]
      fixed <- tryCatch(sf::st_make_valid(one), error = function(e) NULL)
      if (is.null(fixed)) fixed <- tryCatch(sf::st_buffer(one, 0), error = function(e) NULL)
      repaired[[i]] <- if (is.null(fixed) || all(sf::st_is_empty(fixed))) NULL else fixed[[1]]
    }
    keep <- !vapply(repaired, is.null, logical(1))
    if (any(!keep)) message("    geometries that could not be repaired and were left out: ", sum(!keep))
    table <- table[keep, ]
    sf::st_geometry(table) <- sf::st_sfc(repaired[keep], crs = sf::st_crs(geometries))
    is_collection <- sf::st_geometry_type(table) == "GEOMETRYCOLLECTION"
    if (any(is_collection)) {
      polygons <- suppressWarnings(sf::st_collection_extract(table[is_collection, ], "POLYGON"))
      table <- rbind(table[!is_collection, ], polygons)
    }
    table
  }
  
  # clean_with_same_rules() applies the wdpa_clean() rules without grid
  # rounding.
  clean_with_same_rules <- function(table) {
    table <- table[table$STATUS %in% retained_status, ]
    table <- table[!(table$DESIG_ENG %in% "UNESCO-MAB Biosphere Reserve"), ]
    is_point <- sf::st_geometry_type(table) %in% c("POINT", "MULTIPOINT")
    table <- table[!(is_point & !is.finite(table$REP_AREA)), ]
    table <- repair_one_by_one(table)
    table <- sf::st_transform(table, cleaning_crs)
    is_point <- sf::st_geometry_type(table) %in% c("POINT", "MULTIPOINT")
    if (any(is_point)) {
      circles <- sf::st_buffer(table[is_point, ], sqrt((table$REP_AREA[is_point] * 1e6) / pi))
      table <- rbind(table[!is_point, ], circles)
    }
    table
  }
  
  fallback_countries <- character()
  if (requireNamespace("wdpar", quietly = TRUE)) {
    library(wdpar)
    cleaned_tables <- lapply(names(country_tables), function(code) {
      message("  cleaning the protected areas of ", code, " with wdpar::wdpa_clean()")
      cleaned <- tryCatch(
        wdpar::wdpa_clean(country_tables[[code]], erase_overlaps = erase_protected_area_overlaps,
                          verbose = FALSE),
        error = function(problem) {
          message("    wdpa_clean() stopped for ", code, ": ", conditionMessage(problem))
          NULL
        })
      if (is.null(cleaned)) {
        message("    the same selection rules are applied to ", code, " without grid rounding")
        fallback_countries <<- c(fallback_countries, code)
        cleaned <- clean_with_same_rules(country_tables[[code]])
      }
      message("    ", nrow(cleaned), " protected areas kept for ", code)
      cleaned
    })
    protected_areas <- bind_country_tables(cleaned_tables)
    message("  cleaning finished: ", nrow(protected_areas), " protected areas in ", length(cleaned_tables),
            if (length(cleaned_tables) == 1) " table" else " countries")
  } else {
    warning("wdpar is not installed: the same selection rules are applied without wdpa_clean().", call. = FALSE)
    fallback_countries <- names(country_tables)
    protected_areas <- bind_country_tables(lapply(country_tables, clean_with_same_rules))
  }
  
  cleaning_note <- paste0(
    "statuses Designated, Inscribed and Established; UNESCO Biosphere Reserves excluded; ",
    "points without a reported area removed and other points buffered to their reported area; ",
    if (length(fallback_countries) == 0) "all countries cleaned with wdpar::wdpa_clean()" else
      paste0("cleaned with wdpar::wdpa_clean() except ", paste(fallback_countries, collapse = ", "),
             ", for which wdpa_clean() failed on a geometry and the same selection rules were applied ",
             "with geometries repaired at full precision (sf::st_make_valid)"))
  if (length(fallback_countries) > 0) message("  ", cleaning_note)
  
  # Locations are compared in the equal-area system of the protected areas.
  location_points <- st_transform(st_as_sf(locations, coords = c("longitude", "latitude"), crs = 4326),
                                  st_crs(protected_areas))
  inside <- suppressMessages(st_intersects(location_points, protected_areas))
  first_area <- vapply(inside, function(i) if (length(i)) i[1] else NA_integer_, integer(1))
  name_column <- intersect(c("NAME", "NAME_ENG"), names(protected_areas))[1]
  protected_values <- tibble(latitude = locations$latitude_text, longitude = locations$longitude_text,
                             in_protected_area = lengths(inside) > 0,
                             protected_area_name = if (!is.na(name_column))
                               as.character(st_drop_geometry(protected_areas)[[name_column]][first_area]) else NA_character_)
  write_cache("protected_areas", protected_values)
  record_source("protected areas", "World Database on Protected Areas (UNEP-WCMC and IUCN), cleaned with wdpar",
                wdpa_version, cleaning_note)
} else {
  record_source("protected areas", "WDPA (values saved by an earlier run)", "see earlier run")
}


# 10. Slope: Copernicus DEM GLO-30 

message("\nSlope (Copernicus DEM GLO-30):")
slope_cache_name <- "slope_glo30"
slope_values <- read_cache(slope_cache_name, locations)
if (is.null(slope_values)) {
  half_window <- 0.005
  
  # slope_at_points() calculates slope at each location from a small window of the elevation model.
  slope_at_points <- function(dem, longitude, latitude, half_window = 0.005) {
    vapply(seq_along(longitude), function(i) {
      window <- terra::ext(longitude[i] - half_window, longitude[i] + half_window,
                           latitude[i] - half_window, latitude[i] + half_window)
      piece <- tryCatch(terra::crop(dem, window), error = function(e) NULL)
      if (is.null(piece)) return(NA_real_)
      slope <- terra::terrain(piece, v = "slope", unit = "degrees")
      as.numeric(terra::extract(slope, cbind(longitude[i], latitude[i]))[1, 1])
    }, numeric(1))
  }
  
  slope <- rep(NA_real_, nrow(locations))
  if (given(dem_local_source)) {
    dem <- terra::rast(dem_local_source)
    if (!terra::is.lonlat(dem)) stop("The elevation model must be in longitude and latitude (EPSG:4326).", call. = FALSE)
    slope <- slope_at_points(dem, locations$longitude, locations$latitude, half_window)
    dem_version <- paste("local copy:", basename(dem_local_source))
  } else {
    terra::setGDALconfig("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
    
    # tile_name() gives the GLO-30 tile that contains a point.
    tile_name <- function(lat, lon) {
      sprintf("Copernicus_DSM_COG_10_%s%02d_00_%s%03d_00_DEM",
              ifelse(floor(lat) < 0, "S", "N"), abs(floor(lat)),
              ifelse(floor(lon) < 0, "W", "E"), abs(floor(lon)))
    }
    # tile_address() gives the web address of a tile, read directly by GDAL.
    tile_address <- function(tile) {
      paste0("/vsicurl/https://copernicus-dem-30m.s3.amazonaws.com/", tile, "/", tile, ".tif")
    }
    
    corners <- expand.grid(x_side = c(-1, 1), y_side = c(-1, 1))
    tiles_needed <- unique(unlist(lapply(seq_len(nrow(corners)), function(k) {
      tile_name(locations$latitude + corners$y_side[k] * half_window,
                locations$longitude + corners$x_side[k] * half_window)
    })))
    
    tile_exists <- vapply(tiles_needed, function(tile) {
      !is.null(tryCatch(terra::rast(tile_address(tile)), error = function(e) NULL))
    }, logical(1))
    for (tile in tiles_needed[!tile_exists]) message("  tile not available (sea or no access): ", tile)
    
    locations_per_tile <- table(tile_name(locations$latitude, locations$longitude))
    for (tile in names(locations_per_tile)) message("  ", tile, ": ", locations_per_tile[[tile]], " locations")
    
    if (any(tile_exists)) {
      dem <- terra::vrt(tile_address(tiles_needed[tile_exists]))
      slope <- slope_at_points(dem, locations$longitude, locations$latitude, half_window)
    }
    dem_version <- "Copernicus DEM GLO-30 (ESA and Airbus), 2021 release, open AWS archive"
  }
  
  slope_values <- tibble(latitude = locations$latitude_text, longitude = locations$longitude_text,
                         slope_degrees = round(slope, 2))
  if (any(!is.na(slope))) {
    write_cache(slope_cache_name, slope_values)
    record_source("slope", "Slope from a digital elevation model, terra::terrain() over a ~1 km window",
                  dem_version, "degrees")
  } else {
    message("  no slope could be calculated; the slope panel will be empty. Set dem_local_source in 1_Harmonise_tree_data.R.")
  }
} else {
  record_source("slope", "Slope (values saved by an earlier run, Copernicus DEM GLO-30)", "see earlier run")
}


# 11. Combining and saving the covariates 

covariates <- locations %>%
  select(location_id, latitude = latitude_text, longitude = longitude_text, country, records) %>%
  left_join(map_by_location, by = c("latitude", "longitude")) %>%
  left_join(select(mat_by_location, latitude, longitude, temperature_c), by = c("latitude", "longitude")) %>%
  left_join(pet_by_location, by = c("latitude", "longitude")) %>%
  left_join(soil_values, by = c("latitude", "longitude")) %>%
  left_join(protected_values, by = c("latitude", "longitude")) %>%
  left_join(slope_values, by = c("latitude", "longitude")) %>%
  mutate(aridity_index = ifelse(pet_mm > 0, precipitation_mm / pet_mm, NA_real_),
         aridity_class = cut(aridity_index, c(aridity_classes$lower, Inf), labels = aridity_classes$class,
                             right = FALSE))

# The covariate file has one row per unique location, with the number of
# baobabs recorded there (records). Figure 3 and the summary describe the
# baobabs themselves, so each location's row is repeated once for each baobab
# recorded at it: baobabs that share a location share its values.

baobab_covariates <- uncount(covariates, records, .remove = FALSE)
message("Baobabs: ", nrow(baobab_covariates), " at ", nrow(covariates), " unique locations")

continuous <- c(precipitation_mm = "Mean annual precipitation (mm)", temperature_c = "Mean annual temperature (\u00b0C)",
                pet_mm = "Annual potential evapotranspiration (mm)", aridity_index = "Aridity index (MAP / PET)",
                slope_degrees = "Slope (degrees)",
                sand_percent = "Sand (%)", silt_percent = "Silt (%)", clay_percent = "Clay (%)")
# The summary is calculated over baobabs; the number of unique locations
# with a value is given beside the number of baobabs.
covariate_summary <- bind_rows(lapply(names(continuous), function(variable) {
  values <- baobab_covariates[[variable]]
  tibble(covariate = continuous[[variable]], baobabs_with_value = sum(!is.na(values)),
         locations_with_value = sum(!is.na(covariates[[variable]])),
         median = median(values, na.rm = TRUE),
         lower_quartile = quantile(values, 0.25, na.rm = TRUE, names = FALSE),
         upper_quartile = quantile(values, 0.75, na.rm = TRUE, names = FALSE))
}), tibble(covariate = "Inside a protected area (% of baobabs)",
           baobabs_with_value = sum(!is.na(baobab_covariates$in_protected_area)),
           locations_with_value = sum(!is.na(covariates$in_protected_area)),
           median = 100 * mean(baobab_covariates$in_protected_area, na.rm = TRUE))) %>%
  mutate(across(c(median, lower_quartile, upper_quartile), ~ round(.x, 3)))

write_excel_csv(covariates, covariates_output, na = "NA")
write_excel_csv(covariate_summary, covariate_summary_output, na = "NA")
write_excel_csv(sources_used, sources_output, na = "NA")


# 12. Figure 3: environment of the baobabs

text_size_pt <- 7
small_text_size_pt <- 6
line_width <- function(points) points / (ggplot2::.pt * 72 / 96)  # points to ggplot2 linewidth
figure_font <- intersect(c("Arial", "Helvetica", "Liberation Sans", "Arimo"),
                         unique(systemfonts::system_fonts()$family))[1]
if (is.na(figure_font)) figure_font <- "sans"

baobab_colour   <- "#2a78d6"   # baobab colour of Figure 1
texture_colours <- c(Sand = "#eb6834", Silt = "#1baf7a", Clay = "#2a78d6")

theme_publication <- function() {
  theme_classic(base_size = text_size_pt, base_family = figure_font,
                base_line_size = line_width(0.5), base_rect_size = line_width(0.5)) +
    theme(text = element_text(colour = "black"),
          axis.title = element_text(size = text_size_pt, colour = "black"),
          axis.text = element_text(size = small_text_size_pt, colour = "black"),
          axis.line = element_line(linewidth = line_width(0.5), colour = "black"),
          axis.ticks = element_line(linewidth = line_width(0.5), colour = "black"),
          axis.ticks.length = unit(1, "mm"),
          panel.grid = element_blank(),
          legend.text = element_text(size = small_text_size_pt, colour = "black"),
          legend.key.size = unit(3, "mm"),
          plot.tag = element_text(size = 8, face = "bold", colour = "black"),
          plot.subtitle = element_text(size = small_text_size_pt, colour = "black",
                                       margin = margin(b = 1.5, unit = "mm")),
          plot.margin = margin(2, 3, 2, 2, unit = "mm"))
}
format_number <- function(x) trimws(format(round(x), big.mark = ","))

# histogram_panel() draws the number of baobabs by one variable in
# 40 bins. The bins are set to run exactly from the lowest to the highest
# value, and the x-axis has no padding, so the axis ends where the first and
# last bins end.

histogram_panel <- function(variable, x_label, bins = 40) {
  values <- baobab_covariates[[variable]][!is.na(baobab_covariates[[variable]])]
  value_range <- range(values)
  bin_edges <- if (diff(value_range) > 0) seq(value_range[1], value_range[2], length.out = bins + 1) else NULL
  message("  Figure 3, ", x_label, ": axis from ", signif(value_range[1], 4), " to ",
          signif(value_range[2], 4), " (", format_number(length(values)), " baobabs)")
  ggplot(tibble(value = values), aes(x = value)) +
    geom_histogram(breaks = bin_edges, bins = if (is.null(bin_edges)) bins else NULL,
                   fill = baobab_colour, colour = "white", linewidth = line_width(0.25)) +
    scale_x_continuous(labels = scales::label_comma(), expand = expansion(mult = 0)) +
    scale_y_continuous(labels = scales::label_comma(), expand = expansion(mult = c(0, 0.05))) +
    labs(x = x_label, y = "Baobabs") +
    theme_publication()
}

panel_precipitation <- histogram_panel("precipitation_mm", "Mean annual precipitation (mm)")
panel_temperature   <- histogram_panel("temperature_c", "Mean annual temperature (\u00b0C)")
panel_aridity       <- histogram_panel("aridity_index", "Aridity index (MAP/PET)")
panel_slope         <- histogram_panel("slope_degrees", "Slope (\u00b0)")

protected_table <- bind_rows(
  baobab_covariates %>% filter(!is.na(in_protected_area)) %>%
    summarise(group = "All baobabs", n = n(), percent = 100 * mean(in_protected_area)),
  baobab_covariates %>% filter(!is.na(in_protected_area)) %>%
    group_by(group = coalesce(country, "Country not assigned")) %>%
    summarise(n = n(), percent = 100 * mean(in_protected_area), .groups = "drop") %>%
    filter(n >= minimum_baobabs_per_country) %>% arrange(desc(n))
) %>%
  mutate(category = paste0(group, " (n = ", format_number(n), ")"),
         category = factor(category, levels = rev(category)))

# Each bar is split into the shares inside and outside protected areas, which
# add up to 100%. "Inside" is the last level, so it is drawn from zero.
protected_long <- protected_table %>%
  transmute(category, Inside = percent, Outside = 100 - percent) %>%
  pivot_longer(c(Inside, Outside), names_to = "status", values_to = "percent") %>%
  mutate(status = factor(status, levels = c("Outside", "Inside")))
protected_colours <- c(Inside = baobab_colour, Outside = "grey85")

panel_protected <- ggplot(protected_long, aes(x = percent, y = category, fill = status)) +
  geom_col(width = 0.65, position = "stack") +
  scale_fill_manual(values = protected_colours, breaks = c("Inside", "Outside"),
                    labels = c("Inside a protected area", "Outside protected areas"), name = NULL) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Baobabs (%)", y = NULL,
       subtitle = if (minimum_baobabs_per_country > 1)
         paste0("Countries with \u2265", minimum_baobabs_per_country, " baobabs") else NULL) +
  theme_publication() +
  theme(axis.ticks.y = element_blank(), plot.title.position = "plot",
        legend.position = "top", legend.margin = margin(0, 0, 0, 0))

# Panel f: the USDA soil texture triangle. Each baobab is placed by the
# topsoil sand, silt and clay content of its location, over the twelve USDA
# texture classes.
# The triangle is drawn with ggplot2 alone, so that it combines with the
# other panels; no ternary-plot package is needed.

# texture_xy() converts sand, silt and clay (%) into plotting positions. The
# conventional orientation is used: 100% sand at the lower left, 100% silt at
# the lower right and 100% clay at the top. Clay is read on the left side,
# silt on the right side and sand along the base.

texture_xy <- function(clay, silt) {
  tibble(x = silt + clay / 2, y = clay * sqrt(3) / 2)
}

# Corners of the twelve USDA texture classes, as (clay, silt) in per cent;
# sand is 100 - clay - silt. These follow the class rules of
# usda_texture_class() above.

texture_class_corners <- tribble(
  ~class,            ~clay,  ~silt,
  "Sand",                0,      0,   "Sand",  10,  0,   "Sand",   0, 15,
  "Loamy sand",         10,      0,   "Loamy sand", 15, 0,  "Loamy sand", 0, 30,  "Loamy sand", 0, 15,
  "Sandy loam",         15,      0,   "Sandy loam", 20, 0,  "Sandy loam", 20, 28, "Sandy loam", 7, 41,
  "Sandy loam",          7,     50,   "Sandy loam",  0, 50, "Sandy loam",  0, 30,
  "Loam",                7,     41,   "Loam", 20, 28, "Loam", 27, 28, "Loam", 27, 50, "Loam", 7, 50,
  "Silt loam",           0,     50,   "Silt loam", 27, 50, "Silt loam", 27, 73, "Silt loam", 12, 88,
  "Silt loam",          12,     80,   "Silt loam",  0, 80,
  "Silt",                0,     80,   "Silt", 12, 80, "Silt", 12, 88, "Silt", 0, 100,
  "Sandy clay loam",    20,      0,   "Sandy clay loam", 35, 0, "Sandy clay loam", 35, 20,
  "Sandy clay loam",    27,     28,   "Sandy clay loam", 20, 28,
  "Clay loam",          27,     28,   "Clay loam", 40, 15, "Clay loam", 40, 40, "Clay loam", 27, 53,
  "Silty clay loam",    27,     53,   "Silty clay loam", 40, 40, "Silty clay loam", 40, 60,
  "Silty clay loam",    27,     73,
  "Sandy clay",         35,      0,   "Sandy clay", 55, 0, "Sandy clay", 35, 20,
  "Silty clay",         40,     40,   "Silty clay", 60, 40, "Silty clay", 40, 60,
  "Clay",               40,     15,   "Clay", 55, 0, "Clay", 100, 0, "Clay", 60, 40, "Clay", 40, 40
)
texture_class_outlines <- bind_cols(texture_class_corners, texture_xy(texture_class_corners$clay,
                                                                      texture_class_corners$silt))

# Label positions (clay, silt) inside each class; the loamy sand label is
# turned to follow its narrow band. The check below stops the script if a
# label position falls outside its class.

texture_class_labels <- tribble(
  ~class,            ~text,                ~clay, ~silt, ~angle,
  "Sand",            "Sand",                   3,     4,      0,
  "Loamy sand",      "Loamy sand",           6.5,    13,    -30,
  "Sandy loam",      "Sandy loam",            10,    24,      0,
  "Loam",            "Loam",                  17,    40,      0,
  "Silt loam",       "Silt loam",             14,    65,      0,
  "Silt",            "Silt",                   5,    88,      0,
  "Sandy clay loam", "Sandy clay loam",       27,    12,      0,
  "Clay loam",       "Clay loam",             33,    34,      0,
  "Silty clay loam", "Silty clay\nloam",      33,    56,      0,
  "Sandy clay",      "Sandy\nclay",           40,   7.5,      0,
  "Silty clay",      "Silty\nclay",           46,    46,      0,
  "Clay",            "Clay",                  62,    18,      0
)
stopifnot(identical(usda_texture_class(100 - texture_class_labels$clay - texture_class_labels$silt,
                                       texture_class_labels$silt, texture_class_labels$clay),
                    texture_class_labels$class))

# Share of baobabs in each class, saved for every class and added to the
# label of the larger classes.
location_texture <- baobab_covariates %>%
  filter(!is.na(sand_percent), !is.na(silt_percent), !is.na(clay_percent)) %>%
  mutate(total = sand_percent + silt_percent + clay_percent,
         clay = 100 * clay_percent / total, silt = 100 * silt_percent / total,
         class = usda_texture_class(sand_percent, silt_percent, clay_percent))
class_share <- location_texture %>%
  count(class, name = "baobabs") %>%
  right_join(select(texture_class_labels, class), by = "class") %>%
  mutate(baobabs = coalesce(baobabs, 0L),
         percent_of_baobabs = round(100 * baobabs / sum(baobabs), 1)) %>%
  arrange(desc(baobabs))
write_excel_csv(class_share, file.path(environment_results_folder, "texture_class_shares.csv"), na = "NA")
texture_class_labels <- texture_class_labels %>%
  left_join(class_share, by = "class") %>%
  mutate(label = ifelse(percent_of_baobabs >= texture_share_label_minimum,
                        paste0(text, "\n", round(percent_of_baobabs), "%"), text)) %>%
  bind_cols(texture_xy(.$clay, .$silt))
location_points <- bind_cols(location_texture, texture_xy(location_texture$clay, location_texture$silt))


haloed_text <- function(data, size, halo = 0.45, ...) {
  shifts <- tibble(angle = seq(0, 315, by = 45)) %>%
    mutate(dx = halo * cos(angle * pi / 180), dy = halo * sin(angle * pi / 180))
  halo_layers <- lapply(seq_len(nrow(shifts)), function(i) {
    geom_text(data = mutate(data, x = x + shifts$dx[i], y = y + shifts$dy[i]),
              aes(x = x, y = y, label = label, angle = angle), size = size, colour = "white", ...)
  })
  c(halo_layers, list(geom_text(data = data, aes(x = x, y = y, label = label, angle = angle),
                                size = size, colour = "grey10", ...)))
}

# Ticks every 10% on each side, drawn outwards along the lines of equal
# content, as on the conventional USDA triangle.

tick_values  <- seq(0, 100, by = 10)
tick_length  <- 2.5
label_offset <- 5.5
side_ticks <- bind_rows(
  # Base: sand, read from right (0%) to left (100%).
  tibble(axis = "sand", value = tick_values, texture_xy(0, 100 - tick_values),
         dx = 0.5, dy = -sqrt(3) / 2, hjust = 0.5, vjust = 1),
  # Left side: clay, read upwards.
  tibble(axis = "clay", value = tick_values, texture_xy(tick_values, 0),
         dx = -1, dy = 0, hjust = 1, vjust = 0.5),
  # Right side: silt, read downwards from the top.
  tibble(axis = "silt", value = tick_values, texture_xy(100 - tick_values, tick_values),
         dx = 0.5, dy = sqrt(3) / 2, hjust = 0, vjust = 0.5)
) %>%
  mutate(xend = x + tick_length * dx, yend = y + tick_length * dy,
         label_x = x + label_offset * dx, label_y = y + label_offset * dy)
triangle_outline <- texture_xy(c(0, 0, 100, 0), c(0, 100, 0, 0))

side_midpoints <- texture_xy(c(0, 50, 50), c(50, 0, 50))
outward  <- tibble(nx = c(0, -sqrt(3) / 2, sqrt(3) / 2), ny = c(-1, 0.5, 0.5))
increase <- tibble(ux = c(-1, 0.5, 0.5), uy = c(0, sqrt(3) / 2, -sqrt(3) / 2))
title_distance <- c(11, 15, 15)
arrow_distance <- title_distance + 4
arrow_half_length <- 7
side_titles <- tibble(label = c("Sand (%)", "Clay (%)", "Silt (%)"),
                      x = side_midpoints$x + title_distance * outward$nx,
                      y = side_midpoints$y + title_distance * outward$ny,
                      angle = c(0, 60, -60))
side_arrows <- tibble(x    = side_midpoints$x + arrow_distance * outward$nx - arrow_half_length * increase$ux,
                      y    = side_midpoints$y + arrow_distance * outward$ny - arrow_half_length * increase$uy,
                      xend = side_midpoints$x + arrow_distance * outward$nx + arrow_half_length * increase$ux,
                      yend = side_midpoints$y + arrow_distance * outward$ny + arrow_half_length * increase$uy)

label_size <- 5 / ggplot2::.pt   # 5 pt
panel_texture <- ggplot() +
  geom_polygon(data = texture_class_outlines, aes(x = x, y = y, group = class),
               fill = "white", colour = "grey55", linewidth = line_width(0.35)) +
  geom_point(data = location_points, aes(x = x, y = y), colour = baobab_colour,
             size = 0.5, alpha = 0.35, stroke = 0) +
  geom_path(data = triangle_outline, aes(x = x, y = y), linewidth = line_width(0.5), colour = "black") +
  haloed_text(texture_class_labels, size = label_size, family = figure_font, lineheight = 0.9) +
  geom_segment(data = side_ticks, aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = line_width(0.5), colour = "black") +
  geom_text(data = side_ticks, aes(x = label_x, y = label_y, label = value, hjust = hjust, vjust = vjust),
            size = small_text_size_pt / ggplot2::.pt, family = figure_font, colour = "black") +
  geom_text(data = side_titles, aes(x = x, y = y, label = label, angle = angle),
            size = text_size_pt / ggplot2::.pt, family = figure_font, colour = "black") +
  geom_segment(data = side_arrows, aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(1.2, "mm"), type = "closed"),
               linewidth = line_width(0.5), colour = "black") +
  annotate("text", x = -12, y = 94, hjust = 0, vjust = 1, family = figure_font, lineheight = 1,
           size = small_text_size_pt / ggplot2::.pt, colour = "black",
           label = paste0("Topsoil (0–20 cm)\nn = ", format_number(nrow(location_points)), " baobabs")) +
  coord_equal(xlim = c(-13, 113), ylim = c(-18, 96), expand = FALSE, clip = "off") +
  theme_void(base_size = text_size_pt, base_family = figure_font) +
  theme(plot.tag = element_text(size = 8, face = "bold", colour = "black"),
        plot.margin = margin(2, 3, 2, 2, unit = "mm"))


# aridity (a), slope (b), topsoil texture (c), precipitation (d),
# temperature (e) and protected areas (f).

free_if_available <- function(panel) {
  if (exists("free", envir = asNamespace("patchwork"))) patchwork::free(panel) else panel
}

# figure 3

figure_3_panels <- list(
  panel_aridity       + labs(tag = "a"),
  panel_slope         + labs(tag = "b"),
  panel_texture       + labs(tag = "c"),
  panel_precipitation + labs(tag = "d"),
  free_if_available(panel_protected + labs(tag = "f")),
  panel_temperature   + labs(tag = "e"))
# patchwork:: is written out because terra, loaded later, has its own area().
figure_3_cells <- c(patchwork::area(1, 1), patchwork::area(1, 2),   # a, b
                    patchwork::area(2, 1), patchwork::area(2, 2),   # c, d
                    patchwork::area(3, 2), patchwork::area(3, 1))   # f (right), e (left)
figure_3 <- patchwork::wrap_plots(figure_3_panels, design = figure_3_cells, heights = c(1, 1.4, 1.3))

figure_file <- figure_output_stem
ggsave(paste0(figure_file, ".pdf"), figure_3, width = 183, height = figure_3_height_mm, units = "mm",
       bg = "white", device = cairo_pdf)
ggsave(paste0(figure_file, ".png"), figure_3, width = 183, height = figure_3_height_mm, units = "mm",
       dpi = 600, bg = "white", device = ragg::agg_png)

message("\nFigure 3 saved in ", figures_folder, "; tables saved in ", environment_results_folder)
print(as.data.frame(sources_used[, c("covariate", "source", "version_or_period")]))