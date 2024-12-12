library(odbc)
library(glue)
library(dplyr)
library(stringr)


cells_query <- function(survey, ship, datasetid, zonemax, zonemin) {
  sql = glue("SELECT
      macebase2.integration_results.survey,
      macebase2.integration_results.zone,
      macebase2.integration_results.interval,
      macebase2.integration_results.layer,
      macebase2.integration_results.class,
      macebase2.integration_results.data_set_id,
      macebase2.integration_results.prc_nasc
      FROM
      macebase2.integration_results
      WHERE
      macebase2.integration_results.survey ='{survey} '
      AND macebase2.integration_results.data_set_id ='{datasetid}';
      ")
  return(sql)
}

intervals_query <- function(survey, ship, datasetid, transectmax, transectmin) {
  sql = glue("SELECT
      macebase2.intervals.survey,
      macebase2.intervals.interval,
      macebase2.intervals.start_longitude,
      macebase2.intervals.start_latitude,
      macebase2.intervals.start_vessel_log,
      macebase2.intervals.mean_bottom_depth,
      macebase2.intervals.transect
      FROM
      macebase2.intervals
      WHERE
      macebase2.intervals.survey = '{survey}'
      AND macebase2.intervals.data_set_id = '{datasetid}';
      ")
  return(sql)
}

cleanup <- function(df) {
  names(df) <- tolower(names(df))
  return(df)
}

download_scaling <- function(connection, survey, datasetid, analysisid) {
  sql = glue("SELECT *
          FROM macebase2.scaling_key_source_data
          WHERE macebase2.scaling_key_source_data.survey = {survey}
          AND macebase2.scaling_key_source_data.data_set_id ='{datasetid}'
          AND macebase2.scaling_key_source_data.analysis_id ='{analysisid}';")
  print("Downloading scaling key source...")
  scaling = dbGetQuery(connection, sql)

  scaling = select(scaling, ! FILTER_ID)
  scaling = select(scaling, ! SEX)
  scaling = cleanup(scaling)

  scaling = scaling %>%
    left_join(species, by="species_code") %>%
    mutate(w = catch_sampling_expansion * user_defined_expansion * sample_correction_scalar * haul_weight)
  return(scaling)
}

download_scaling_gap <- function(connection, survey) {
  print("Dowloading GAP specimen data...")
  year = floor(survey / 100)
  query = glue("SELECT
    racebase.specimen.hauljoin,
    racebase.specimen.cruise,
    racebase.specimen.vessel,
    racebase.specimen.haul,
    racebase.specimen.specimenid,
    racebase.specimen.species_code,
    racebase.specimen.length
    FROM racebase.specimen
    WHERE floor(cruise/100) = {year};")
  scaling_gap <- dbGetQuery(connection, query)
  scaling_gap <- cleanup(scaling_gap)
  return(scaling_gap)
}

download_trawl_locations <- function(connection, survey) {
  print("Downloading trawl locations...")
  trawl_locations = dbGetQuery(connection,
                               glue("SELECT *
          FROM clamsbase2.event_data
          WHERE
          clamsbase2.event_data.survey = {survey};"))
  trawl_locations = cleanup(trawl_locations)
  trawl_locations = trawl_locations %>%
    filter(event_parameter %in% c("EQLongitude", "EQLatitude")) %>%
    select(survey, event_id, event_parameter, parameter_value) %>%
    tidyr::pivot_wider(names_from=event_parameter, values_from=parameter_value) %>%
    mutate(latitude = EQLatitude, longitude = EQLongitude) %>%
    select(survey, event_id, latitude, longitude)
  return(trawl_locations)
}

download_trawl_locations_gap <- function(connection, survey) {
  print("Dowloading GAP trawl locations...")
  year = floor(survey/100)
  query = glue("SELECT 
    racebase.haul.hauljoin, 
    racebase.haul.cruise, 
    racebase.haul.haul,
    racebase.haul.duration, 
    racebase.haul.distance_fished, 
    racebase.haul.net_width, 
    racebase.haul.start_latitude, 
    racebase.haul.start_longitude
    FROM racebase.haul
    WHERE performance >= 0
    AND haul_type in (3, 13)
    AND region = 'BS'
    AND stratum is not null
    AND stationid is not null
    and floor(cruise/100) = {year};")
  trawl_locations_gap <- dbGetQuery(connection, query)
  trawl_locations_gap <- cleanup(trawl_locations_gap) %>%
    mutate(survey = cruise, event_id = -haul,
      latitude = start_latitude, longitude=start_longitude) %>%
    select(survey, event_id, latitude, longitude)
  return(trawl_locations_gap)
}

download_length_weight_measurements <- function(connection, survey) {
  print("Downloading length-weight measurements...")
  length_weight = dbGetQuery(connection,
                            glue("SELECT *
          FROM clamsbase2.measurements
          WHERE
          clamsbase2.measurements.survey = {survey}
          AND (clamsbase2.measurements.measurement_type = 'fork_length'
            OR clamsbase2.measurements.measurement_type = 'organism_weight');"))
  species = dbGetQuery(connection,
                        glue("SELECT
          clamsbase2.samples.sample_id, 
          clamsbase2.samples.species_code
          FROM clamsbase2.samples
          WHERE 
          clamsbase2.samples.survey = {survey};"))
  length_weight <- cleanup(length_weight)
  species <- cleanup(species)
  length_weight <- left_join(length_weight, species, by="sample_id")
  return(length_weight)
}

download_age_length_key <- function(connection, survey, data_set_id, analysis_id) {
  print("Downloading age-length keys...")
  age_length_key = dbGetQuery(connection,
                               glue("SELECT *
          FROM macebase2.age_length_key_data
          WHERE
          macebase2.age_length_key_data.survey = {survey}
          AND macebase2.age_length_key_data.data_set_id = {data_set_id}
          AND macebase2.age_length_key_data.analysis_id = {analysis_id};"))
  age_length_key <- cleanup(age_length_key)
  totals = age_length_key %>%
    group_by(survey, length) %>%
    summarize(ptotal = sum(proportion_in_key), .groups="drop")

  age_length_key = age_length_key %>%
    group_by(survey, age, length) %>%
    summarize(proportion = sum(proportion_in_key), .groups="drop") %>%
    left_join(totals, by=c("survey", "length")) %>%
    mutate(proportion = proportion / ptotal) %>%
    select(survey, length, age, proportion)

  return(age_length_key)
}

download_age_length_measurements <- function(connection, survey) {
  query = glue("select a.ship, a.survey, a.event_id, a.species_code, a.specimen_id,
      a.fork_length,a.organism_weight,a.age, a.sex, b.eq_longitude 
    FROM
    clamsbase2.v_age_length a 
    JOIN clamsbase2.v_event_data b ON a.ship=b.ship
    AND
    a.survey=b.survey
    AND
    a.event_id=b.event_id
    WHERE
    a.survey = {survey} 
    AND
    a.species_code=21740;")
  age_length_measurements = dbGetQuery(connection, query)
  age_length_measurements <- cleanup(age_length_measurements)
  return(age_length_measurements)
}

fixlongitude <- function(lon) ifelse(lon > 0, -(360-lon), lon)
lookup_survey_ship <- function(survey) ifelse(survey %in% c(200608, 200408), 21, 157)
ispollock <- function(class) substr(toupper(class), 1, 2) == "PK"


download_acoustics <- function(connection, survey, ebs=TRUE) {
  print("Downloading acoustic data...")

  datasetid = 1
  zonemin = 1
  zonemax = 3
  transectmin = 1
  transectmax = 999
  ship = lookup_survey_ship(survey)

  qc = cells_query(survey, ship, datasetid, zonemax, zonemin)
  cellsdata = dbGetQuery(connection, qc)
  cellsdata = cleanup(cellsdata)

  qi = intervals_query(survey, ship, datasetid, transectmax, transectmin)
  intervalsdata = dbGetQuery(connection, qi)
  intervalsdata = cleanup(intervalsdata)
  intervalsdata = mutate(intervalsdata, start_longitute = fixlongitude(start_longitude))

  cellsdata = left_join(cellsdata, intervalsdata, by=c("survey", "interval"))
  # if this is an EBS survey, create a new virtual scaling stratum for the bottom 3 m
  if (ebs) {
    cellsdata$class[cellsdata$zone %in% c(2, 3)] = "BT"
  }

  integrated = cellsdata %>%
    # filter(zone <= zonemax) %>%
    group_by(survey, transect, interval, class) %>%
    summarize(nasc = sum(prc_nasc), .groups="drop") %>%
    left_join(intervalsdata, by=c("survey", "interval", "transect")) %>%
    mutate(transect = round(transect))
  integrated$nasc[is.na(integrated$nasc)] <- 0

  return(integrated)
}

download_survey <- function(connection, survey, data_set_id, analysis_id, ebs=TRUE) {
  print(glue("Fetching data from survey {survey}"))
  trawl_locations <- download_trawl_locations(afsc, survey)
  scaling <- download_scaling(afsc, survey, data_set_id, analysis_id)
  length_weight <- download_length_weight_measurements(afsc, survey)
  if (survey == 201608) {
    # hard-coded dataset 1 here because there is no dataset 2 for age-length data in 2016.
    # This is an ugly hack (obviously...)
    age_length_key <- download_age_length_key(afsc, survey, 1, analysis_id)
  } else {
    age_length_key <- download_age_length_key(afsc, survey, data_set_id, analysis_id)
  }
  age_length <- download_age_length_measurements(afsc, survey)
  acoustics <- download_acoustics(afsc, survey)

  surveydir = paste0("surveydata/", survey)
  dir.create(surveydir, showWarnings = FALSE, recursive = TRUE)
  write.csv(trawl_locations, paste0(surveydir, "/trawl_locations_mace.csv"))
  write.csv(scaling, paste0(surveydir, "/scaling_mace.csv"))
  write.csv(length_weight, paste0(surveydir, "/measurements.csv"))
  write.csv(age_length_key, paste0(surveydir, "/age_length_key.csv"))
  write.csv(age_length, paste0(surveydir, "/age_length.csv"))
  write.csv(acoustics, paste0(surveydir, "/acoustics.csv"))
  if (ebs) {
    scaling_gap <- download_scaling_gap(afsc, survey)
    trawl_locations_gap <- download_trawl_locations_gap(afsc, survey)
    write.csv(trawl_locations_gap, paste0(surveydir, "/trawl_locations_gap.csv"))
    write.csv(scaling_gap, paste0(surveydir, "/scaling_gap.csv"))
  }
  print(paste0("Downloaded survey ", survey, "!"))
}


cat("Enter user ID: ")
uid <- readLines("stdin", n=1)
cat(paste0("Enter password for user ", uid, ": "))
pwd <- readLines("stdin", n=1)
afsc <- dbConnect(odbc(), "AFSC", UID=uid, PWD=pwd)

species = dbGetQuery(afsc,
                    glue("SELECT
      clamsbase2.species.species_code,
      clamsbase2.species.scientific_name,
      clamsbase2.species.common_name
      FROM
      clamsbase2.species;"))
species = cleanup(species)
write.csv(species, "surveydata/species.csv")

# survey.specs <-data.frame(
#   survey = c(200707, 200809, 200909, 201006, 201207, 201407, 201608, 201807, 202207),
#   data_set_id = c(1,      1,      1,      1,      1,      1,      2,      1,      1),
#   analysis_id = c(3,      4,      4,      3,      5,      4,      1,      7,      1)
# )

# for (i in 1:nrow(survey.specs)) {
#   download_survey(afsc, 
#     survey.specs[i, "survey"], 
#     survey.specs[i, "data_set_id"], 
#     survey.specs[i, "analysis_id"])
# }

# Uncomment the line below and fill in the parameters to download a survey one-off
download_survey(afsc, survey=202104, data_set_id=3, analysis_id=1)
