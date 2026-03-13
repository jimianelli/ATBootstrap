source("R/funs.R")

survey <- "202207"
source_dir <- file.path("surveydata", survey)
surveydir <- file.path(tempdir(), survey)
log_ranges <- list(c(250, 2190), c(2752.5, 2778), c(2400, 2752.49), c(2780, 6100))

dir.create(surveydir, recursive = TRUE, showWarnings = FALSE)
source_files <- list.files(source_dir, recursive = TRUE, full.names = TRUE)

for (src in source_files) {
  rel <- sub(paste0("^", source_dir, "/"), "", src)
  dest <- file.path(surveydir, rel)
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  file.copy(src, dest, overwrite = TRUE)
}

preprocess_survey_data(
  surveydir,
  dx = 10,
  ebs = TRUE,
  log_ranges = log_ranges,
  grid_method = TransectRibbons(width = 40)
)

sfiles <- read_survey_files(surveydir)
stopifnot(nrow(sfiles$acoustics) > 0)
stopifnot(nrow(sfiles$trawl_locations) > 0)
stopifnot(nrow(sfiles$surveygrid) > 0)
stopifnot(nrow(sfiles$length_weight) > 0)

surveydata <- ATSurveyData(
  acoustics = subset(sfiles$acoustics, class %in% c("SS1", "SS2", "BT")),
  scaling = sfiles$scaling,
  age_length = sfiles$age_length,
  length_weight = sfiles$length_weight,
  trawl_locations = sfiles$trawl_locations,
  grid = sfiles$surveygrid,
  dA = (10 * km2nmi)^2
)

grid_n <- nrow(sfiles$surveygrid)
cp <- ScalingClassProblem(
  class_name = "SS1",
  age_max = 10,
  aged_species = c(21740),
  simdomain = sfiles$surveygrid,
  simulator = function(scp) rep(100, grid_n)
)

res <- simulate_class_iteration(
  cp,
  surveydata,
  bs = BootSpecs(
    selectivity = FALSE,
    predict_ts = FALSE,
    resample_scaling = FALSE,
    nearbottom_coefs = FALSE,
    age_length = FALSE,
    weights_at_age = FALSE,
    trawl_assignments = FALSE,
    simulate_nasc = TRUE,
    calibration = FALSE
  ),
  i = 1
)

stopifnot(nrow(res) > 0)
stopifnot(sum(res$n, na.rm = TRUE) > 0)

cat("Non-spatial R smoke test passed.\n")
