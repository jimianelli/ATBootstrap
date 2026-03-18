source("R/funs.R")

sfiles <- read_survey_files("surveydata/202207")
surveydata <- ATSurveyData(
  acoustics = subset(sfiles$acoustics, class %in% c("SS1", "SS2", "BT")),
  scaling = sfiles$scaling,
  age_length = sfiles$age_length,
  length_weight = sfiles$length_weight,
  trawl_locations = sfiles$trawl_locations,
  grid = sfiles$surveygrid,
  dA = (10 * km2nmi)^2
)

atbp <- ATBootstrapProblem(
  surveydata,
  scaling_classes = "SS1",
  zdist_nreplicates = 2
)

stopifnot(length(atbp$class_problems) == 1)
cp <- atbp$class_problems[[1]]

sim <- simulate_nasc(cp)
stopifnot(length(sim) == nrow(sfiles$surveygrid))
stopifnot(all(is.finite(sim)))
stopifnot(all(sim >= 0))

common_bs <- BootSpecs(
  selectivity = FALSE,
  predict_ts = FALSE,
  resample_scaling = FALSE,
  nearbottom_coefs = FALSE,
  age_length = FALSE,
  weights_at_age = FALSE,
  trawl_assignments = FALSE,
  calibration = FALSE
)

res_sim <- simulate_class_iteration(
  cp,
  surveydata,
  bs = modifyList(common_bs, list(simulate_nasc = TRUE)),
  i = 1
)
res_det <- simulate_class_iteration(
  cp,
  surveydata,
  bs = modifyList(common_bs, list(simulate_nasc = FALSE)),
  i = 1
)

stopifnot(nrow(res_sim) > 0, nrow(res_det) > 0)
stopifnot(sum(res_sim$n, na.rm = TRUE) > 0, sum(res_det$n, na.rm = TRUE) > 0)

cat("Spatial R smoke test passed.\n")
