source("R/funs.R")

survey <- Sys.getenv("ATB_SURVEY", unset = "202207")
nsim <- as.integer(Sys.getenv("ATB_SPATIAL_NSIM", unset = "40"))
out_path <- Sys.getenv(
  "ATB_SPATIAL_PROBE_OUT",
  unset = file.path("docs", "artifacts", "spatial_backend_r.csv")
)

set.seed(as.integer(survey))

surveydir <- file.path("surveydata", survey)
sfiles <- read_survey_files(surveydir)
scaling_classes <- c("SS1", "SS2", "BT")

surveydata <- ATSurveyData(
  acoustics = subset(sfiles$acoustics, class %in% scaling_classes),
  scaling = sfiles$scaling,
  age_length = sfiles$age_length,
  length_weight = sfiles$length_weight,
  trawl_locations = sfiles$trawl_locations,
  grid = sfiles$surveygrid,
  dA = (10 * km2nmi)^2
)

atbp <- ATBootstrapProblem(
  surveydata,
  scaling_classes = scaling_classes,
  zdist_nreplicates = max(20L, min(100L, nsim))
)

det_bs <- BootSpecs(
  selectivity = FALSE,
  predict_ts = FALSE,
  resample_scaling = FALSE,
  nearbottom_coefs = FALSE,
  age_length = FALSE,
  weights_at_age = FALSE,
  trawl_assignments = FALSE,
  simulate_nasc = FALSE,
  calibration = FALSE
)

spatial_only_bs <- BootSpecs(
  selectivity = FALSE,
  predict_ts = FALSE,
  resample_scaling = FALSE,
  nearbottom_coefs = FALSE,
  age_length = FALSE,
  weights_at_age = FALSE,
  trawl_assignments = FALSE,
  simulate_nasc = TRUE,
  calibration = FALSE
)

lag_points <- seq(0, 200, by = 20)

metrics <- bind_rows(lapply(atbp$class_problems, function(cp) {
  obs_nasc <- subset(surveydata$acoustics, class == cp$class)$nasc
  det_field <- nonneg_lumult(cp$params, .zdist_mean(cp$zdists))
  sim_fields <- replicate(nsim, simulate_nasc(cp), simplify = FALSE)
  sim_means <- vapply(sim_fields, mean, numeric(1))
  sim_sds <- vapply(sim_fields, sd, numeric(1))

  det_res <- simulate_class_iteration(cp, surveydata, bs = det_bs, i = 1)
  sim_res <- lapply(seq_len(nsim), function(i) simulate_class_iteration(cp, surveydata, bs = spatial_only_bs, i = i))
  sim_total_n <- vapply(sim_res, function(df) sum(df$n, na.rm = TRUE), numeric(1))
  sim_total_biomass <- vapply(sim_res, function(df) sum(df$biomass, na.rm = TRUE), numeric(1))

  row <- tibble(
    backend = "r",
    survey = survey,
    class = cp$class,
    zdist = as.character(cp$zfamily),
    obs_mean_nasc = mean(obs_nasc),
    obs_sd_nasc = sd(obs_nasc),
    field_det_mean = mean(det_field),
    field_det_sd = sd(det_field),
    field_sim_mean_mean = mean(sim_means),
    field_sim_mean_sd = sd(sim_means),
    field_sim_sd_mean = mean(sim_sds),
    field_sim_sd_sd = sd(sim_sds),
    det_total_n = sum(det_res$n, na.rm = TRUE),
    det_total_biomass = sum(det_res$biomass, na.rm = TRUE),
    sim_total_n_mean = mean(sim_total_n),
    sim_total_n_sd = sd(sim_total_n),
    sim_total_biomass_mean = mean(sim_total_biomass),
    sim_total_biomass_sd = sd(sim_total_biomass)
  )

  model_fun <- if (is.function(cp$variogram$model)) cp$variogram$model else NULL
  if (!is.null(model_fun)) {
    gamma_vals <- model_fun(lag_points)
  } else {
    gamma_vals <- approx(cp$variogram$model$lag, cp$variogram$model$gamma, xout = lag_points, rule = 2)$y
  }

  for (i in seq_along(lag_points)) {
    row[[paste0("variogram_gamma_", lag_points[[i]])]] <- gamma_vals[[i]]
  }

  row
}))

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(metrics, out_path)
cat("Wrote R spatial backend metrics to ", out_path, "\n", sep = "")
