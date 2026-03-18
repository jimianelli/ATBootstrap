source("R/funs.R")

survey <- Sys.getenv("ATB_SURVEY", unset = "202207")
out_path <- Sys.getenv(
  "ATB_SPATIAL_CONDITIONAL_OUT",
  unset = file.path("docs", "artifacts", "spatial_conditional_r.csv")
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

atbp <- ATBootstrapProblem(surveydata, scaling_classes = scaling_classes)

.head_string <- function(x, n = 8) {
  x <- as.integer(x)
  if (length(x) == 0) {
    return("")
  }
  paste(utils::head(x, n), collapse = " ")
}

.summary_stats <- function(x, prefix) {
  x <- as.numeric(x)
  tibble(
    !!paste0(prefix, "_n") := length(x),
    !!paste0(prefix, "_mean") := mean(x),
    !!paste0(prefix, "_sd") := stats::sd(x),
    !!paste0(prefix, "_min") := min(x),
    !!paste0(prefix, "_max") := max(x),
    !!paste0(prefix, "_nonpos_n") := sum(x <= 0)
  )
}

.derive_mu_z <- function(params, eps = .Machine$double.eps^(1 / 3)) {
  shift <- pmax(params$mu_x, eps)
  target_mean <- mean(params$data, na.rm = TRUE)
  if (!is.finite(target_mean) || target_mean <= 0) {
    target_mean <- mean(shift, na.rm = TRUE)
  }
  scaled_shift <- shift * target_mean / mean(shift, na.rm = TRUE)
  mu_z <- as.vector(forwardsolve(params$L, scaled_shift, upper.tri = FALSE, transpose = FALSE))
  mu_z[mu_z <= eps] <- eps
  list(scaled_shift = scaled_shift, mu_z = mu_z)
}

.nearest_domain_index <- function(domain_xy, obs_xy) {
  d <- .cross_distance_matrix(obs_xy, domain_xy)
  max.col(-d, ties.method = "first")
}

metrics <- bind_rows(lapply(atbp$class_problems, function(cp) {
  params <- cp$params
  obs_xy <- as.matrix(cp$geosetup$observed[, c("x", "y"), drop = FALSE])
  dom_xy <- as.matrix(cp$geosetup$domain[, c("x", "y"), drop = FALSE])
  nn_idx <- .nearest_domain_index(dom_xy, obs_xy)
  derived <- .derive_mu_z(params)
  ldiag <- diag(params$L)
  implied_var <- rowSums(params$L * params$L)

  bind_cols(
    tibble(
      backend = "r",
      survey = survey,
      class = cp$class,
      domain_n = nrow(dom_xy),
      obs_n = nrow(obs_xy),
      mapped_dloc_n = length(unique(nn_idx)),
      mapped_collision_n = nrow(obs_xy) - length(unique(nn_idx)),
      params_data_n = length(params$data),
      params_dloc_n = length(params$dlocs),
      params_sloc_n = length(params$slocs),
      dloc_head = .head_string(params$dlocs),
      sloc_head = .head_string(params$slocs),
      data_mean = mean(params$data),
      data_sd = stats::sd(params$data),
      mu = params$mu
    ),
    .summary_stats(params$mu_x, "prep_shift"),
    .summary_stats(derived$scaled_shift, "scaled_shift"),
    .summary_stats(derived$mu_z, "mu_z"),
    .summary_stats(ldiag, "L_diag"),
    .summary_stats(implied_var, "implied_var")
  )
}))

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(metrics, out_path)
cat("Wrote R spatial conditional metrics to ", out_path, "\n", sep = "")
