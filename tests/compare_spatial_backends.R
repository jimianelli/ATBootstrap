suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

artifact_dir <- file.path("docs", "artifacts")
dir.create(artifact_dir, recursive = TRUE, showWarnings = FALSE)

survey <- Sys.getenv("ATB_SURVEY", unset = "202207")
nsim <- as.integer(Sys.getenv("ATB_SPATIAL_NSIM", unset = "40"))
skip_probes <- identical(Sys.getenv("ATB_SKIP_PROBES", unset = "0"), "1")
r_out <- file.path(artifact_dir, "spatial_backend_r.csv")
julia_out <- file.path(artifact_dir, "spatial_backend_julia.csv")
compare_out <- file.path(artifact_dir, "spatial_backend_parity.csv")

run_cmd <- function(command, args) {
  status <- system2(command, args = args)
  if (!identical(status, 0L)) {
    stop("Command failed: ", paste(c(command, args), collapse = " "), call. = FALSE)
  }
}

Sys.setenv(ATB_SURVEY = survey, ATB_SPATIAL_NSIM = nsim)

if (!skip_probes) {
  run_cmd(
    "Rscript",
    c("tests/spatial_backend_probe.R")
  )

  run_cmd(
    "julia",
    c("--project=.", "tests/spatial_backend_probe.jl")
  )
}

r_metrics <- read_csv(r_out, show_col_types = FALSE)
julia_metrics <- read_csv(julia_out, show_col_types = FALSE)

metric_cols <- setdiff(intersect(names(r_metrics), names(julia_metrics)), c("backend", "survey", "class", "zdist"))
variogram_cols <- grep("^variogram_gamma_", metric_cols, value = TRUE)
variogram_cols <- setdiff(variogram_cols, "variogram_gamma_0")
core_metric_cols <- setdiff(metric_cols, variogram_cols)

normalize_zdist <- function(x) {
  x <- tolower(x)
  x <- sub("^distributions\\.", "", x)
  x <- gsub("[^a-z]", "", x)
  dplyr::case_when(
    x == "gamma" ~ "gamma",
    x == "inversegaussian" ~ "inversegaussian",
    x == "inversegamma" ~ "inversegamma",
    x == "lognormal" ~ "lognormal",
    TRUE ~ x
  )
}

joined <- r_metrics %>%
  rename(zdist_r = zdist) %>%
  inner_join(julia_metrics %>% rename(zdist_julia = zdist), by = c("survey", "class"), suffix = c("_r", "_julia"))

rel_diff <- function(a, b) {
  denom <- pmax(abs(b), 1e-8)
  abs(a - b) / denom
}

comparison <- joined %>%
  transmute(
    survey,
    class,
    zdist_r,
    zdist_julia,
    zdist_match = normalize_zdist(zdist_r) == normalize_zdist(zdist_julia),
    obs_mean_nasc_rel_diff = rel_diff(obs_mean_nasc_r, obs_mean_nasc_julia),
    obs_sd_nasc_rel_diff = rel_diff(obs_sd_nasc_r, obs_sd_nasc_julia),
    field_det_mean_rel_diff = rel_diff(field_det_mean_r, field_det_mean_julia),
    field_det_sd_rel_diff = rel_diff(field_det_sd_r, field_det_sd_julia),
    field_sim_mean_mean_rel_diff = rel_diff(field_sim_mean_mean_r, field_sim_mean_mean_julia),
    field_sim_sd_mean_rel_diff = rel_diff(field_sim_sd_mean_r, field_sim_sd_mean_julia),
    det_total_n_rel_diff = rel_diff(det_total_n_r, det_total_n_julia),
    det_total_biomass_rel_diff = rel_diff(det_total_biomass_r, det_total_biomass_julia),
    sim_total_n_mean_rel_diff = rel_diff(sim_total_n_mean_r, sim_total_n_mean_julia),
    sim_total_biomass_mean_rel_diff = rel_diff(sim_total_biomass_mean_r, sim_total_biomass_mean_julia),
    variogram_max_rel_diff = do.call(
      pmax,
      c(
        lapply(variogram_cols, function(col) rel_diff(joined[[paste0(col, "_r")]], joined[[paste0(col, "_julia")]])),
        list(na.rm = TRUE)
      )
    ),
    max_rel_diff = pmax(
      field_det_mean_rel_diff,
      field_det_sd_rel_diff,
      field_sim_mean_mean_rel_diff,
      field_sim_sd_mean_rel_diff,
      det_total_n_rel_diff,
      det_total_biomass_rel_diff,
      sim_total_n_mean_rel_diff,
      sim_total_biomass_mean_rel_diff,
      variogram_max_rel_diff
    )
  ) %>%
  mutate(
    status = case_when(
      !zdist_match ~ "DIFF_ZDIST",
      max_rel_diff <= 0.1 ~ "CLOSE",
      max_rel_diff <= 0.3 ~ "PARTIAL",
      TRUE ~ "DIVERGENT"
    )
  )

write_csv(comparison, compare_out)
cat("Wrote spatial backend parity comparison to ", compare_out, "\n", sep = "")
print(comparison)
