suppressPackageStartupMessages({
  library(readr)
  library(MASS)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(here)
  library(sf)
  library(FNN)
  library(purrr)
  library(stringr)
  library(broom)
  library(concaveman)
})

km2nmi <- 1 / 1.852
TS_SE_DEFAULT <- 3.0
nearbottom_intercept <- 3.43

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.atb_root <- function() {
  tryCatch(here::here(), error = function(...) getwd())
}

.atb_path <- function(...) {
  file.path(.atb_root(), ...)
}

.drop_helper_cols <- function(df) {
  df %>% select(-any_of(c("X", "Column1")))
}

.read_csv_local <- function(path, na = c(".", "NA")) {
  readr::read_csv(path, na = na, show_col_types = FALSE, name_repair = "unique_quiet") %>%
    .drop_helper_cols()
}

.require_columns <- function(df, cols, label) {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop(label, " is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

.as_utm_km <- function(df, lon_col, lat_col, zone = 3) {
  .require_columns(df, c(lon_col, lat_col), "coordinate input")
  sf_points <- st_as_sf(df, coords = c(lon_col, lat_col), crs = 4326, remove = FALSE)
  utm_epsg <- 32600 + zone
  utm_points <- st_transform(sf_points, utm_epsg)
  coords <- st_coordinates(utm_points) / 1000
  df$x <- coords[, 1]
  df$y <- coords[, 2]
  df
}

.normalize_lengths <- function(x, target) {
  if (length(x) == target) {
    x
  } else if (length(x) == 1L) {
    rep(x, target)
  } else {
    stop("Inputs cannot be recycled to a common length.", call. = FALSE)
  }
}

TransectRibbons <- function(width = 20, buffer = 0.1) {
  structure(list(width = width, buffer = buffer), class = "TransectRibbons")
}

SurveyHull <- function(k = 10) {
  structure(list(k = as.integer(k)), class = "SurveyHull")
}

make_haul_id <- function(prefix, ship, event_id) {
  paste0(prefix, "-", ship, "-", event_id)
}

in_intervals <- function(x, intervals) {
  result <- rep(FALSE, length(x))
  for (interval in intervals) {
    result <- result | (x >= interval[1] & x <= interval[2])
  }
  result
}

ATSurveyData <- function(acoustics, scaling, age_length, length_weight, trawl_locations, grid, dA) {
  structure(
    list(
      acoustics = acoustics,
      scaling = scaling,
      age_length = age_length,
      length_weight = length_weight,
      trawl_locations = trawl_locations,
      grid = grid,
      dA = dA
    ),
    class = "ATSurveyData"
  )
}

ScalingClassProblem <- function(class_name,
                                variogram = NULL,
                                geosetup = NULL,
                                params = NULL,
                                zfamily = NA_character_,
                                zdists = NULL,
                                cal_error = 0.1,
                                age_max = 10,
                                aged_species = c(21740),
                                simdomain = NULL,
                                simulator = NULL) {
  structure(
    list(
      class = class_name,
      variogram = variogram,
      geosetup = geosetup,
      params = params,
      zfamily = zfamily,
      zdists = zdists,
      cal_error = cal_error,
      age_max = age_max,
      aged_species = aged_species,
      simdomain = simdomain,
      simulator = simulator
    ),
    class = "ScalingClassProblem"
  )
}

ATBootstrapProblem <- function(surveydata = NULL,
                               class_problems = NULL,
                               age_max = 10,
                               aged_species = c(21740),
                               scaling_classes = NULL,
                               ...) {
  if (is.null(class_problems)) {
    stop(
      paste(
        "Spatial setup is not ported to R in this repository.",
        "Construct `class_problems` with an external spatial backend or continue using Julia for that layer."
      ),
      call. = FALSE
    )
  }

  if (is.null(scaling_classes)) {
    scaling_classes <- vapply(class_problems, function(cp) cp$class, character(1))
  }

  structure(
    list(
      class_problems = class_problems,
      scaling_classes = scaling_classes,
      age_max = age_max,
      aged_species = aged_species
    ),
    class = "ATBootstrapProblem"
  )
}

read_survey_files <- function(surveydir) {
  acoustics <- .read_csv_local(file.path(surveydir, "acoustics_projected.csv"))
  trawl_locations <- .read_csv_local(file.path(surveydir, "trawl_locations_projected.csv"))
  scaling <- .read_csv_local(file.path(surveydir, "scaling.csv")) %>%
    mutate(sample_correction_scalar = as.numeric(sample_correction_scalar))
  age_length <- .read_csv_local(file.path(surveydir, "age_length.csv"))
  length_weight <- .read_csv_local(file.path(surveydir, "length_weight.csv"))
  surveygrid <- .read_csv_local(file.path(surveydir, "surveygrid.csv"))
  surveygrid <- surveygrid[sample.int(nrow(surveygrid)), , drop = FALSE]

  list(
    acoustics = acoustics,
    scaling = scaling,
    age_length = age_length,
    length_weight = length_weight,
    trawl_locations = trawl_locations,
    surveygrid = surveygrid
  )
}

nearest_length <- function(L, unique_lengths) {
  which.min(abs(L - unique_lengths))
}

simulate_cal_error <- function(cal_error, stochastic = TRUE) {
  if (stochastic) {
    10^(cal_error * rnorm(1) / 10)
  } else {
    1
  }
}

to_linear <- function(x) {
  10^(x / 10)
}

make_age_length_function <- function(age_length, age_max = 10, stochastic = TRUE) {
  key <- age_length %>%
    mutate(
      fork_length = round(fork_length),
      age = ifelse(age > age_max, age_max, age)
    ) %>%
    count(fork_length, age, name = "n") %>%
    pivot_wider(names_from = age, values_from = n, values_fill = 0) %>%
    arrange(fork_length)

  age_cols <- setdiff(names(key), "fork_length")
  age_mat <- as.matrix(key[, age_cols, drop = FALSE])
  age_totals <- rowSums(age_mat)
  age_mat[age_totals > 0, ] <- age_mat[age_totals > 0, , drop = FALSE] / age_totals[age_totals > 0]
  lengths <- key$fork_length
  ages <- as.integer(age_cols)

  function(L) {
    vapply(L, function(len) {
      idx <- nearest_length(len, lengths)
      probs <- age_mat[idx, ]
      if (stochastic) {
        sample(ages, size = 1, prob = probs)
      } else {
        ages[which.max(probs)]
      }
    }, numeric(1))
  }
}

make_weight_function <- function(length_weight, stochastic = TRUE, nmin = 5) {
  n <- nrow(length_weight)
  idx <- if (stochastic) sample.int(n, n, replace = TRUE) else seq_len(n)
  length_weight_boot <- length_weight[idx, , drop = FALSE]

  model <- lm(log(organism_weight) ~ log(fork_length), data = length_weight_boot)
  coefs <- broom::tidy(model)$estimate
  a <- exp(coefs[1])
  b <- coefs[2]

  Lmax <- 100
  all_lengths <- tibble(fork_length = seq_len(Lmax))
  binned <- length_weight_boot %>%
    right_join(all_lengths, by = "fork_length") %>%
    group_by(fork_length) %>%
    summarise(
      mean_weight = mean(organism_weight, na.rm = TRUE),
      n = sum(!is.na(organism_weight)),
      .groups = "drop"
    ) %>%
    mutate(weight = ifelse(n >= nmin, mean_weight, a * fork_length^b))

  weight_dict <- setNames(binned$weight, binned$fork_length)

  function(L) {
    L1 <- pmin(pmax(round(L), 1), Lmax)
    unname(weight_dict[as.character(L1)])
  }
}

age0_pollock_ts <- function(L) 20 * log10(L) - 64.86
arctic_cod_ts <- function(L) 8.03 * log10(L) - 60.78
capelin_ts <- function(L) 20 * log10(L) - 70.3
chrysaora_melanaster_ts <- function(L) 10 * log10(pi * (2 * L)^2) - 86.8
eulachon_ts <- function(L) 20 * log10(L) - 84.5
generalized_physoclist_ts <- function(L) 20 * log10(L) - 67.4
generic_fish_no_swimbladder_ts <- function(L) 20 * log10(L) - 83.2
herring_75m_v2_ts <- function(L) 20 * log10(L) - log10(1 + 75 / 10) - 65.4
myctophids_sleucopsarus_ts <- function(L) 32.1 * log10(log10(L)) - 64.1
pacific_hake_ts <- function(L) 20 * log10(L) - 68.0
sandlance_ts <- function(L) 56.5 * log10(L) - 125.1
squids_ts <- function(L) 20 * log10(L) - 75.4
standard_pollock_ts <- function(L) 20 * log10(L) - 66

euphausiids_15_65mm_38khz <- function(L) {
  A <- -9.30429983e2
  B <- 3.21027896
  C <- 1.74003785
  D <- 1.36133896e-8
  E <- -2.26958555e-6
  F <- 1.50291244e-4
  G <- -4.86306872e-3
  H <- 7.38748423e-2
  I <- -4.08004891e-1
  J <- -7.39078690e1
  Lo <- 3.835e-2
  c <- 1470
  nu <- 38
  k <- 2 * pi * (nu * 10^3) / c

  out <- numeric(length(L))
  out[L < 1.5] <- -105
  out[L > 6.5] <- -73

  mid <- L >= 1.5 & L <= 6.5
  Lm <- L[mid] / 100
  out[mid] <- A * (log10(B * k * Lm) / (B * k * Lm))^C +
    D * (k * Lm)^6 +
    E * (k * Lm)^5 +
    F * (k * Lm)^4 +
    G * (k * Lm)^3 +
    H * (k * Lm)^2 +
    I * (k * Lm) +
    J +
    20 * log10(Lm / Lo)
  out
}

ts_lookup <- list(
  age0_pollock = list(func = age0_pollock_ts, se = TS_SE_DEFAULT),
  arctic_cod = list(func = arctic_cod_ts, se = TS_SE_DEFAULT),
  capelin = list(func = capelin_ts, se = TS_SE_DEFAULT),
  chrysaora_melanaster = list(func = chrysaora_melanaster_ts, se = TS_SE_DEFAULT),
  eulachon = list(func = eulachon_ts, se = TS_SE_DEFAULT),
  eulachon_new = list(func = eulachon_ts, se = TS_SE_DEFAULT),
  euphausiids_15_65mm_38khz = list(func = euphausiids_15_65mm_38khz, se = TS_SE_DEFAULT),
  generalized_physoclist = list(func = generalized_physoclist_ts, se = TS_SE_DEFAULT),
  generic_fish_no_swimbladder = list(func = generic_fish_no_swimbladder_ts, se = TS_SE_DEFAULT),
  generic_swimbladder_fish = list(func = generalized_physoclist_ts, se = TS_SE_DEFAULT),
  herring_75m_v2 = list(func = herring_75m_v2_ts, se = TS_SE_DEFAULT),
  myctophids_sleucopsarus = list(func = myctophids_sleucopsarus_ts, se = TS_SE_DEFAULT),
  pacific_hake = list(func = pacific_hake_ts, se = TS_SE_DEFAULT),
  sandlance = list(func = sandlance_ts, se = TS_SE_DEFAULT),
  squids = list(func = squids_ts, se = TS_SE_DEFAULT),
  standard_pollock = list(func = standard_pollock_ts, se = 0.14)
)

make_ts_function <- function() {
  error_dict <- map_dbl(names(ts_lookup), ~ rnorm(1) * ts_lookup[[.x]]$se)
  names(error_dict) <- names(ts_lookup)

  function(relationship, L, stochastic = FALSE) {
    n <- max(length(relationship), length(L))
    relationship <- .normalize_lengths(relationship, n)
    L <- .normalize_lengths(L, n)
    stochastic <- .normalize_lengths(stochastic, n)

    vapply(seq_len(n), function(i) {
      ts_spec <- ts_lookup[[relationship[i]]]
      if (is.null(ts_spec)) {
        stop("Unknown TS relationship: ", relationship[i], call. = FALSE)
      }
      err <- if (stochastic[i]) error_dict[[relationship[i]]] else 0
      ts_spec$func(L[i]) + err
    }, numeric(1))
  }
}

make_selectivity_function <- function(stochastic = TRUE) {
  mu_lfs <- c(-2.4664253, 0.2698528)
  sigma_lfs <- matrix(c(0.41987457, -0.0233237, -0.0233237, 0.001459591), nrow = 2)
  mu_awt <- c(-1.0558410, 0.1741619)
  sigma_awt <- matrix(c(0.49771768, -0.01810365, -0.01810365, 0.0008656043), nrow = 2)

  beta_lfs <- if (stochastic) MASS::mvrnorm(1, mu_lfs, sigma_lfs) else mu_lfs
  beta_awt <- if (stochastic) MASS::mvrnorm(1, mu_awt, sigma_awt) else mu_awt

  function(L, survey) {
    lfs <- exp(beta_lfs[1] + L * beta_lfs[2]) / (1 + exp(beta_lfs[1] + L * beta_lfs[2]))
    awt <- exp(beta_awt[1] + L * beta_awt[2]) / (1 + exp(beta_awt[1] + L * beta_awt[2]))
    ifelse(survey < 202001, awt, lfs)
  }
}

apply_selectivity <- function(scaling, selectivity_function) {
  scaling %>%
    mutate(
      selectivity = ifelse(species_code == 21740 & class != "BT", selectivity_function(primary_length, survey), 1),
      sample_correction_scalar = ifelse(species_code == 21740 & class != "BT", 1 / selectivity, sample_correction_scalar),
      w = catch_sampling_expansion * user_defined_expansion * sample_correction_scalar * haul_weight
    ) %>%
    select(-selectivity)
}

.nearbottom_species <- function() {
  .read_csv_local(.atb_path("surveydata", "species.csv")) %>%
    mutate(
      nearbottom_group = case_when(
        species_code == 21740 ~ "Pollock",
        species_code == 21725 ~ "Arctic cod",
        between(species_code, 10110, 10120) ~ "Large flatfish",
        between(species_code, 30050, 30535) ~ "Rockfish",
        TRUE ~ "Misc"
      )
    )
}

make_nearbottom_dict <- function(stochastic = TRUE) {
  species <- .nearbottom_species()
  nearbottom_coefs <- tibble(
    nearbottom_group = c("Pollock", "Arctic cod", "Large flatfish", "Rockfish", "Misc"),
    A = c(2.52, 16.39, 0.85, 93.59, 11.63),
    lower = c(2.21, 2.84, 0.16, 8.63, 3.92),
    upper = c(2.86, 62.26, 1.67, 343.43, 22.09),
    b = c(66, 67.4, 67.4, 67.4, 67.4)
  ) %>%
    mutate(sd_A = ((A - lower) + (upper - A)) / 4) %>%
    left_join(species, by = "nearbottom_group")

  aa <- if (stochastic) {
    pmap_dbl(
      list(nearbottom_coefs$A, nearbottom_coefs$sd_A, nearbottom_coefs$lower, nearbottom_coefs$upper),
      function(mu, sigma, lower, upper) {
        repeat {
          draw <- rnorm(1, mu, sigma)
          if (draw >= lower && draw <= upper) {
            return(draw)
          }
        }
      }
    )
  } else {
    nearbottom_coefs$A
  }

  aa <- aa / 10^(-nearbottom_coefs$b / 10) / 1852^2 / (4 * pi)
  setNames(aa, nearbottom_coefs$species_code)
}

apply_nearbottom_coefficient <- function(scaling, nearbottom_dict) {
  scaling %>%
    mutate(
      user_defined_expansion = ifelse(
        as.character(species_code) %in% names(nearbottom_dict),
        unname(nearbottom_dict[as.character(species_code)]),
        unname(nearbottom_dict[["24166"]])
      ),
      w = catch_sampling_expansion * user_defined_expansion * sample_correction_scalar * haul_weight
    )
}

BootSpecs <- function(selectivity = TRUE,
                      predict_ts = TRUE,
                      resample_scaling = TRUE,
                      nearbottom_coefs = TRUE,
                      age_length = TRUE,
                      weights_at_age = TRUE,
                      trawl_assignments = TRUE,
                      simulate_nasc = TRUE,
                      calibration = TRUE) {
  list(
    selectivity = selectivity,
    predict_ts = predict_ts,
    resample_scaling = resample_scaling,
    nearbottom_coefs = nearbottom_coefs,
    age_length = age_length,
    weights_at_age = weights_at_age,
    trawl_assignments = trawl_assignments,
    simulate_nasc = simulate_nasc,
    calibration = calibration
  )
}

error_labels <- tibble(
  added_error = c(
    "calibration",
    "simulate_nasc",
    "selectivity",
    "resample_scaling",
    "nearbottom_coefs",
    "trawl_assignments",
    "predict_ts",
    "age_length",
    "weights_at_age",
    "All"
  ),
  error_label = factor(
    c(
      "Calibration",
      "Spatial sampling",
      "Selectivity",
      "Resample catches",
      "Nearbottom coefs",
      "Trawl assignment",
      "TS models",
      "Age-length",
      "Length-weight",
      "All"
    ),
    levels = c(
      "Calibration",
      "Spatial sampling",
      "Selectivity",
      "Resample catches",
      "Nearbottom coefs",
      "Trawl assignment",
      "TS models",
      "Age-length",
      "Length-weight",
      "All"
    )
  )
)

resample_df <- function(df, stochastic = TRUE) {
  n <- nrow(df)
  if (stochastic) df[sample.int(n, n, replace = TRUE), , drop = FALSE] else df
}

resample_scaling <- function(df, stochastic = TRUE) {
  df %>%
    group_by(haul_id, class) %>%
    group_modify(~ resample_df(.x, stochastic)) %>%
    ungroup()
}

trawl_assignments <- function(pixel_coords, trawl_coords, stochastic = TRUE, nneighbors = 4, a = 2) {
  pixel_coords <- as.matrix(pixel_coords[, c("x", "y"), drop = FALSE])
  trawl_coords <- as.matrix(trawl_coords[, c("x", "y"), drop = FALSE])
  nneighbors <- min(nneighbors, nrow(trawl_coords))
  nn_result <- FNN::get.knnx(trawl_coords, pixel_coords, k = nneighbors)

  assignments <- integer(nrow(pixel_coords))
  for (i in seq_len(nrow(pixel_coords))) {
    idx <- nn_result$nn.index[i, ]
    dists <- nn_result$nn.dist[i, ]
    if (stochastic) {
      if (any(dists == 0)) {
        choices <- idx[dists == 0]
        assignments[i] <- sample(choices, 1)
      } else {
        assignments[i] <- sample(idx, 1, prob = dists^(-a))
      }
    } else {
      assignments[i] <- idx[which.min(dists)]
    }
  }
  assignments
}

get_trawl_category_means <- function(scaling, aged_species, predict_weight) {
  scaling %>%
    mutate(
      category = ifelse(species_code %in% aged_species, paste0(species_code, "@", age), paste0(species_code, "@-1")),
      weight = predict_weight(primary_length),
      p_nasc = sigma_bs * w
    ) %>%
    filter(w > 0) %>%
    group_by(haul_id, category) %>%
    summarise(
      sigma_bs = mean(sigma_bs),
      p_nasc = sum(p_nasc),
      weight = mean(weight),
      .groups = "drop"
    ) %>%
    group_by(haul_id) %>%
    mutate(p_nasc = p_nasc / sum(p_nasc)) %>%
    ungroup()
}

summarize_bootstrap <- function(results, variable = "n", species_codes = 21740) {
  results %>%
    filter(species_code %in% species_codes) %>%
    arrange(age) %>%
    group_by(age) %>%
    summarise(
      mean = mean(.data[[variable]], na.rm = TRUE),
      std = sd(.data[[variable]], na.rm = TRUE),
      cv = std / mean * 100,
      .groups = "drop"
    ) %>%
    rename(!!variable := mean)
}

summarize_stepwise_bootstrap <- function(results, variable = "n", species_codes = 21740) {
  results %>%
    filter(species_code %in% species_codes) %>%
    group_by(added_error, age) %>%
    summarise(
      mean = mean(.data[[variable]], na.rm = TRUE),
      std = sd(.data[[variable]], na.rm = TRUE),
      cv = std / mean * 100,
      .groups = "drop"
    ) %>%
    rename(!!variable := mean)
}

merge_results <- function(results, results_step) {
  results_all <- results %>%
    group_by(i) %>%
    summarise(n = sum(n, na.rm = TRUE), biomass = sum(biomass, na.rm = TRUE), .groups = "drop") %>%
    mutate(added_error = "All")

  stepwise_totals <- results_step %>%
    group_by(added_error, i) %>%
    summarise(n = sum(n, na.rm = TRUE), biomass = sum(biomass, na.rm = TRUE), .groups = "drop")

  bind_rows(stepwise_totals, results_all) %>%
    left_join(error_labels, by = "added_error")
}

.zdist_mean <- function(z) {
  if (is.numeric(z) && length(z) == 1L) {
    z
  } else if (is.list(z) && !is.null(z$mean)) {
    z$mean
  } else if (!is.null(attr(z, "mean"))) {
    attr(z, "mean")
  } else {
    stop("Cannot infer z-distribution mean from supplied `zdists`.", call. = FALSE)
  }
}

nonneg_lumult <- function(params, z) {
  required <- c("data", "L", "dlocs", "slocs")
  if (!all(required %in% names(params))) {
    stop("`params` must contain: ", paste(required, collapse = ", "), call. = FALSE)
  }

  npts <- length(params$dlocs) + length(params$slocs)
  x <- numeric(npts)
  x[params$slocs] <- as.vector(params$L %*% z)
  x[params$dlocs] <- params$data
  x
}

nonneg_lusim <- function(params, zdists) {
  z <- vapply(zdists, function(zdist) {
    if (is.function(zdist)) {
      zdist()
    } else if (is.numeric(zdist) && length(zdist) == 1L) {
      zdist
    } else if (is.list(zdist) && is.function(zdist$draw)) {
      zdist$draw()
    } else {
      stop("Unsupported z-distribution representation.", call. = FALSE)
    }
  }, numeric(1))
  nonneg_lumult(params, z)
}

simulate_nasc <- function(scp) {
  if (is.function(scp$simulator)) {
    return(scp$simulator(scp))
  }
  if (!is.null(scp$params) && !is.null(scp$zdists)) {
    return(nonneg_lusim(scp$params, scp$zdists))
  }
  stop("No R spatial simulator is available for this `ScalingClassProblem`.", call. = FALSE)
}

solution_domain <- function(scp) {
  domain <- scp$simdomain %||% scp$geosetup$domain
  if (is.null(domain)) {
    stop("No simulation domain is attached to this `ScalingClassProblem`.", call. = FALSE)
  }

  if (inherits(domain, "sf")) {
    coords <- st_coordinates(domain)
    return(tibble(x = coords[, 1], y = coords[, 2]))
  }
  if (is.data.frame(domain)) {
    .require_columns(domain, c("x", "y"), "simulation domain")
    return(domain[, c("x", "y")])
  }
  if (is.matrix(domain) && ncol(domain) >= 2) {
    return(tibble(x = domain[, 1], y = domain[, 2]))
  }

  stop("Unsupported simulation domain representation.", call. = FALSE)
}

zdists <- function(atbp) {
  tibble(
    class = vapply(atbp$class_problems, function(cp) cp$class, character(1)),
    zdist = vapply(atbp$class_problems, function(cp) as.character(cp$zfamily %||% NA_character_), character(1))
  )
}

plot_class_variograms <- function(atbp) {
  plots <- lapply(atbp$class_problems, function(cp) {
    empirical <- cp$variogram$empirical %||% cp$variogram
    if (is.null(empirical) || !is.data.frame(empirical)) {
      stop("`plot_class_variograms` requires empirical variogram data frames from a spatial backend.", call. = FALSE)
    }

    lag_col <- intersect(c("lag", "distance", "x"), names(empirical))[1]
    gamma_col <- intersect(c("gamma", "ordinate", "y"), names(empirical))[1]
    if (is.na(lag_col) || is.na(gamma_col)) {
      stop("Empirical variogram data must include lag and gamma columns.", call. = FALSE)
    }

    p <- ggplot(empirical, aes(x = .data[[lag_col]], y = .data[[gamma_col]])) +
      geom_point() +
      labs(title = cp$class, x = "Lag (km)", y = "\u03b3") +
      theme_minimal()

    model <- cp$variogram$model %||% NULL
    if (is.function(model)) {
      xseq <- seq(min(empirical[[lag_col]], na.rm = TRUE), max(empirical[[lag_col]], na.rm = TRUE), length.out = 200)
      p <- p + geom_line(data = tibble(x = xseq, y = model(xseq)), aes(x = x, y = y), inherit.aes = FALSE)
    } else if (is.data.frame(model)) {
      mx <- intersect(c("lag", "distance", "x"), names(model))[1]
      my <- intersect(c("gamma", "ordinate", "y"), names(model))[1]
      p <- p + geom_line(data = model, aes(x = .data[[mx]], y = .data[[my]]), inherit.aes = FALSE)
    }
    p
  })
  plots
}

plot_simulated_nasc <- function(x, surveydata, ...) {
  if (inherits(x, "ATBootstrapProblem")) {
    return(lapply(x$class_problems, plot_simulated_nasc, surveydata = surveydata, ...))
  }

  dom <- solution_domain(x)
  sim_field <- simulate_nasc(x)
  obs <- surveydata$acoustics %>% filter(class == x$class)

  ggplot(tibble(x = dom$x, y = dom$y, nasc = sim_field), aes(x = x, y = y, color = nasc)) +
    geom_point(shape = 15, size = 1.2) +
    geom_point(data = obs, aes(x = x, y = y, size = nasc), inherit.aes = FALSE, color = "white", alpha = 0.35) +
    scale_color_viridis_c() +
    coord_equal() +
    labs(title = x$class, x = "Easting (km)", y = "Northing (km)") +
    theme_minimal()
}

plot_geosim_stats <- function(atbp, surveydata, n = 200) {
  stats <- map_dfr(atbp$class_problems, function(cp) {
    obs_nasc <- surveydata$acoustics %>% filter(class == cp$class) %>% pull(nasc)
    sims <- replicate(n, simulate_nasc(cp), simplify = FALSE)
    tibble(
      class = cp$class,
      obs_mean = mean(obs_nasc),
      obs_sd = sd(obs_nasc),
      sim_mean = vapply(sims, mean, numeric(1)),
      sim_sd = vapply(sims, sd, numeric(1))
    )
  })

  stats_long <- bind_rows(
    stats %>% transmute(class, metric = "Mean", observed = obs_mean, simulated = sim_mean),
    stats %>% transmute(class, metric = "SD", observed = obs_sd, simulated = sim_sd)
  )

  ggplot(stats_long, aes(x = simulated)) +
    geom_histogram(bins = 30, fill = "#5a8f7b", color = "white") +
    geom_vline(aes(xintercept = observed), color = "#c2410c", linewidth = 1) +
    facet_grid(metric ~ class, scales = "free_x") +
    theme_minimal() +
    labs(x = "Simulated statistic", y = "Count")
}

plot_boot_results <- function(results) {
  results %>%
    filter(species_code == 21740) %>%
    pivot_longer(c(n, biomass), names_to = "variable", values_to = "value") %>%
    ggplot(aes(x = factor(age), y = value)) +
    geom_violin(fill = "#0f766e", alpha = 0.55, color = NA) +
    facet_wrap(~ variable, scales = "free_y") +
    theme_minimal() +
    labs(x = "Age", y = "Bootstrap value")
}

plot_error_source_by_age <- function(results_step, results, variable = "n", species_codes = 21740) {
  baseline <- summarize_bootstrap(results, variable = variable, species_codes = species_codes) %>%
    mutate(added_error = "All")
  stepwise <- summarize_stepwise_bootstrap(results_step, variable = variable, species_codes = species_codes)

  bind_rows(stepwise, baseline) %>%
    ggplot(aes(x = age, y = cv, color = added_error)) +
    geom_line() +
    geom_point(size = 1.5) +
    theme_minimal() +
    labs(x = "Age", y = "CV (%)", color = "Error source")
}

plot_error_sources <- function(results_totals, xlims = NULL, plot_title = NULL) {
  summary <- results_totals %>%
    pivot_longer(c(n, biomass), names_to = "variable", values_to = "value") %>%
    group_by(added_error, error_label, variable) %>%
    summarise(cv = sd(value, na.rm = TRUE) / mean(value, na.rm = TRUE), .groups = "drop")

  p <- ggplot(summary, aes(x = cv, y = forcats::fct_reorder(error_label, cv), color = variable)) +
    geom_point(size = 2.2) +
    theme_minimal() +
    labs(x = "CV", y = NULL, title = plot_title, color = NULL)

  if (!is.null(xlims)) {
    p <- p + coord_cartesian(xlim = xlims)
  }
  p
}

simulate_class_iteration <- function(scp, surveydata, bs = BootSpecs(), i = 1) {
  scaling_sub <- surveydata$scaling %>% filter(class == scp$class)

  selectivity_function <- make_selectivity_function(bs$selectivity)
  scaling_boot <- resample_scaling(scaling_sub, bs$resample_scaling)
  scaling_boot <- apply_selectivity(scaling_boot, selectivity_function)

  predict_ts <- make_ts_function()
  predict_age <- make_age_length_function(surveydata$age_length, scp$age_max, bs$age_length)
  scaling_boot <- scaling_boot %>%
    mutate(
      sigma_bs = to_linear(predict_ts(ts_relationship, ts_length, bs$predict_ts)),
      age = predict_age(primary_length)
    )

  geotrawls <- scaling_boot %>%
    group_by(haul_id) %>%
    summarise(class = first(class), .groups = "drop") %>%
    inner_join(surveydata$trawl_locations, by = "haul_id") %>%
    select(haul_id, x, y)

  ii <- trawl_assignments(surveydata$grid, geotrawls, bs$trawl_assignments)

  nasc_val <- if (bs$simulate_nasc) {
    simulate_nasc(scp)
  } else {
    if (is.null(scp$zdists) || is.null(scp$params)) {
      stop("Deterministic spatial path requires `params` and `zdists` on the class problem.", call. = FALSE)
    }
    z0 <- vapply(scp$zdists, .zdist_mean, numeric(1))
    nonneg_lumult(scp$params, z0)
  }

  cal_error_sim <- simulate_cal_error(scp$cal_error, bs$calibration)
  nasc_df <- tibble(
    nasc = nasc_val * cal_error_sim,
    haul_id = geotrawls$haul_id[ii]
  )

  if (scp$class == "BT") {
    nearbottom_dict <- make_nearbottom_dict(bs$nearbottom_coefs)
    scaling_boot <- apply_nearbottom_coefficient(scaling_boot, nearbottom_dict) %>%
      mutate(
        sigma_bs = ifelse(
          species_code == 21740,
          sigma_bs,
          to_linear(predict_ts(ts_relationship, ts_length, FALSE))
        )
      )
    nasc_df <- nasc_df %>%
      mutate(nasc = pmax(nasc - nearbottom_intercept, 0))
  }

  predict_weight <- make_weight_function(surveydata$length_weight, stochastic = bs$weights_at_age)
  trawl_means_cat <- get_trawl_category_means(scaling_boot, scp$aged_species, predict_weight)

  category_nasc <- trawl_means_cat %>%
    select(haul_id, category, p_nasc) %>%
    pivot_wider(names_from = category, values_from = p_nasc, values_fill = 0)
  category_sigma <- trawl_means_cat %>% select(haul_id, category, sigma_bs)
  category_weight <- trawl_means_cat %>% select(haul_id, category, weight)

  nasc_df %>%
    left_join(category_nasc, by = "haul_id") %>%
    pivot_longer(-c(haul_id, nasc), names_to = "category", values_to = "p_nasc") %>%
    left_join(category_sigma, by = c("haul_id", "category")) %>%
    left_join(category_weight, by = c("haul_id", "category")) %>%
    mutate(
      n = nasc * p_nasc / (4 * pi * sigma_bs) * surveydata$dA,
      biomass = n * weight
    ) %>%
    group_by(category) %>%
    summarise(n = sum(n, na.rm = TRUE), biomass = sum(biomass, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      species_code = as.integer(str_extract(category, "^[0-9]+")),
      age = as.integer(str_extract(category, "-?[0-9]+$")),
      i = i,
      class = scp$class
    ) %>%
    select(class, i, species_code, age, category, n, biomass) %>%
    arrange(i, species_code, age)
}

simulate_class <- function(scp, surveydata, nreplicates = 500, bs = BootSpecs()) {
  bind_rows(lapply(seq_len(nreplicates), function(i) simulate_class_iteration(scp, surveydata, bs, i)))
}

simulate <- function(atbp,
                     surveydata,
                     nreplicates = 500,
                     bs = BootSpecs(),
                     report_species = c(21740),
                     report_ages = seq_len(atbp$age_max)) {
  class_results <- lapply(atbp$class_problems, simulate_class, surveydata = surveydata, nreplicates = nreplicates, bs = bs)

  if (length(report_species) == 0) {
    report_species <- unique(surveydata$scaling$species_code)
  }
  if (length(report_ages) == 0) {
    report_ages <- c(-1, seq_len(atbp$age_max))
  }

  bind_rows(class_results) %>%
    filter(age %in% report_ages, species_code %in% report_species) %>%
    group_by(i, species_code, age) %>%
    summarise(n = sum(n, na.rm = TRUE), biomass = sum(biomass, na.rm = TRUE), .groups = "drop") %>%
    mutate(age = ifelse(age == -1, NA_integer_, age))
}

stepwise_error <- function(atbp, surveydata, nreplicates = 500, remove = FALSE) {
  error_sources <- names(BootSpecs())
  colname <- if (remove) "eliminated_error" else "added_error"

  bind_rows(lapply(seq_along(error_sources), function(i) {
    errs <- rep(remove, length(error_sources))
    errs[i] <- !remove
    names(errs) <- error_sources
    bs <- do.call(BootSpecs, as.list(errs))
    res <- simulate(atbp, surveydata, nreplicates = nreplicates, bs = bs)
    res[[colname]] <- error_sources[i]
    res
  }))
}

.transect_ribbon_geometry <- function(transect, method, order_col, dx) {
  ord_values <- round(transect[[order_col]] / dx) * dx
  line_df <- tibble(order_value = ord_values, x = transect$x, y = transect$y) %>%
    group_by(order_value) %>%
    summarise(x = mean(x), y = mean(y), .groups = "drop") %>%
    arrange(order_value)

  geom <- if (nrow(line_df) > 1) {
    st_sfc(st_linestring(as.matrix(line_df[, c("x", "y")])), crs = 32603)
  } else {
    st_as_sf(line_df, coords = c("x", "y"), crs = 32603) %>% st_geometry()
  }

  st_buffer(geom, dist = method$width / 2 * 1.852 * (1 + method$buffer), endCapStyle = "ROUND", joinStyle = "ROUND")
}

survey_domain <- function(acoustics, method = TransectRibbons(), order = "y", dx = 10, dy = dx) {
  UseMethod("survey_domain", method)
}

survey_domain.TransectRibbons <- function(acoustics, method = TransectRibbons(), order = "y", dx = 10, dy = dx) {
  geoms <- lapply(split(acoustics, acoustics$transect), .transect_ribbon_geometry, method = method, order_col = order, dx = dx)
  st_as_sf(st_union(do.call(c, geoms)))
}

survey_domain.SurveyHull <- function(acoustics, method = SurveyHull(), order = "y", dx = 10, dy = dx) {
  points <- st_as_sf(acoustics %>% distinct(x, y), coords = c("x", "y"), crs = 32603)
  hull <- tryCatch(concaveman::concaveman(points), error = function(...) st_convex_hull(st_union(points)))
  st_as_sf(hull)
}

grid_domain <- function(domain_sf, dx, dy = dx) {
  bbox <- st_bbox(domain_sf)
  grid <- expand.grid(
    x = seq(bbox[["xmin"]], bbox[["xmax"]], by = dx),
    y = seq(bbox[["ymin"]], bbox[["ymax"]], by = dy)
  )
  grid_sf <- st_as_sf(grid, coords = c("x", "y"), crs = st_crs(domain_sf))
  inside <- lengths(st_within(grid_sf, domain_sf)) > 0
  coords <- st_coordinates(grid_sf[inside, ])
  tibble(x = coords[, 1], y = coords[, 2])
}

get_survey_grid <- function(acoustics, method = TransectRibbons(), dx = 10, dy = dx, order = "y") {
  if (nrow(acoustics) <= 3) {
    stop("Not enough locations to define a survey grid.", call. = FALSE)
  }
  domain <- survey_domain(acoustics, method = method, order = order, dx = dx, dy = dy)
  list(grid = grid_domain(domain, dx = dx, dy = dy), domain = domain)
}

merge_scaling <- function(scaling_mace, scaling_gap) {
  ts_key <- scaling_mace %>%
    group_by(species_code) %>%
    summarise(ts_relationship = first(ts_relationship), .groups = "drop")

  scaling_mace1 <- scaling_mace %>%
    mutate(haul_id = make_haul_id("MACE", 1, event_id)) %>%
    select(
      survey, ship, haul_id, class, species_code, primary_length, ts_length, ts_relationship,
      catch_sampling_expansion, user_defined_expansion, sample_correction_scalar, haul_weight, w
    )

  scaling_gap1 <- scaling_gap %>%
    mutate(
      survey = cruise,
      ship = vessel,
      haul_id = make_haul_id("GAP", vessel, haul),
      class = "BT",
      primary_length = length / 10,
      ts_length = length / 10,
      catch_sampling_expansion = 1,
      user_defined_expansion = 1,
      sample_correction_scalar = 1,
      haul_weight = 1,
      w = 1
    ) %>%
    left_join(ts_key, by = "species_code") %>%
    mutate(ts_relationship = replace_na(ts_relationship, "generic_swimbladder_fish")) %>%
    select(
      survey, ship, haul_id, class, species_code, primary_length, ts_length, ts_relationship,
      catch_sampling_expansion, user_defined_expansion, sample_correction_scalar, haul_weight, w
    ) %>%
    drop_na()

  bind_rows(scaling_mace1, scaling_gap1)
}

merge_trawl_locations <- function(trawl_locations_mace, trawl_locations_gap) {
  survey <- unique(trawl_locations_mace$survey)
  if (length(survey) != 1L) {
    stop("Expected a single survey in trawl location inputs.", call. = FALSE)
  }

  tl_mace <- trawl_locations_mace %>%
    mutate(haul_id = make_haul_id("MACE", 1, event_id)) %>%
    select(survey, haul_id, latitude, longitude)

  tl_gap <- trawl_locations_gap %>%
    mutate(survey = survey[[1]], haul_id = make_haul_id("GAP", vessel, event_id)) %>%
    select(survey, haul_id, latitude, longitude)

  bind_rows(tl_mace, tl_gap)
}

preprocess_survey_data <- function(surveydir,
                                   ebs = TRUE,
                                   log_ranges = NULL,
                                   dx = 10,
                                   dy = dx,
                                   missingstring = c(".", "NA"),
                                   grid_method = TransectRibbons(),
                                   transect_order = "y") {
  scaling_mace <- .read_csv_local(file.path(surveydir, "scaling_mace.csv"), na = missingstring)
  trawl_locations_mace <- .read_csv_local(file.path(surveydir, "trawl_locations_mace.csv"), na = missingstring)

  if (ebs) {
    scaling_gap <- .read_csv_local(file.path(surveydir, "scaling_gap.csv"), na = missingstring)
    trawl_locations_gap <- .read_csv_local(file.path(surveydir, "trawl_locations_gap.csv"), na = missingstring)
    scaling <- merge_scaling(scaling_mace, scaling_gap)
    trawl_locations <- merge_trawl_locations(trawl_locations_mace, trawl_locations_gap)
  } else {
    scaling <- scaling_mace
    trawl_locations <- trawl_locations_mace
  }

  acoustics <- .read_csv_local(file.path(surveydir, "acoustics.csv"), na = missingstring)

  scaling$class <- str_replace_all(scaling$class, "_FILTERED", "")
  acoustics$class <- str_replace_all(acoustics$class, "_FILTERED", "")
  scaling_classes <- unique(scaling$class)

  if (is.null(log_ranges)) {
    log_ranges <- list(range(acoustics$start_vessel_log, na.rm = TRUE))
  }

  acoustics <- acoustics %>%
    transmute(
      transect = transect,
      interval = interval,
      class = class,
      lat = start_latitude,
      lon = start_longitude,
      log = start_vessel_log,
      nasc = nasc
    ) %>%
    filter(
      class %in% scaling_classes,
      in_intervals(log, log_ranges),
      abs(lon) < 360,
      abs(lat) < 360
    ) %>%
    group_by(transect, interval, class, lat, lon, log) %>%
    summarise(nasc = sum(nasc), .groups = "drop") %>%
    pivot_wider(names_from = class, values_from = nasc, values_fill = 0) %>%
    pivot_longer(-c(transect, interval, lat, lon, log), names_to = "class", values_to = "nasc")

  acoustics <- .as_utm_km(acoustics, "lon", "lat", zone = 3)
  grid_info <- get_survey_grid(acoustics, method = grid_method, dx = dx, dy = dy, order = transect_order)
  surveygrid <- grid_info$grid
  surveyhull <- grid_info$domain

  acoustics <- acoustics %>%
    mutate(
      x = round(x / dx) * dx,
      y = round(y / dy) * dy
    ) %>%
    group_by(transect, class, x, y) %>%
    summarise(
      lon = mean(lon),
      lat = mean(lat),
      log = mean(log),
      nasc = mean(nasc),
      .groups = "drop"
    )

  trawl_locations <- .as_utm_km(trawl_locations, "longitude", "latitude", zone = 3)
  trawl_sf <- st_as_sf(trawl_locations, coords = c("x", "y"), crs = 32603, remove = FALSE)
  inside <- lengths(st_within(trawl_sf, surveyhull)) > 0
  trawl_locations <- trawl_locations[inside, , drop = FALSE]

  length_weight <- .read_csv_local(file.path(surveydir, "measurements.csv"), na = missingstring)
  names(length_weight) <- tolower(names(length_weight))
  length_weight <- length_weight %>%
    filter(!is.na(measurement_type)) %>%
    mutate(measurement_value = suppressWarnings(as.numeric(measurement_value))) %>%
    pivot_wider(
      id_cols = c(specimen_id, species_code, event_id),
      names_from = measurement_type,
      values_from = measurement_value,
      values_fn = mean
    ) %>%
    mutate(
      fork_length = as.numeric(fork_length),
      organism_weight = as.numeric(organism_weight)
    ) %>%
    drop_na()

  readr::write_csv(scaling, file.path(surveydir, "scaling.csv"))
  readr::write_csv(acoustics, file.path(surveydir, "acoustics_projected.csv"))
  readr::write_csv(length_weight, file.path(surveydir, "length_weight.csv"))
  readr::write_csv(trawl_locations, file.path(surveydir, "trawl_locations_projected.csv"))
  readr::write_csv(surveygrid, file.path(surveydir, "surveygrid.csv"))

  invisible(
    list(
      scaling = scaling,
      acoustics = acoustics,
      length_weight = length_weight,
      trawl_locations = trawl_locations,
      surveygrid = surveygrid,
      surveyhull = surveyhull
    )
  )
}
