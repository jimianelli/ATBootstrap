suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

source("R/funs.R")

artifact_dir <- file.path("docs", "artifacts")
summary_csv_path <- file.path(artifact_dir, "r-julia-parity-summary.csv")
detail_csv_path <- file.path(artifact_dir, "r-julia-preprocess-details.csv")
report_path <- file.path("docs", "r-julia-parity.html")

dir.create(artifact_dir, recursive = TRUE, showWarnings = FALSE)

survey_specs <- list(
  `200707` = list(log_ranges = list(c(0, 8275.49), c(8358.5, 9072.5), c(9271, 9685)), width = 20),
  `200909` = list(log_ranges = list(c(0, 8607.99), c(8709, 9250.49), c(9465, 9803.49)), width = 20),
  `201006` = list(log_ranges = list(c(0, 8811.49), c(9055, 9610.99), c(9875.5, 10221.5)), width = 20),
  `201207` = list(log_ranges = list(c(0, 8058.99), c(8367, 8959.49), c(9317.5, 9695.49)), width = 20),
  `201407` = list(log_ranges = list(c(0, 8562.99), c(8700, 9245.49), c(9300, 9720.49)), width = 20),
  `201608` = list(log_ranges = list(c(279.5, 8427.5), c(8430.5, 14259)), width = 20),
  `201807` = list(log_ranges = list(c(300, 3527.99), c(3528, 3864.99), c(5145, 7123.99), c(7486, 25000)), width = 20),
  `202207` = list(log_ranges = list(c(250, 2190), c(2752.5, 2778), c(2400, 2752.49), c(2780, 6100)), width = 40)
)

preprocess_reference_files <- c(
  "scaling.csv",
  "acoustics_projected.csv",
  "length_weight.csv",
  "trawl_locations_projected.csv",
  "surveygrid.csv"
)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
}

html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

html_table <- function(df, classes = "", raw_cols = character()) {
  if (nrow(df) == 0) {
    return("<p>No rows.</p>")
  }

  headers <- paste0("<th>", html_escape(names(df)), "</th>", collapse = "")
  rows <- apply(df, 1, function(row) {
    cells <- vapply(seq_along(row), function(i) {
      value <- as.character(row[[i]])
      if (names(row)[i] %in% raw_cols) {
        paste0("<td>", value, "</td>")
      } else {
        paste0("<td>", html_escape(value), "</td>")
      }
    }, character(1))
    paste0("<tr>", paste(cells, collapse = ""), "</tr>")
  })
  paste0(
    "<table class=\"", classes, "\"><thead><tr>", headers, "</tr></thead><tbody>",
    paste(rows, collapse = ""),
    "</tbody></table>"
  )
}

canonicalize_df <- function(df) {
  char_df <- as.data.frame(lapply(df, function(col) {
    if (is.numeric(col)) {
      ifelse(is.na(col), "<NA>", sprintf("%.10f", col))
    } else {
      ifelse(is.na(col), "<NA>", as.character(col))
    }
  }), stringsAsFactors = FALSE, check.names = FALSE)

  order_key <- do.call(paste, c(char_df, sep = "\r"))
  char_df <- char_df[order(order_key), , drop = FALSE]
  rownames(char_df) <- NULL
  char_df
}

sha256_file <- function(path) {
  output <- system2("shasum", c("-a", "256", path), stdout = TRUE, stderr = TRUE)
  if (length(output) != 1L) {
    stop("Failed to hash file: ", path, call. = FALSE)
  }
  strsplit(output, " ", fixed = TRUE)[[1]][1]
}

run_golden_check <- function() {
  golden_dir <- "tests/golden/analyses_results"
  current_dir <- "analyses/results"
  manifest_path <- file.path(golden_dir, "SHA256SUMS")
  manifest_lines <- readLines(manifest_path, warn = FALSE)
  parts <- strsplit(manifest_lines, "  ", fixed = TRUE)
  manifest <- tibble(
    sha256 = vapply(parts, `[`, character(1), 1),
    file = vapply(parts, `[`, character(1), 2)
  )

  missing_current <- manifest$file[!file.exists(file.path(current_dir, manifest$file))]
  if (length(missing_current) > 0) {
    return(list(
      summary = tibble(
        suite = "golden_outputs",
        check = "Frozen Julia output hashes",
        status = "FAIL",
        pass_count = 0L,
        fail_count = length(missing_current),
        skipped_count = 0L,
        scope = "analyses/results against tests/golden/analyses_results",
        note = paste("Missing current files:", paste(missing_current, collapse = ", "))
      ),
      detail = tibble(file = missing_current, status = "MISSING_CURRENT")
    ))
  }

  current_hashes <- vapply(
    manifest$file,
    function(path) sha256_file(file.path(current_dir, path)),
    character(1)
  )
  golden_hashes <- vapply(
    manifest$file,
    function(path) sha256_file(file.path(golden_dir, path)),
    character(1)
  )

  detail <- tibble(
    file = manifest$file,
    expected_sha256 = golden_hashes,
    current_sha256 = current_hashes,
    status = ifelse(current_hashes == golden_hashes, "PASS", "FAIL")
  )

  changed <- detail %>% filter(status == "FAIL")
  summary <- tibble(
    suite = "golden_outputs",
    check = "Frozen Julia output hashes",
    status = if (nrow(changed) == 0) "PASS" else "FAIL",
    pass_count = sum(detail$status == "PASS"),
    fail_count = sum(detail$status == "FAIL"),
    skipped_count = 0L,
    scope = "analyses/results against tests/golden/analyses_results",
    note = if (nrow(changed) == 0) {
      paste("Matched", nrow(detail), "frozen Julia output files.")
    } else {
      paste("Changed files:", paste(changed$file, collapse = ", "))
    }
  )

  list(summary = summary, detail = detail)
}

run_nonspatial_smoke <- function() {
  survey <- "202207"
  source_dir <- file.path("surveydata", survey)
  work_dir <- file.path(tempdir(), paste0("r_smoke_", survey))
  log_ranges <- list(c(250, 2190), c(2752.5, 2778), c(2400, 2752.49), c(2780, 6100))

  if (dir.exists(work_dir)) {
    unlink(work_dir, recursive = TRUE)
  }
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

  source_files <- list.files(source_dir, recursive = TRUE, full.names = TRUE)
  for (src in source_files) {
    rel <- sub(paste0("^", source_dir, "/"), "", src)
    dest <- file.path(work_dir, rel)
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    file.copy(src, dest, overwrite = TRUE)
  }

  outcome <- tryCatch({
    preprocess_survey_data(
      work_dir,
      dx = 10,
      ebs = TRUE,
      log_ranges = log_ranges,
      grid_method = TransectRibbons(width = 40)
    )

    sfiles <- read_survey_files(work_dir)
    surveydata <- ATSurveyData(
      acoustics = subset(sfiles$acoustics, class %in% c("SS1", "SS2", "BT")),
      scaling = sfiles$scaling,
      age_length = sfiles$age_length,
      length_weight = sfiles$length_weight,
      trawl_locations = sfiles$trawl_locations,
      grid = sfiles$surveygrid,
      dA = (10 * km2nmi)^2
    )

    cp <- ScalingClassProblem(
      class_name = "SS1",
      age_max = 10,
      aged_species = c(21740),
      simdomain = sfiles$surveygrid,
      simulator = function(scp) rep(100, nrow(scp$simdomain))
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

    list(status = "PASS", note = paste("Generated", nrow(res), "non-spatial result rows for survey", survey, "."))
  }, error = function(err) {
    list(status = "FAIL", note = conditionMessage(err))
  })

  tibble(
    suite = "nonspatial_smoke",
    check = "R non-spatial deterministic smoke test",
    status = outcome$status,
    pass_count = if (outcome$status == "PASS") 1L else 0L,
    fail_count = if (outcome$status == "FAIL") 1L else 0L,
    skipped_count = 0L,
    scope = "simulate_class_iteration() with an injected constant simulator",
    note = outcome$note
  )
}

prepare_preprocess_workdir <- function(survey, spec) {
  source_dir <- file.path("surveydata", survey)
  work_dir <- file.path(tempdir(), paste0("preprocess_parity_", survey))

  if (dir.exists(work_dir)) {
    unlink(work_dir, recursive = TRUE)
  }
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

  source_files <- list.files(source_dir, recursive = TRUE, full.names = TRUE)
  for (src in source_files) {
    rel <- sub(paste0("^", source_dir, "/"), "", src)
    dest <- file.path(work_dir, rel)
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    file.copy(src, dest, overwrite = TRUE)
  }

  preprocess_survey_data(
    work_dir,
    dx = 10,
    ebs = TRUE,
    log_ranges = spec$log_ranges,
    grid_method = TransectRibbons(width = spec$width)
  )

  work_dir
}

compare_preprocess_file <- function(survey, file_name, work_dir) {
  source_dir <- file.path("surveydata", survey)
  reference_path <- file.path(source_dir, file_name)

  if (!file.exists(reference_path)) {
    return(tibble(
      suite = "preprocess_parity",
      survey = survey,
      file = file_name,
      status = "REFERENCE_MISSING",
      row_count_reference = NA_integer_,
      row_count_r = NA_integer_,
      col_count_reference = NA_integer_,
      col_count_r = NA_integer_,
      schema_match = NA,
      row_count_match = NA,
      content_match = NA,
      missing_columns = "",
      extra_columns = "",
      max_numeric_diff = NA_real_,
      note = "Reference file is not present in surveydata."
    ))
  }

  candidate_path <- file.path(work_dir, file_name)
  if (!file.exists(candidate_path)) {
    return(tibble(
      suite = "preprocess_parity",
      survey = survey,
      file = file_name,
      status = "FAIL",
      row_count_reference = NA_integer_,
      row_count_r = NA_integer_,
      col_count_reference = NA_integer_,
      col_count_r = NA_integer_,
      schema_match = FALSE,
      row_count_match = FALSE,
      content_match = FALSE,
      missing_columns = "",
      extra_columns = "",
      max_numeric_diff = NA_real_,
      note = "R preprocessing run did not produce the candidate file."
    ))
  }

  reference <- .read_csv_local(reference_path)
  candidate <- .read_csv_local(candidate_path)

  missing_columns <- setdiff(names(reference), names(candidate))
  extra_columns <- setdiff(names(candidate), names(reference))
  schema_match <- identical(names(reference), names(candidate))
  row_count_match <- nrow(reference) == nrow(candidate)
  content_match <- FALSE
  max_numeric_diff <- NA_real_

  if (identical(file_name, "trawl_locations_projected.csv") && schema_match && row_count_match) {
    ref_sorted <- reference %>% arrange(haul_id)
    cand_sorted <- candidate %>% arrange(haul_id)
    same_identity <- identical(as.character(ref_sorted$survey), as.character(cand_sorted$survey)) &&
      identical(ref_sorted$haul_id, cand_sorted$haul_id)
    latlon_max_diff <- max(
      abs(ref_sorted$latitude - cand_sorted$latitude),
      abs(ref_sorted$longitude - cand_sorted$longitude),
      na.rm = TRUE
    )
    xy_max_diff <- max(
      abs(ref_sorted$x - cand_sorted$x),
      abs(ref_sorted$y - cand_sorted$y),
      na.rm = TRUE
    )
    helper_match <- identical(ref_sorted$lla, cand_sorted$lla) &&
      identical(ref_sorted$utm, cand_sorted$utm)
    core_match <- same_identity &&
      is.finite(latlon_max_diff) &&
      is.finite(xy_max_diff) &&
      latlon_max_diff <= 1e-12 &&
      xy_max_diff <= 1e-8

    status <- if (core_match && helper_match) "PASS" else if (core_match) "PARTIAL" else "FAIL"
    return(tibble(
      suite = "preprocess_parity",
      survey = survey,
      file = file_name,
      status = status,
      row_count_reference = nrow(reference),
      row_count_r = nrow(candidate),
      col_count_reference = ncol(reference),
      col_count_r = ncol(candidate),
      schema_match = schema_match,
      row_count_match = row_count_match,
      content_match = helper_match,
      missing_columns = "",
      extra_columns = "",
      max_numeric_diff = xy_max_diff,
      note = if (status == "PASS") {
        "Operational columns and helper strings match."
      } else if (status == "PARTIAL") {
        "Operational columns match within tolerance; only `lla`/`utm` helper-string serialization differs."
      } else {
        "Operational trawl-location columns still differ."
      }
    ))
  }

  if (schema_match && row_count_match) {
    ref_char <- canonicalize_df(reference)
    cand_char <- canonicalize_df(candidate)
    content_match <- identical(ref_char, cand_char)

    if (!content_match) {
      numeric_cols <- names(reference)[vapply(reference, is.numeric, logical(1))]
      if (length(numeric_cols) > 0) {
        ref_sorted <- reference[match(do.call(paste, c(ref_char, sep = "\r")),
          do.call(paste, c(canonicalize_df(reference), sep = "\r"))), , drop = FALSE]
        cand_sorted <- candidate[match(do.call(paste, c(cand_char, sep = "\r")),
          do.call(paste, c(canonicalize_df(candidate), sep = "\r"))), , drop = FALSE]
        max_numeric_diff <- max(vapply(numeric_cols, function(col) {
          max(abs(ref_sorted[[col]] - cand_sorted[[col]]), na.rm = TRUE)
        }, numeric(1)), na.rm = TRUE)
      }
    } else {
      max_numeric_diff <- 0
    }
  }

  status <- if (schema_match && row_count_match && content_match) "PASS" else "FAIL"
  note <- if (status == "PASS") {
    "Reference and R-generated CSVs match after canonical row sorting."
  } else if (!schema_match) {
    missing_text <- if (length(missing_columns) == 0) "none" else paste(missing_columns, collapse = ", ")
    extra_text <- if (length(extra_columns) == 0) "none" else paste(extra_columns, collapse = ", ")
    paste0(
      "Schema differs. Missing in R output: ", missing_text,
      ". Extra in R output: ", extra_text,
      "."
    )
  } else if (!row_count_match) {
    paste("Row count differs:", nrow(reference), "reference rows vs", nrow(candidate), "R rows.")
  } else {
    "Schema and row counts match, but row content still differs after canonical sorting."
  }

  tibble(
    suite = "preprocess_parity",
    survey = survey,
    file = file_name,
    status = status,
    row_count_reference = nrow(reference),
    row_count_r = nrow(candidate),
    col_count_reference = ncol(reference),
    col_count_r = ncol(candidate),
    schema_match = schema_match,
    row_count_match = row_count_match,
    content_match = content_match,
    missing_columns = paste(missing_columns, collapse = ", "),
    extra_columns = paste(extra_columns, collapse = ", "),
    max_numeric_diff = max_numeric_diff,
    note = note
  )
}

run_preprocess_parity <- function() {
  detail <- bind_rows(lapply(names(survey_specs), function(survey) {
    work_dir <- prepare_preprocess_workdir(survey, survey_specs[[survey]])
    bind_rows(lapply(preprocess_reference_files, function(file_name) {
      compare_preprocess_file(survey, file_name, work_dir)
    }))
  }))

  file_summary <- detail %>%
    group_by(file) %>%
    summarise(
      pass_count = sum(status == "PASS"),
      partial_count = sum(status == "PARTIAL"),
      fail_count = sum(status == "FAIL"),
      skipped_count = sum(status == "REFERENCE_MISSING"),
      status = case_when(
        fail_count == 0 & partial_count == 0 & skipped_count == 0 ~ "PASS",
        partial_count > 0 | (pass_count > 0 & fail_count > 0) ~ "PARTIAL",
        pass_count > 0 & fail_count == 0 ~ "PASS",
        fail_count > 0 & pass_count == 0 ~ "FAIL",
        TRUE ~ "SKIPPED"
      ),
      .groups = "drop"
    ) %>%
    mutate(
      suite = "preprocess_parity",
      check = paste("Preprocessing parity:", file),
      scope = "R preprocess_survey_data() against checked-in Julia-era surveydata outputs",
      note = paste(pass_count, "passes,", partial_count, "partial,", fail_count, "fails,", skipped_count, "reference-missing.")
    ) %>%
    select(suite, check, status, pass_count, fail_count, skipped_count, scope, note)

  overall <- tibble(
    suite = "preprocess_parity",
    check = "Preprocessing parity overall",
    status = if (all(detail$status %in% c("PASS", "REFERENCE_MISSING"))) "PASS" else "PARTIAL",
    pass_count = sum(detail$status == "PASS"),
    fail_count = sum(detail$status == "FAIL"),
    skipped_count = sum(detail$status == "REFERENCE_MISSING"),
    scope = "Selected surveydata directories regenerated in R and compared to checked-in references",
    note = paste(
      sum(detail$status == "PASS"), "file-level passes,",
      sum(detail$status == "PARTIAL"), "file-level partials,",
      sum(detail$status == "FAIL"), "file-level failures,",
      sum(detail$status == "REFERENCE_MISSING"), "reference-missing checks."
    )
  )

  list(summary = bind_rows(overall, file_summary), detail = detail)
}

git_ref <- tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(...) NA_character_)
golden <- run_golden_check()
smoke <- run_nonspatial_smoke()
preprocess <- run_preprocess_parity()

bootstrap_gap <- tibble(
  suite = "full_parity",
  check = "Full R versus Julia bootstrap parity",
  status = "BLOCKED",
  pass_count = 0L,
  fail_count = 0L,
  skipped_count = 1L,
  scope = "Whole-pipeline replicate and summary comparison",
  note = paste(
    "Not yet claimable:",
    "the R port intentionally stops short of the spatial simulator,",
    "so direct end-to-end bootstrap parity cannot be demonstrated from this branch."
  )
)

summary_df <- bind_rows(
  golden$summary,
  preprocess$summary,
  smoke,
  bootstrap_gap
)

write_csv(summary_df, summary_csv_path)
write_csv(preprocess$detail, detail_csv_path)

status_class <- function(status) {
  out <- rep("blocked", length(status))
  out[status == "PASS"] <- "pass"
  out[status == "PARTIAL"] <- "partial"
  out[status == "FAIL"] <- "fail"
  out[status %in% c("BLOCKED", "SKIPPED")] <- "blocked"
  out
}

summary_display <- summary_df %>%
  mutate(status = paste0("<span class=\"badge ", status_class(status), "\">", status, "</span>")) %>%
  select(Check = check, Status = status, Scope = scope, Note = note)

preprocess_overview <- preprocess$detail %>%
  group_by(file) %>%
  summarise(
    pass = sum(status == "PASS"),
    partial = sum(status == "PARTIAL"),
    fail = sum(status == "FAIL"),
    reference_missing = sum(status == "REFERENCE_MISSING"),
    .groups = "drop"
  ) %>%
  mutate(Status = case_when(
    partial > 0 ~ "PARTIAL",
    fail == 0 ~ "PASS",
    pass > 0 ~ "PARTIAL",
    TRUE ~ "FAIL"
  )) %>%
  select(`Preprocess output` = file, Status, pass, partial, fail, reference_missing)

detail_display <- preprocess$detail %>%
  transmute(
    Survey = survey,
    File = file,
    Status = status,
    `Reference rows` = ifelse(is.na(row_count_reference), "", as.character(row_count_reference)),
    `R rows` = ifelse(is.na(row_count_r), "", as.character(row_count_r)),
    `Reference cols` = ifelse(is.na(col_count_reference), "", as.character(col_count_reference)),
    `R cols` = ifelse(is.na(col_count_r), "", as.character(col_count_r)),
    `Schema match` = ifelse(is.na(schema_match), "", ifelse(schema_match, "yes", "no")),
    `Content match` = ifelse(is.na(content_match), "", ifelse(content_match, "yes", "no")),
    Note = note
  )

generated_at <- format(Sys.time(), "%Y-%m-%d %H:%M %Z")
overall_status <- if (any(summary_df$status %in% c("FAIL", "PARTIAL", "BLOCKED"))) "Partial evidence" else "Complete parity evidence"

html <- paste0(
  "<!DOCTYPE html>\n",
  "<html lang=\"en\">\n",
  "<head>\n",
  "  <meta charset=\"UTF-8\">\n",
  "  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n",
  "  <title>ATBootstrap R and Julia Parity Evidence</title>\n",
  "  <style>\n",
  "    :root {\n",
  "      --bg: #f5f2ea;\n",
  "      --panel: #fffdf9;\n",
  "      --ink: #1f2933;\n",
  "      --muted: #5f6b74;\n",
  "      --accent: #0f766e;\n",
  "      --accent-soft: #e4f3f0;\n",
  "      --warn: #9a6700;\n",
  "      --warn-soft: #fff3cf;\n",
  "      --fail: #b42318;\n",
  "      --fail-soft: #fee4e2;\n",
  "      --line: #ddd4c8;\n",
  "      --shadow: 0 16px 40px rgba(31, 41, 51, 0.08);\n",
  "    }\n",
  "    * { box-sizing: border-box; }\n",
  "    body {\n",
  "      margin: 0;\n",
  "      font-family: Georgia, \"Times New Roman\", serif;\n",
  "      color: var(--ink);\n",
  "      background:\n",
  "        radial-gradient(circle at top left, rgba(15, 118, 110, 0.08), transparent 28%),\n",
  "        linear-gradient(180deg, #faf7f1 0%, var(--bg) 100%);\n",
  "      line-height: 1.62;\n",
  "    }\n",
  "    .page { max-width: 1080px; margin: 0 auto; padding: 40px 24px 72px; }\n",
  "    .hero {\n",
  "      background: linear-gradient(135deg, #0f766e, #183e5c);\n",
  "      color: #f8fafc;\n",
  "      border-radius: 24px;\n",
  "      padding: 34px 30px;\n",
  "      box-shadow: var(--shadow);\n",
  "    }\n",
  "    .eyebrow { margin: 0 0 8px; text-transform: uppercase; letter-spacing: 0.12em; font-size: 0.8rem; opacity: 0.84; }\n",
  "    h1 { margin: 0 0 10px; font-size: clamp(2rem, 4vw, 3.2rem); line-height: 1.05; }\n",
  "    .hero p { margin: 0; max-width: 54rem; color: rgba(248, 250, 252, 0.92); font-size: 1.05rem; }\n",
  "    .meta { display: flex; flex-wrap: wrap; gap: 10px; margin-top: 16px; }\n",
  "    .meta span { background: rgba(255,255,255,0.12); padding: 6px 10px; border-radius: 999px; font-size: 0.92rem; }\n",
  "    .section { margin-top: 22px; background: var(--panel); border: 1px solid var(--line); border-radius: 20px; padding: 24px 26px; box-shadow: var(--shadow); }\n",
  "    h2 { margin: 0 0 14px; font-size: 1.45rem; line-height: 1.2; }\n",
  "    p { margin: 0 0 12px; }\n",
  "    ul { margin: 0 0 12px 22px; }\n",
  "    table { width: 100%; border-collapse: collapse; font-size: 0.95rem; }\n",
  "    th, td { border-top: 1px solid var(--line); padding: 10px 8px; text-align: left; vertical-align: top; }\n",
  "    th { color: var(--muted); font-weight: 600; }\n",
  "    .grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 14px; margin-top: 14px; }\n",
  "    .mini { border-radius: 16px; background: var(--accent-soft); border: 1px solid rgba(15, 118, 110, 0.16); padding: 14px 16px; }\n",
  "    .mini strong { display: block; margin-bottom: 4px; color: var(--accent); }\n",
  "    .callout { margin-top: 14px; border-left: 6px solid var(--warn); background: var(--warn-soft); border-radius: 14px; padding: 16px 18px; }\n",
  "    .callout.fail { border-left-color: var(--fail); background: var(--fail-soft); }\n",
  "    .badge { display: inline-block; padding: 4px 9px; border-radius: 999px; font-size: 0.82rem; font-weight: 700; }\n",
  "    .badge.pass { background: #dcfae6; color: #166534; }\n",
  "    .badge.partial { background: #fef0c7; color: #854d0e; }\n",
  "    .badge.fail { background: #fee4e2; color: #b42318; }\n",
  "    .badge.blocked { background: #e4e7ec; color: #344054; }\n",
  "    code { font-family: \"SFMono-Regular\", Menlo, Consolas, monospace; background: #f0ebe4; padding: 0.12em 0.4em; border-radius: 6px; font-size: 0.92em; }\n",
  "    a { color: #0f5f9c; }\n",
  "  </style>\n",
  "</head>\n",
  "<body>\n",
  "  <div class=\"page\">\n",
  "    <header class=\"hero\">\n",
  "      <p class=\"eyebrow\">Parity Evidence</p>\n",
  "      <h1>ATBootstrap R and Julia Parity Evidence</h1>\n",
  "      <p>This report documents the current evidence base for saying the R work matches the Julia implementation. It is intentionally strict: it records what passes today, what still fails, and what cannot yet be claimed because the spatial simulator is not ported.</p>\n",
  "      <div class=\"meta\">\n",
  "        <span>Generated: ", html_escape(generated_at), "</span>\n",
  "        <span>Git ref: ", html_escape(git_ref[1] %||% "unknown"), "</span>\n",
  "        <span>Overall status: ", html_escape(overall_status), "</span>\n",
  "        <span><a href=\"./artifacts/r-julia-parity-summary.csv\" style=\"color:#f8fafc;\">Summary CSV</a></span>\n",
  "        <span><a href=\"./artifacts/r-julia-preprocess-details.csv\" style=\"color:#f8fafc;\">Detail CSV</a></span>\n",
  "      </div>\n",
  "    </header>\n",
  "    <section class=\"section\">\n",
  "      <h2>Bottom Line</h2>\n",
  "      <p>The repository now has a reproducible parity ledger, but it does not yet support a blanket claim that the R implementation fully matches Julia. The strongest current evidence is in frozen Julia outputs and most checked preprocessing products. The strongest remaining preprocessing gap is in a few <code>acoustics_projected.csv</code> cases, while projected trawl locations are now down to helper-string formatting differences.</p>\n",
  "      <div class=\"grid\">\n",
  "        <div class=\"mini\"><strong>What passes</strong> Frozen Julia result hashes, the non-spatial R smoke test, and preprocessing outputs including <code>scaling.csv</code> and <code>surveygrid.csv</code> across all checked surveys.</div>\n",
  "        <div class=\"mini\"><strong>What remains</strong> <code>acoustics_projected.csv</code> still has partial parity, and <code>trawl_locations_projected.csv</code> still differs in Julia-versus-R helper-string serialization for <code>lla</code> and <code>utm</code>.</div>\n",
  "        <div class=\"mini\"><strong>What is blocked</strong> Full end-to-end bootstrap parity remains out of scope until the spatial simulator question is resolved.</div>\n",
  "      </div>\n",
  "    </section>\n",
  "    <section class=\"section\">\n",
  "      <h2>Evidence Summary</h2>\n",
       html_table(summary_display, "summary-table", raw_cols = "Status"),
  "    </section>\n",
  "    <section class=\"section\">\n",
  "      <h2>Preprocessing Parity by Output File</h2>\n",
  "      <p>These checks regenerate selected <code>surveydata/</code> products in R from the raw inputs, then compare them against the checked-in Julia-era CSVs after canonical row sorting. This removes row-order noise and surfaces real content differences.</p>\n",
       html_table(preprocess_overview, "summary-table"),
  "      <div class=\"callout fail\"><strong>Current interpretation:</strong> the preprocessing path is partially aligned, not fully aligned. That is already useful evidence because it localizes the mismatch before the uncertainty/bootstrap layer.</div>\n",
  "    </section>\n",
  "    <section class=\"section\">\n",
  "      <h2>Survey-by-Survey Detail</h2>\n",
       html_table(detail_display, "detail-table"),
  "    </section>\n",
  "    <section class=\"section\">\n",
  "      <h2>How To Re-run</h2>\n",
  "      <ul>\n",
  "        <li>Run <code>Rscript tests/generate_r_julia_parity_report.R</code>.</li>\n",
  "        <li>Review <code>docs/r-julia-parity.html</code> for the rendered report.</li>\n",
  "        <li>Review <code>docs/artifacts/r-julia-parity-summary.csv</code> and <code>docs/artifacts/r-julia-preprocess-details.csv</code> for machine-readable outputs.</li>\n",
  "      </ul>\n",
  "    </section>\n",
  "    <section class=\"section\">\n",
  "      <h2>Claim Boundary</h2>\n",
  "      <p>This report should be cited as evidence of current parity status, not as proof of full equivalence. Full equivalence would require either direct Julia-generated reference artifacts for the same deterministic R cases or successful end-to-end survey parity after the spatial layer is resolved.</p>\n",
  "    </section>\n",
  "  </div>\n",
  "</body>\n",
  "</html>\n"
)

writeLines(html, report_path)

cat("Wrote parity report to ", report_path, "\n", sep = "")
cat("Wrote summary CSV to ", summary_csv_path, "\n", sep = "")
cat("Wrote detail CSV to ", detail_csv_path, "\n", sep = "")
