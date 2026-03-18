source("R/funs.R")

survey <- Sys.getenv("ATB_SURVEY", unset = "202207")
r_out <- Sys.getenv(
  "ATB_SPATIAL_CONDITIONAL_R_OUT",
  unset = file.path("docs", "artifacts", "spatial_conditional_r.csv")
)
j_out <- Sys.getenv(
  "ATB_SPATIAL_CONDITIONAL_J_OUT",
  unset = file.path("docs", "artifacts", "spatial_conditional_julia.csv")
)
cmp_out <- Sys.getenv(
  "ATB_SPATIAL_CONDITIONAL_CMP_OUT",
  unset = file.path("docs", "artifacts", "spatial_conditional_parity.csv")
)

status <- system2("Rscript", c("tests/spatial_conditional_probe.R"), stdout = TRUE, stderr = TRUE)
cat(paste(status, collapse = "\n"), "\n")

status <- system2("julia", c("--project=.", "tests/spatial_conditional_probe.jl"), stdout = TRUE, stderr = TRUE)
cat(paste(status, collapse = "\n"), "\n")

r_tbl <- readr::read_csv(r_out, show_col_types = FALSE)
j_tbl <- readr::read_csv(j_out, show_col_types = FALSE)

shared_num <- intersect(
  names(select(r_tbl, where(is.numeric))),
  names(select(j_tbl, where(is.numeric)))
)
shared_num <- setdiff(shared_num, c("survey"))

cmp <- full_join(r_tbl, j_tbl, by = c("survey", "class"), suffix = c("_r", "_julia"))

for (col in shared_num) {
  rcol <- paste0(col, "_r")
  jcol <- paste0(col, "_julia")
  diff_col <- paste0(col, "_diff")
  rel_col <- paste0(col, "_rel_diff")
  denom <- pmax(abs(cmp[[jcol]]), 1e-8)
  cmp[[diff_col]] <- cmp[[rcol]] - cmp[[jcol]]
  cmp[[rel_col]] <- abs(cmp[[diff_col]]) / denom
}

cmp$structural_match <- with(
  cmp,
  params_dloc_n_r == params_dloc_n_julia &
    params_sloc_n_r == params_sloc_n_julia &
    params_data_n_r == params_data_n_julia &
    prep_shift_n_r == prep_shift_n_julia
)

cmp$status <- ifelse(cmp$structural_match, "NUMERIC_ONLY", "STRUCTURAL_MISMATCH")

readr::write_csv(cmp, cmp_out)
cat("Wrote spatial conditional parity metrics to ", cmp_out, "\n", sep = "")
