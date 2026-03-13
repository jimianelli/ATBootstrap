suppressPackageStartupMessages(library(tools))

golden_dir <- "tests/golden/analyses_results"
current_dir <- "analyses/results"
manifest_path <- file.path(golden_dir, "SHA256SUMS")

manifest_lines <- readLines(manifest_path, warn = FALSE)
parts <- strsplit(manifest_lines, "  ", fixed = TRUE)
manifest <- data.frame(
  sha256 = vapply(parts, `[`, character(1), 1),
  file = vapply(parts, `[`, character(1), 2),
  stringsAsFactors = FALSE
)

sha256_file <- function(path) {
  output <- system2("shasum", c("-a", "256", path), stdout = TRUE, stderr = TRUE)
  if (length(output) != 1L) {
    stop("Failed to hash file: ", path, call. = FALSE)
  }
  strsplit(output, " ", fixed = TRUE)[[1]][1]
}

missing_current <- manifest$file[!file.exists(file.path(current_dir, manifest$file))]
if (length(missing_current) > 0) {
  cat("Missing current results files:\n")
  cat(paste0(" - ", missing_current, "\n"), sep = "")
  quit(status = 1)
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

changed <- manifest$file[current_hashes != golden_hashes]
if (length(changed) > 0) {
  cat("Current analyses/results files differ from frozen golden references:\n")
  cat(paste0(" - ", changed, "\n"), sep = "")
  quit(status = 1)
}

cat("Golden-reference check passed for ", length(manifest$file), " files.\n", sep = "")
