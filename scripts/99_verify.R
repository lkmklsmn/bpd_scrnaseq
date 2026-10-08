# 99 - Verify a re-run against the committed reference outputs.
#
# Why this is not a file-hash check -------------------------------------------
# Re-running the pipeline never produces byte-identical PDFs: fonts, embedded
# timestamps and graphics-device versions all differ. So verification compares
# NUMBERS, not files. Every panel writes the data it plotted next to it
# (<panel>_data.csv, see R/figure_io.R), and this script diffs those plus the
# supplemental tables against reference copies in reference/.
#
# Usage
#   Rscript scripts/99_verify.R              compare results/ to reference/
#   Rscript scripts/99_verify.R --snapshot   populate reference/ from results/
#
# Take a snapshot once, from a run you trust, and commit it. After that this
# script is a regression test: it reports any panel or table whose numbers
# have moved, and exits non-zero so it can be wired into CI.

source(file.path("config", "paths.R"))

args <- commandArgs(trailingOnly = TRUE)
snapshot <- "--snapshot" %in% args
ref_dir <- file.path(getwd(), "reference")
tol <- 1e-6

collect <- function(root) {
  f <- c(
    list.files(file.path(root, "figures"), pattern = "_data[.]csv$",
               full.names = TRUE),
    list.files(file.path(root, "tables"), pattern = "[.]csv$",
               full.names = TRUE)
  )
  stats::setNames(f, sub(paste0("^", root, "/?"), "",
                         gsub("\\\\", "/", f)))
}

# ---- snapshot mode ---------------------------------------------------------
if (snapshot) {
  src <- collect(results_dir)
  if (!length(src)) {
    stop("Nothing to snapshot: run scripts 02-12 first.", call. = FALSE)
  }
  for (rel in names(src)) {
    dest <- file.path(ref_dir, rel)
    dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)
    file.copy(src[[rel]], dest, overwrite = TRUE)
  }
  message(sprintf("snapshot: copied %d files into reference/", length(src)))
  message("Commit reference/ so future runs can be checked against it.")
  quit(save = "no", status = 0)
}

# ---- compare mode ----------------------------------------------------------
if (!dir.exists(ref_dir)) {
  stop("No reference/ directory. Create one with:\n",
       "  Rscript scripts/99_verify.R --snapshot", call. = FALSE)
}

ref <- collect(ref_dir)
cur <- collect(results_dir)

only_ref <- setdiff(names(ref), names(cur))
only_cur <- setdiff(names(cur), names(ref))
both <- intersect(names(ref), names(cur))

#' Compare two tables: same columns, same rows, numerics within tolerance.
compare_one <- function(a_path, b_path) {
  a <- utils::read.csv(a_path, stringsAsFactors = FALSE, check.names = FALSE)
  b <- utils::read.csv(b_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!identical(dim(a), dim(b))) {
    return(sprintf("dimensions %dx%d vs %dx%d",
                   nrow(a), ncol(a), nrow(b), ncol(b)))
  }
  if (!identical(names(a), names(b))) {
    return(paste("columns differ:",
                 paste(setdiff(union(names(a), names(b)),
                               intersect(names(a), names(b))),
                       collapse = ", ")))
  }
  worst <- 0
  worst_col <- NA_character_
  for (nm in names(a)) {
    if (is.numeric(a[[nm]]) && is.numeric(b[[nm]])) {
      d <- suppressWarnings(max(abs(a[[nm]] - b[[nm]]), na.rm = TRUE))
      if (is.finite(d) && d > worst) {
        worst <- d
        worst_col <- nm
      }
    } else if (!identical(as.character(a[[nm]]), as.character(b[[nm]]))) {
      n_diff <- sum(as.character(a[[nm]]) != as.character(b[[nm]]),
                    na.rm = TRUE)
      return(sprintf("column '%s' differs in %d row(s)", nm, n_diff))
    }
  }
  if (worst > tol) {
    return(sprintf("numeric drift up to %.3g in column '%s'", worst,
                   worst_col))
  }
  NA_character_
}

results <- data.frame(file = both, status = NA_character_,
                      detail = NA_character_, stringsAsFactors = FALSE)
for (i in seq_along(both)) {
  d <- compare_one(ref[[both[i]]], cur[[both[i]]])
  results$status[i] <- if (is.na(d)) "match" else "DIFFERS"
  results$detail[i] <- if (is.na(d)) "" else d
}

cat("\n", strrep("-", 72), "\n", sep = "")
cat(sprintf("%-52s %s\n", "file", "status"))
cat(strrep("-", 72), "\n", sep = "")
for (i in seq_len(nrow(results))) {
  cat(sprintf("%-52s %s%s\n", results$file[i], results$status[i],
              if (nzchar(results$detail[i])) {
                paste0("  (", results$detail[i], ")")
              } else ""))
}
for (f in only_ref) cat(sprintf("%-52s MISSING from results/\n", f))
for (f in only_cur) cat(sprintf("%-52s not in reference/\n", f))
cat(strrep("-", 72), "\n", sep = "")

n_diff <- sum(results$status == "DIFFERS")
cat(sprintf("%d compared | %d match | %d differ | %d missing | %d new\n",
            nrow(results), sum(results$status == "match"), n_diff,
            length(only_ref), length(only_cur)))

failed <- n_diff > 0 || length(only_ref) > 0
if (failed) {
  cat("\nVERIFICATION FAILED\n")
  quit(save = "no", status = 1)
}
cat("\nAll reference outputs reproduced.\n")
