# Figure and table output helpers.
#
# Why this file exists: re-running the pipeline will never produce
# byte-identical PDFs (fonts, embedded timestamps and graphics-device
# versions all differ), so "reproducible" cannot be checked on the files
# themselves. Every panel therefore writes BOTH
#
#   results/figures/<panel>.pdf        the figure
#   results/figures/<panel>_data.csv   the numbers that were plotted
#
# scripts/99_verify.R diffs the CSVs (and Tables S1-S6) against the committed
# reference copies, which is a check that can actually pass or fail.

suppressMessages({
  library(ggplot2)
})

#' Save one manuscript panel plus the data behind it.
#'
#' @param panel Manuscript panel id, e.g. "Fig1C" or "FigS3A". Becomes the
#'   filename stem, so every output maps 1:1 onto a panel in the paper.
#' @param plot  A ggplot, pheatmap or grob.
#' @param data  Data frame of the plotted values. Omit only for panels with no
#'   tabular representation (UMAP embeddings pass the coordinates).
#' @param width,height Device size in inches.
save_panel <- function(panel, plot, data = NULL, width = 6, height = 5,
                       slug = NULL) {
  stem <- if (is.null(slug)) panel else paste0(panel, "_", slug)
  pdf_path <- file.path(fig_dir, paste0(stem, ".pdf"))

  if (inherits(plot, "pheatmap")) {
    ggplot2::ggsave(pdf_path, plot$gtable, width = width, height = height)
  } else {
    ggplot2::ggsave(pdf_path, plot, width = width, height = height)
  }

  if (!is.null(data)) {
    csv_path <- file.path(fig_dir, paste0(stem, "_data.csv"))
    utils::write.csv(round_numeric(as.data.frame(data)), csv_path,
                     row.names = FALSE)
  }
  message("  ", panel, " -> ", basename(pdf_path),
          if (!is.null(data)) " (+ _data.csv)" else "")
  invisible(pdf_path)
}

#' Round numeric columns so verification is not defeated by float noise.
round_numeric <- function(df, digits = 6) {
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], round, digits = digits)
  df
}

#' Save a supplemental table.
save_table <- function(table_id, data, path = NULL) {
  p <- if (is.null(path)) {
    file.path(tab_dir, paste0(table_id, ".csv"))
  } else {
    path
  }
  utils::write.csv(round_numeric(as.data.frame(data)), p, row.names = FALSE)
  message("  ", table_id, " -> ", basename(p), " (", nrow(data), " rows)")
  invisible(p)
}

#' Assert a number reported in the manuscript, and fail loudly if it drifts.
#'
#' This is the guard that stops a package upgrade from silently producing a
#' different paper.
#' @param on_fail "stop" for structural invariants (wrong object loaded, wrong
#'   cluster count) where continuing is meaningless; "warn" for values that can
#'   legitimately shift with an external resource version, such as MSigDB
#'   gene-set membership.
expect_value <- function(label, observed, expected, tol = 0,
                         on_fail = c("stop", "warn")) {
  on_fail <- match.arg(on_fail)
  ok <- if (is.character(expected)) {
    identical(as.character(observed), expected)
  } else {
    isTRUE(abs(observed - expected) <= tol)
  }
  if (!ok) {
    msg <- sprintf("MANUSCRIPT VALUE CHANGED: %s\n  expected: %s\n  observed: %s",
                   label, paste(expected, collapse = ", "),
                   paste(observed, collapse = ", "))
    if (on_fail == "stop") stop(msg, call. = FALSE)
    warning(msg, call. = FALSE, immediate. = TRUE)
    return(invisible(FALSE))
  }
  message(sprintf("  [ok] %s = %s", label, paste(observed, collapse = ", ")))
  invisible(TRUE)
}

#' Load the frozen object and attach cell-type annotations.
#'
#' Returns the object rather than assigning into the caller's environment, so
#' that `seu` is always an explicit argument downstream.
load_annotated <- function() {
  e <- new.env(parent = emptyenv())
  load(seurat_obj_file, envir = e)
  nm <- ls(e)
  if (length(nm) != 1L) {
    stop("Expected exactly one object in ", basename(seurat_obj_file),
         "; found: ", paste(nm, collapse = ", "))
  }
  seu <- e[[nm]]
  seu@meta.data$celltype <- unname(
    cluster_annotations[as.character(seu$seurat_clusters)]
  )
  if (anyNA(seu@meta.data$celltype)) {
    stop("Unannotated clusters: ",
         paste(sort(unique(as.character(
           seu$seurat_clusters[is.na(seu@meta.data$celltype)]
         ))), collapse = ", "),
         "\nThe frozen object does not match config/paths.R.")
  }
  # objects built before the gender -> sex rename carry `gender`
  if (is.null(seu@meta.data$sex) && !is.null(seu@meta.data$gender)) {
    seu@meta.data$sex <- seu@meta.data$gender
  }
  expect_value("cells in atlas", ncol(seu), 32780)
  expect_value("cell types", length(unique(seu@meta.data$celltype)), 9)
  seu
}
