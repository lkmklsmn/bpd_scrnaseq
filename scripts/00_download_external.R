# 00 - Download and verify the data this pipeline needs.
#
# Reproduction starts from a frozen Seurat object plus four external datasets
# (see data/README.md). They total 8.2 GB on disk, 5.28 GB as a gzipped
# archive, and are distributed as a
# single archive rather than fetched from each original source, because the
# LungMAP and GEO objects were processed before deposition and the exact
# intermediate files cannot be regenerated from the raw accessions.
#
# Usage
#   Rscript scripts/00_download_external.R           download, then verify
#   Rscript scripts/00_download_external.R --verify  verify what is on disk
#
# Set BPD_DATA_DIR to place the data outside the repository, which is
# recommended if the repository lives in a synced folder such as OneDrive.

source(file.path("config", "paths.R"))

# ---------------------------------------------------------------------------
# Deposited at Zenodo: https://doi.org/10.5281/zenodo.23249196
#
# Zenodo serves stable direct downloads of the form
#   https://zenodo.org/records/<RECORD_ID>/files/<FILENAME>?download=1
# which need no interstitial handling. Set archive_url to the archive's own
# such URL once the record is public, or override it with BPD_ARCHIVE_URL.
# ---------------------------------------------------------------------------
zenodo_doi <- "10.5281/zenodo.23249196"
zenodo_record <- "23249196"
archive_url <- Sys.getenv("BPD_ARCHIVE_URL", unset = "")
archive_sha256 <- "f4721ced89371c5c156774366e3b8310533e0bed751840aba627a935e45d318c"
archive_name <- "bpd_scrnaseq_data.tar.gz"

# Files the pipeline requires, with sizes and checksums for verification.
# sha256 values are filled in by --checksum below once the archive is built.
manifest <- data.frame(
  path = c(
    file.path("raw", "Seurat_object_scrublet.RData"),
    file.path("external", "GSE151974", "Seurat_object_GSE151974.RData"),
    file.path("external", "lungmap", "LMEX0000004400",
              "BPD-adata_combined.h5ad"),
    file.path("external", "lungmap", "LMEX0000004400",
              "BPD_RNA_author-clusters.txt"),
    file.path("external", "GSE32472", "gse32472.RData")
  ),
  bytes = c(1600926947, 2577472660, 4522574254, 16204931, 80364579),
  sha256 = c(
    "a7b2aeae7c9bf628d74fb2386bad845093739e9f6bc96571ce443bc8daeaba71",
    "b14fe3af379814ce4ba53c262d8344f44209dfbecaeddf01ff104b23287f9696",
    "e1c21358c293f69ade480a038af214c8c1a511ad44af700d5a868549d997ae68",
    "90f08675a1941492f375d97036b9970f545409efe2afcc645e98994358e380df",
    "8cf306852320a77cecb799242f8266cab3183a3d00ed1b913a8fd48e5d1bf194"
  ),
  needed_by = c(
    "all scripts (frozen entry point)",
    "Fig 2C, 5E, 7B, S1",
    "Fig 7C",
    "Fig 7C",
    "Fig 7D1"
  ),
  stringsAsFactors = FALSE
)

args <- commandArgs(trailingOnly = TRUE)
verify_only <- "--verify" %in% args
emit_checksums <- "--checksum" %in% args

sha256_of <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    return(NA_character_)
  }
  digest::digest(path, algo = "sha256", file = TRUE)
}

# ---- emit checksums for the manifest (run once, when building the archive)
if (emit_checksums) {
  for (i in seq_len(nrow(manifest))) {
    f <- file.path(data_dir, manifest$path[i])
    if (file.exists(f)) {
      cat(sprintf('  "%s",  # %s\n', sha256_of(f), manifest$path[i]))
    } else {
      cat(sprintf('  "",  # MISSING %s\n', manifest$path[i]))
    }
  }
  quit(save = "no", status = 0)
}

# ---- download --------------------------------------------------------------
if (!verify_only) {
  present <- vapply(manifest$path, function(p) {
    file.exists(file.path(data_dir, p))
  }, logical(1))
  if (all(present)) {
    message("All data files already present under ", data_dir,
            "; skipping download. Use --verify to check them.")
  } else {
    if (!nzchar(archive_url)) {
      stop("No download URL configured.\n",
           "Download the data from https://doi.org/", zenodo_doi,
           " and place it as listed in data/README.md, then re-run with ",
           "--verify.\nAlternatively set archive_url in this script, or ",
           "export BPD_ARCHIVE_URL.", call. = FALSE)
    }
    dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
    dest <- file.path(data_dir, archive_name)
    message("Downloading 5.28 GB to ", dest, " ...")
    options(timeout = max(7200, getOption("timeout")))
    utils::download.file(archive_url, dest, mode = "wb")

    if (nzchar(archive_sha256)) {
      got <- sha256_of(dest)
      if (!identical(got, archive_sha256)) {
        stop("Archive checksum mismatch - the download is incomplete or ",
             "corrupt.\n  expected ", archive_sha256, "\n  observed ", got,
             call. = FALSE)
      }
      message("  archive checksum OK")
    } else {
      warning("No archive checksum configured; skipping that check.",
              call. = FALSE, immediate. = TRUE)
    }

    message("Extracting ...")
    utils::untar(dest, exdir = data_dir)
    unlink(dest)
  }
}

# ---- verify ----------------------------------------------------------------
cat("\n", strrep("-", 76), "\n", sep = "")
cat(sprintf("%-58s %s\n", "file", "status"))
cat(strrep("-", 76), "\n", sep = "")
ok <- TRUE
for (i in seq_len(nrow(manifest))) {
  f <- file.path(data_dir, manifest$path[i])
  if (!file.exists(f)) {
    cat(sprintf("%-58s MISSING\n", manifest$path[i]))
    ok <- FALSE
    next
  }
  size <- file.info(f)$size
  if (size != manifest$bytes[i]) {
    cat(sprintf("%-58s WRONG SIZE (%s vs %s expected)\n", manifest$path[i],
                format(size, big.mark = ","),
                format(manifest$bytes[i], big.mark = ",")))
    ok <- FALSE
    next
  }
  if (nzchar(manifest$sha256[i])) {
    got <- sha256_of(f)
    if (!is.na(got) && !identical(got, manifest$sha256[i])) {
      cat(sprintf("%-58s CHECKSUM MISMATCH\n", manifest$path[i]))
      ok <- FALSE
      next
    }
    cat(sprintf("%-58s ok (size + sha256)\n", manifest$path[i]))
  } else {
    cat(sprintf("%-58s ok (size only)\n", manifest$path[i]))
  }
}
cat(strrep("-", 76), "\n", sep = "")

# frozen inputs ship with the repository and must not come from the archive
for (p in c(file.path("frozen", "drug_target_annotation_2026-06-30.csv"),
            file.path("frozen", "GSE220135_sample_metadata.csv"))) {
  f <- file.path(data_dir, p)
  cat(sprintf("%-58s %s\n", p,
              if (file.exists(f)) "ok (version-controlled)" else "MISSING"))
  if (!file.exists(f)) ok <- FALSE
}

if (!ok) {
  cat("\nData verification FAILED. See data/README.md.\n")
  quit(save = "no", status = 1)
}
cat("\nAll required data present. Run scripts 02-12 in order.\n")
