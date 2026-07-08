# Central path configuration for the pipeline.
# Edit these two roots to match your machine; every script sources this file.

# Repo root (this file lives in <repo>/config/)
base_repo <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."), mustWork = FALSE)
if (is.na(base_repo) || base_repo == "..") base_repo <- getwd()

# Where the discovery Seurat object and intermediate objects live.
# By default kept OUTSIDE the repo (large; git-ignored). Point this wherever you keep them.
data_dir     <- file.path(base_repo, "data")
# raw_dir / external_dir default to the repo's data/, but can be pointed elsewhere
# (e.g., a data drive) via environment variables without editing this file.
raw_dir      <- Sys.getenv("BPD_RAW_DIR",      file.path(data_dir, "raw"))       # this study: h5ad, Seurat_object_scrublet.RData
external_dir <- Sys.getenv("BPD_EXTERNAL_DIR", file.path(data_dir, "external"))  # GEO / LungMAP downloads
results_dir  <- file.path(base_repo, "results")
fig_dir      <- file.path(results_dir, "figures")
tab_dir      <- file.path(results_dir, "tables")

for (d in c(fig_dir, tab_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# Key intermediate objects (produced by upstream scripts)
seurat_obj_file <- file.path(raw_dir, "Seurat_object_scrublet.RData")   # <- 02 output
deseq_xlsx      <- file.path(tab_dir, "condition_deg_by_celltype_deseq.xlsx")  # <- 05 output

# MSigDB / gene-set files (download or via msigdbr)
seed <- 1234
set.seed(seed)
