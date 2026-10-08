# 01 - decontX, QC filtering, normalisation and clustering.
#
# Ported from code/seurat_processing_v2.R in the analysis repository.
#
# Input   filtered_feature_bc_matrix.counts_full.h5ad  (scanpy + scrublet;
#         see scripts/aux_qc_scrublet.py and GEO accession GSE346853)
# Output  data/raw/Seurat_object_scrublet.RData
#
# PROVENANCE ONLY - NOT ON THE CRITICAL PATH ---------------------------------
# This script documents how the frozen analysis object was built. Reproduction
# of the manuscript starts from that object, not from here, because
# SCTransform, PCA, UMAP and Louvain clustering are not stable across package
# versions: a re-run can renumber the 14 clusters, and every downstream script
# relies on the cluster -> cell type mapping in config/paths.R being fixed.
#
# Running this script will OVERWRITE the frozen object. It refuses to do so
# unless BPD_ALLOW_PREPROCESS=1 is set, so that it cannot be triggered by
# accident while reproducing figures.
#
# If you do re-run it, expect the cluster numbering to differ and verify the
# annotation against marker expression before trusting anything downstream.

source(file.path("config", "paths.R"))
suppressMessages({
  library(celda)
  library(Seurat)
  library(SingleCellExperiment)
  library(zellkonverter)
  library(ggplot2)
  library(future)
})

if (!identical(Sys.getenv("BPD_ALLOW_PREPROCESS"), "1")) {
  stop("01_preprocess.R would overwrite the frozen analysis object:\n  ",
       seurat_obj_file,
       "\nReproduction starts from that object; see the header of this file.",
       "\nSet BPD_ALLOW_PREPROCESS=1 if you really intend to rebuild it.",
       call. = FALSE)
}

h5ad <- file.path(raw_dir, "filtered_feature_bc_matrix.counts_full.h5ad")
if (!file.exists(h5ad)) {
  stop("Raw counts not found:\n  ", h5ad,
       "\nDeposited at GEO accession GSE346853.", call. = FALSE)
}

# ---- decontX ambient-RNA correction ----------------------------------------
sce <- readH5AD(h5ad)
sce <- decontX(sce)
seu <- CreateSeuratObject(counts = decontXcounts(sce), project = "bpd")

# ---- normalise and compute QC metrics --------------------------------------
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst")
seu <- PercentageFeatureSet(seu, pattern = "^mt-", col.name = "percent.mt")
seu@meta.data$orig.ident <- vapply(
  colnames(seu), function(x) strsplit(x, "-", fixed = TRUE)[[1]][2],
  character(1)
)

# ---- QC filter: < 20,000 UMIs and < 10% mitochondrial reads ----------------
keep <- which(seu@meta.data$nCount_RNA < 20000 & seu@meta.data$percent.mt < 10)
message(sprintf("  QC: retaining %d of %d cells", length(keep), ncol(seu)))
seu <- seu[, keep]

# ---- cell-cycle scoring (human gene lists mapped to mouse symbols) ---------
to_mouse <- function(x) {
  x <- strsplit(tolower(x), "", fixed = TRUE)[[1]]
  paste(c(toupper(x[1]), x[-1]), collapse = "")
}
s_genes <- vapply(cc.genes$s.genes, to_mouse, character(1))
g2m_genes <- vapply(cc.genes$g2m.genes, to_mouse, character(1))
seu <- CellCycleScoring(seu,
  s.features = s_genes, g2m.features = g2m_genes, set.ident = TRUE
)

# ---- SCTransform, regressing out cell cycle and mitochondrial content ------
plan(sequential)
options(future.globals.maxSize = Inf)
seu <- SCTransform(seu,
  vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), verbose = FALSE
)

# ---- library -> sample, sex and exposure -----------------------------------
conversions <- c("1" = "94", "2" = "93", "3" = "92", "4" = "95")
seu@meta.data$sample <- conversions[seu@meta.data$orig.ident]
seu@meta.data$sex <- ifelse(
  seu@meta.data$sample %in% c("93", "95"), "male", "female"
)
seu@meta.data$condition <- ifelse(
  seu@meta.data$sample %in% c("94", "95"), "O2", "Air"
)

# ---- dimension reduction and clustering ------------------------------------
seu <- RunPCA(seu, verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:20, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.2, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:20, verbose = FALSE)

message(sprintf("  clusters at resolution 0.2: %d",
                length(unique(seu$seurat_clusters))))
if (length(unique(seu$seurat_clusters)) != 14) {
  warning("This run produced a different number of clusters than the frozen ",
          "object (14). The cluster -> cell type map in config/paths.R will ",
          "NOT be valid for this object.", call. = FALSE, immediate. = TRUE)
}

save(seu, file = seurat_obj_file)
message("01_preprocess.R complete -> ", seurat_obj_file)
