# 02 - decontX, QC filtering, normalization, clustering.
# Source: seurat_processing_v2.R
# Input : data/raw/filtered_feature_bc_matrix.counts_full.h5ad  (from 01)
# Output: data/raw/Seurat_object_scrublet.RData                  (-> 03-10)
# Figures: QC (Fig S3), elbow (Fig S3)

source(file.path("config", "paths.R"))
suppressMessages({
  library(celda); library(Seurat); library(SingleCellExperiment)
  library(zellkonverter); library(ggplot2); library(future)
})

# ---- load raw counts (scanpy/scrublet output) & run decontX -----------------
sce <- readH5AD(file.path(raw_dir, "filtered_feature_bc_matrix.counts_full.h5ad"))
sce <- decontX(sce)
seu <- CreateSeuratObject(counts = decontXcounts(sce), project = "bpd")

# ---- normalize + QC metrics -------------------------------------------------
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst")
seu <- PercentageFeatureSet(seu, pattern = "^mt-", col.name = "percent.mt")

seu@meta.data$orig.ident <- vapply(colnames(seu),
  function(x) strsplit(x, "-", fixed = TRUE)[[1]][2], character(1))

# QC scatter (Fig S3): thresholds < 20,000 UMIs and < 10% mito
ggplot(seu@meta.data, aes(nCount_RNA, percent.mt, color = orig.ident)) +
  facet_wrap(~ orig.ident) + geom_point() +
  geom_vline(xintercept = 20000, linetype = "dashed") +
  geom_hline(yintercept = 10, linetype = "dashed") + theme_minimal()
ggsave(file.path(fig_dir, "S3_qc_scatter.pdf"), height = 6, width = 7)

seu <- seu[, which(seu@meta.data$nCount_RNA < 20000 & seu@meta.data$percent.mt < 10)]

# ---- cell-cycle scoring (human -> mouse gene symbols) -----------------------
to_mouse <- function(x){ x <- tolower(x); x <- strsplit(x, "")[[1]]
  paste(c(toupper(x[1]), x[-1]), collapse = "") }
s.genes   <- vapply(cc.genes$s.genes,   to_mouse, character(1))
g2m.genes <- vapply(cc.genes$g2m.genes, to_mouse, character(1))
seu <- CellCycleScoring(seu, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# ---- SCTransform ------------------------------------------------------------
plan(sequential); options(future.globals.maxSize = Inf)
seu <- SCTransform(seu, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), verbose = FALSE)

# ---- sample metadata (library -> sample/sex/condition) ----------------------
conversions <- c("1" = "94", "2" = "93", "3" = "92", "4" = "95")
seu@meta.data$sample    <- conversions[seu@meta.data$orig.ident]
seu@meta.data$sex       <- ifelse(seu@meta.data$sample %in% c("93", "95"), "male", "female")
seu@meta.data$condition <- ifelse(seu@meta.data$sample %in% c("94", "95"), "O2", "Air")

# ---- dimension reduction + clustering (resolution 0.2 -> 14 clusters) -------
seu <- RunPCA(seu, verbose = FALSE)
var_prop <- seu@reductions$pca@stdev / sum(seu@reductions$pca@stdev)
ggplot(data.frame(PC = factor(seq_along(var_prop)), variance = var_prop),
       aes(PC, variance)) + geom_col(fill = "grey30") +
  ylab("Proportion of variance") + xlab("Principal component") + theme_classic()
ggsave(file.path(fig_dir, "S3_elbow_pca.pdf"), height = 5, width = 5)

seu <- FindNeighbors(seu, dims = 1:20, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.2, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:20, verbose = FALSE)

save(seu, file = seurat_obj_file)
