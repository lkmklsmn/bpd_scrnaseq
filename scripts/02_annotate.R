# 02 - Cell-type annotation and the immune atlas.
#
# Ported from code/figures_for_paper_dblt_removed.R (lines 239-336) in the
# analysis repository, which is the code that produced the submitted figures.
#
# Input   config/paths.R -> seurat_obj_file (frozen object, see data/README.md)
# Outputs Fig 1A  UMAP coloured by cell type
#         Fig 1B  UMAP coloured by exposure
#         Fig 1C  heatmap, top 100 markers per cell type
#         Fig 1D  dot plot, one canonical marker per cell type
#         Fig S2  immunoglobulin variable genes across the B-cell clusters
#         Table S1  every significant marker gene
#
# Marker table: FindAllMarkers (Wilcoxon, positive only, min.pct = 0.1,
# logfc.threshold = 0.25); genes with Bonferroni-adjusted P < 0.05 are kept and
# ranked within each cell type by detection difference (pct.1 - pct.2). The top
# 100 per cell type go into Fig 1C; all of them go into Table S1.

source(file.path("config", "paths.R"))
source(file.path("R", "figure_io.R"))
suppressMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(pheatmap)
  library(ggplot2)
})

seu <- load_annotated()

# ---- Fig 1A/B: UMAPs --------------------------------------------------------
umap_xy <- as.data.frame(Embeddings(seu, "umap"))
names(umap_xy) <- c("umap_1", "umap_2")
umap_xy$celltype <- seu@meta.data$celltype
umap_xy$condition <- seu@meta.data$condition

save_panel(
  "Fig1A",
  DimPlot(seu, group.by = "celltype", label = FALSE, repel = TRUE, raster = TRUE),
  data = umap_xy[, c("umap_1", "umap_2", "celltype")],
  width = 9, height = 7, slug = "umap_celltype"
)
save_panel(
  "Fig1B",
  DimPlot(seu, group.by = "condition", label = FALSE, repel = TRUE, raster = TRUE),
  data = umap_xy[, c("umap_1", "umap_2", "condition")],
  width = 7.5, height = 7, slug = "umap_condition"
)

# ---- marker discovery -------------------------------------------------------
Idents(seu) <- "celltype"
des <- FindAllMarkers(seu,
  only.pos = TRUE, min.pct = min_pct,
  logfc.threshold = 0.25, verbose = FALSE
)
des$diff <- des$pct.1 - des$pct.2

sig <- des |>
  filter(p_val_adj < padj_marker) |>
  group_by(cluster) |>
  arrange(desc(diff), .by_group = TRUE) |>
  ungroup()

expect_value("significant marker genes (Table S1)", nrow(sig), 8388)
save_table("TableS1_celltype_markers", sig, path = markers_csv)

top_n_markers <- sig |>
  group_by(cluster) |>
  slice_max(order_by = diff, n = n_heatmap, with_ties = FALSE) |>
  ungroup()

# ---- Fig 1C: top-100 marker heatmap, row-scaled mean expression ------------
sct <- GetAssayData(seu, assay = "SCT", layer = "data")
cells_by_type <- split(seq_len(ncol(seu)), seu@meta.data$celltype)
means <- do.call(cbind, lapply(cells_by_type, function(idx) {
  Matrix::rowMeans(sct[match(top_n_markers$gene, rownames(sct)), idx])
}))
means <- means[, match(unique(top_n_markers$cluster), colnames(means))]

p1c <- pheatmap::pheatmap(means,
  cluster_rows = FALSE, cluster_cols = FALSE,
  scale = "row", show_rownames = FALSE,
  color = viridisLite::inferno(100),
  silent = TRUE
)
save_panel("Fig1C", p1c,
  data = data.frame(
    gene = top_n_markers$gene, cell_type = top_n_markers$cluster,
    avg_log2FC = top_n_markers$avg_log2FC, pct_in_type = top_n_markers$pct.1,
    pct_other_types = top_n_markers$pct.2, pct_difference = top_n_markers$diff
  ),
  width = 7, height = 5, slug = "marker_heatmap"
)

# ---- Fig 1D: one canonical marker per cell type -----------------------------
# Defined in config/paths.R so the figure, Table S1 and the manuscript cannot
# drift apart. Every marker must sit inside the top 100 of its own cell type.
stopifnot(all(names(fig1d_markers) %in% unique(sig$cluster)))
in_top <- vapply(names(fig1d_markers), function(ct) {
  fig1d_markers[[ct]] %in% top_n_markers$gene[top_n_markers$cluster == ct]
}, logical(1))
if (!all(in_top)) {
  stop("Fig 1D markers outside the top ", n_heatmap, " of their cell type: ",
       paste(sprintf("%s (%s)", fig1d_markers[!in_top],
                     names(fig1d_markers)[!in_top]), collapse = ", "))
}

p1d <- DotPlot(seu, features = unname(fig1d_markers)) +
  labs(x = "Marker gene", y = "Cell type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_panel("Fig1D", p1d, data = p1d$data, width = 8, height = 5,
           slug = "dotplot_markers")

# ---- Fig S2: B-cell immunoglobulin variable genes ---------------------------
# Supports the statement that the six B-cell clusters separate on clonal Ig
# variable-gene diversity rather than distinct cell states.
bcells <- seu[, seu@meta.data$celltype == "B cells"]
xy <- Embeddings(bcells, "umap")
keep <- which(xy[, 2] > 2.5 & xy[, 2] < 17 & xy[, 1] > -10 & xy[, 1] < 5)
ig_genes <- c("Igkc", "Igkv1-117", "Igkv1-135", "Igkv9-120", "Igkv1-110", "Iglv1")
save_panel("FigS2", FeaturePlot(bcells[, keep], features = ig_genes),
           width = 8, height = 10, slug = "bcell_ig")

message("02_annotate.R complete")
