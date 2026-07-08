# 03 - Cell-type annotation and the immune atlas (Figure 1; Figs S2, S4).
# Source: figures_for_paper_dblt_removed.R (annotation, UMAPs, markers)
# Input : Seurat_object_scrublet.RData
# Output: Fig 1 (UMAPs, marker heatmap, dot plot), marker tables, Fig S2 (B-cell Ig)

source(file.path("config", "paths.R"))
source(file.path("scripts", "utils", "functions.R"))
source(file.path("scripts", "utils", "load_annotate.R"))   # -> seu (+ celltype)
library(dplyr); library(Seurat); library(ggplot2); library(pheatmap)

# ---- Fig 1A/B: UMAPs --------------------------------------------------------
DimPlot(seu, group.by = "celltype", raster = TRUE)
ggsave(file.path(fig_dir, "Fig1A_umap_celltype.pdf"), height = 7, width = 9)
DimPlot(seu, group.by = "condition", raster = TRUE)
ggsave(file.path(fig_dir, "Fig1B_umap_condition.pdf"), height = 7, width = 7.5)

# ---- marker discovery + tables ----------------------------------------------
Idents(seu) <- "celltype"
des <- FindAllMarkers(seu, only.pos = TRUE)
top100 <- des %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 100)
write.csv(top100, file.path(tab_dir, "top100_markers.csv"), row.names = FALSE)
top10 <- des %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 10)
write.csv(top10, file.path(tab_dir, "top10_markers.csv"), row.names = FALSE)

# ---- Fig 1C: top-10 marker heatmap (row-scaled mean expression) -------------
asplit <- split(1:ncol(seu), seu@meta.data$celltype)
means <- do.call(cbind, lapply(asplit, function(x)
  rowMeans(seu@assays$SCT@data[match(top10$gene, rownames(seu@assays$SCT@data)), x])))
means <- means[, match(unique(top10$cluster), colnames(means))]
p <- pheatmap(means, cluster_rows = FALSE, cluster_cols = FALSE, scale = "row",
              show_rownames = FALSE, color = viridisLite::inferno(100))
ggsave(p, filename = file.path(fig_dir, "Fig1C_heatmap_markers.pdf"), height = 5, width = 7)

# ---- Fig 1D: canonical-marker dot plot --------------------------------------
group_markers <- c("Cd79a","Cd3e","C1qb","S100a9","Car4","Ace","Nkg7","Wdfy4","Mki67")
DotPlot(seu, features = group_markers) +
  labs(x = "Marker gene", y = "Cell type") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig_dir, "Fig1D_dotplot_markers.pdf"), height = 5, width = 8)

# ---- QC metrics per cell type (Fig S3 companion) ----------------------------
VlnPlot(seu, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
        group.by = "celltype", pt.size = 0, ncol = 1)
ggsave(file.path(fig_dir, "S3_vln_metrics.pdf"), height = 8, width = 5)

# ---- Fig S2: B-cell immunoglobulin variable genes ---------------------------
tmp <- seu[, seu$celltype == "B cells"]
ok <- which(tmp@reductions$umap@cell.embeddings[, 2] > 2.5 &
            tmp@reductions$umap@cell.embeddings[, 2] < 17 &
            tmp@reductions$umap@cell.embeddings[, 1] > -10 &
            tmp@reductions$umap@cell.embeddings[, 1] < 5)
FeaturePlot(tmp[, ok],
  features = c("Igkc","Igkv1-117","Igkv1-135","Igkv9-120","Igkv1-110","Iglv1"))
ggsave(file.path(fig_dir, "S2_bcell_ig_featureplot.pdf"), height = 10, width = 8)
