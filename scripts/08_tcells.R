# 08 - T-cell subclustering and gamma-delta expansion (Figure 5).
# Source: figures_for_paper_dblt_removed.R ("Focus on T cells")
# Input : Seurat_object_scrublet.RData

source(file.path("config", "paths.R"))
source(file.path("scripts", "utils", "functions.R"))
source(file.path("scripts", "utils", "load_annotate.R"))   # -> seu
library(Seurat); library(ggplot2); library(DESeq2)

sub <- seu[, seu$celltype == "T cells"]
DefaultAssay(sub) <- "RNA"
sub <- NormalizeData(sub, verbose = FALSE)
sub <- FindVariableFeatures(sub, nfeatures = 3000, verbose = FALSE)
sub <- ScaleData(sub, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), verbose = FALSE)
sub <- RunPCA(sub, npcs = 20, verbose = FALSE)
sub <- FindNeighbors(sub, dims = 1:20, verbose = FALSE)
sub <- RunUMAP(sub, dims = 1:20, verbose = FALSE)
sub <- FindClusters(sub, resolution = 0.1, verbose = FALSE)

tcell_annotations <- c("0" = "CD4 T cells", "1" = "CD8 T cells",
                       "2" = "dg T cells", "3" = "Th2 cells", "4" = "Tregs")
sub@meta.data$celltype <- unname(tcell_annotations[as.character(sub$seurat_clusters)])

# ---- Fig 5A: subtype UMAP ---------------------------------------------------
DimPlot(sub, group.by = "celltype", label = TRUE)
ggsave(file.path(fig_dir, "Fig5A_umap_tcells.pdf"), height = 7, width = 8)

# ---- Fig 5B: subtype marker dot plot ----------------------------------------
tcell_markers <- c("Lef1", "Cd8b1", "Trdc", "Il1rl1", "Foxp3")
DotPlot(sub, features = tcell_markers, group.by = "celltype") +
  labs(x = "Marker gene", y = "Cell type") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig_dir, "Fig5B_dotplot_tcell_markers.pdf"), height = 5, width = 6)

# ---- Fig 5C: differential abundance across subtypes -------------------------
res <- run_differential_frequency(sub)
plot_differential_frequency(res)
ggsave(file.path(fig_dir, "Fig5C_diff_freq_tcell.pdf"), height = 5, width = 6)

# ---- Fig 5D: gamma-delta proportion by condition ----------------------------
cell <- "dg T cells"
fr <- table(sub@meta.data$orig.ident[sub@meta.data$celltype == cell]) /
      table(sub@meta.data$orig.ident)
subm <- data.frame(fr, sub@meta.data[match(names(fr), sub@meta.data$orig.ident), ])
ggplot(subm, aes(condition, Freq * 100, color = condition)) +
  labs(y = "gd T-cell frequency [% T cells]", x = "Condition") +
  geom_boxplot() + geom_point() +
  scale_color_manual(values = c("black", "red")) + theme_classic()
ggsave(file.path(fig_dir, "Fig5D_gd_frequency.pdf"), height = 5, width = 3)

# ---- Fig 5E/F: IL-17 module -------------------------------------------------
il17_genes <- c("Rorc", "Il23r", "Il1r1", "Il17a", "Il17f")
DefaultAssay(sub) <- "RNA"
sub <- AddModuleScore(sub, features = list(il17_genes), name = "IL17_module", assay = "RNA")
FeaturePlot(sub, features = "IL17_module1", min.cutoff = "q05", max.cutoff = "q95") +
  labs(title = "IL-17 module score") + theme_classic()
ggsave(file.path(fig_dir, "Fig5E_il17_module_umap.pdf"), height = 7, width = 8)

# per-cell-type pseudobulk module score, gd T cells, by condition
module_box <- function(genes, subtype){
  ss <- sub[, sub@meta.data$celltype == subtype]
  bsplit <- split(1:ncol(ss), ss@meta.data$sample)
  sums <- do.call(cbind, lapply(bsplit, function(k) rowSums(ss@assays$SCT@counts[, k])))
  sums <- sums[rowSums(sums) > 0, ]
  treat <- ss@meta.data$condition[match(colnames(sums), ss@meta.data$sample)]
  dds <- DESeqDataSetFromMatrix(sums, data.frame(treat), ~ treat); dds <- DESeq(dds, quiet = TRUE)
  norm <- counts(dds, normalized = TRUE); genes <- intersect(genes, rownames(norm))
  module <- colMeans(log1p(norm[genes, , drop = FALSE]))
  subm <- data.frame(module = module, ss@meta.data[match(names(module), ss@meta.data$sample), ])
  ggplot(subm, aes(condition, module, color = condition)) +
    labs(title = subtype, y = "IL-17 module score (pseudobulk)") +
    geom_boxplot() + geom_point() + ggpubr::stat_compare_means() +
    scale_color_manual(values = c("grey", "red")) + theme_classic()
}
module_box(il17_genes, "dg T cells")
ggsave(file.path(fig_dir, "Fig5F_il17_module_gd_boxplot.pdf"), height = 5, width = 3)
