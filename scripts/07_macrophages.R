# 07 - Alveolar vs interstitial macrophage programs (Figure 4; Fig S5, S8).
# Source: figures_for_paper_dblt_removed.R ("Focus on Macrophages")
# Input : Seurat_object_scrublet.RData, des_deseq.rds

source(file.path("config", "paths.R"))
source(file.path("scripts", "utils", "functions.R"))
source(file.path("scripts", "utils", "load_annotate.R"))   # -> seu
library(Seurat); library(ggplot2); library(fgsea); library(msigdbr)

des_deseq <- readRDS(file.path(tab_dir, "des_deseq.rds"))
paths <- msigdbr(species = "Mus musculus", category = "H")
paths <- split(paths$gene_symbol, paths$gs_name)

reembed <- function(cell){
  sub <- seu[, seu$celltype == cell]
  DefaultAssay(sub) <- "RNA"
  sub <- NormalizeData(sub, verbose = FALSE)
  sub <- FindVariableFeatures(sub, nfeatures = 3000, verbose = FALSE)
  sub <- ScaleData(sub, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), verbose = FALSE)
  sub <- RunPCA(sub, npcs = 50, verbose = FALSE)
  sub <- FindNeighbors(sub, dims = 1:20, verbose = FALSE)
  RunUMAP(sub, dims = 1:20, verbose = FALSE)
}

# ---- Fig 4A/B: condition-colored UMAPs (air = red, O2 = blue in the figure) -
DimPlot(reembed("Alveolar macrophages"), group.by = "condition")
ggsave(file.path(fig_dir, "Fig4A_umap_alv_mac.pdf"), height = 7, width = 8)
DimPlot(reembed("Interstitial macrophages"), group.by = "condition")
ggsave(file.path(fig_dir, "Fig4B_umap_int_mac.pdf"), height = 7, width = 8)

# ---- Fig 4C/D: volcanoes ----------------------------------------------------
create_volcano2(des_deseq[["Alveolar macrophages"]])
ggsave(file.path(fig_dir, "Fig4C_volcano_alv_mac.pdf"), height = 5, width = 4)
create_volcano2(des_deseq[["Interstitial macrophages"]])
ggsave(file.path(fig_dir, "Fig4D_volcano_int_mac.pdf"), height = 5, width = 4)

# ---- Fig 4E: representative genes (Fabp1/Cd63 AM; Il1b/Cdkn1a IM) -----------
p <- gridExtra::grid.arrange(
  create_box_plot("Fabp1",  "Alveolar macrophages",    seu),
  create_box_plot("Cd63",   "Alveolar macrophages",    seu),
  create_box_plot("Il1b",   "Interstitial macrophages", seu),
  create_box_plot("Cdkn1a", "Interstitial macrophages", seu), nrow = 1)
ggsave(p, filename = file.path(fig_dir, "Fig4E_representative_genes.pdf"), height = 5, width = 10)

# ---- Fig 4F: interferon enrichment (union of Hallmark IFN-a/-g) -------------
genes <- union(paths[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]],
               paths[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]])
enr <- function(cell, title){
  x <- na.omit(des_deseq[[cell]]); stats <- x$log2FoldChange; names(stats) <- x$gene
  fgsea::plotEnrichment(pathway = genes, stats = stats) + labs(title = title) + theme_classic()
}
p <- gridExtra::grid.arrange(enr("Alveolar macrophages", "Alveolar macrophages"),
                             enr("Interstitial macrophages", "Interstitial macrophages"), ncol = 1)
ggsave(p, filename = file.path(fig_dir, "Fig4F_interferon_enrichment.pdf"), height = 6, width = 4)

# ---- Fig 4G / S8: representative interferon gene Ifi30 -----------------------
p <- gridExtra::grid.arrange(
  create_box_plot("Ifi30", "Alveolar macrophages",    seu),
  create_box_plot("Ifi30", "Interstitial macrophages", seu), nrow = 1)
ggsave(p, filename = file.path(fig_dir, "Fig4G_ifi30_boxplot.pdf"), height = 5, width = 5)

# ---- Fig S5: Fabp1 / Fabp4 cell-type specificity ----------------------------
VlnPlot(seu, features = c("Fabp1", "Fabp4"), group.by = "celltype", pt.size = 0, ncol = 1) &
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig_dir, "S5_fabp_specificity.pdf"), height = 8, width = 7)
