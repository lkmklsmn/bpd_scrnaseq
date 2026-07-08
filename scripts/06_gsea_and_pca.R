# 06 - Hallmark GSEA, interferon dichotomy, fold-change PCA (Fig 3B/C/D; Table S2).
# Source: figures_for_paper_dblt_removed.R (enrichment + interferon + PCA)
# Input : des_deseq.rds (from 05)

source(file.path("config", "paths.R"))
source(file.path("scripts", "utils", "functions.R"))
library(msigdbr); library(fgsea); library(ggplot2)

des_deseq <- readRDS(file.path(tab_dir, "des_deseq.rds"))

# ---- Hallmark gene sets (pin msigdbr version for reproducibility) -----------
paths <- msigdbr(species = "Mus musculus", category = "H")
paths <- split(paths$gene_symbol, paths$gs_name)

enrich_h <- run_enrich(des_deseq, paths)                    # Table S2
tmp <- lapply(enrich_h, function(x){ x$leadingEdge <- as.character(x$leadingEdge); x })
writexl::write_xlsx(tmp, path = file.path(tab_dir, "gsea_hallmark_by_celltype.xlsx"))

# ---- Fig 3B: Hallmark NES heatmap -------------------------------------------
p <- plot_enrich(enrich_h)
ggsave(p, filename = file.path(fig_dir, "Fig3B_hallmark_nes_heatmap.pdf"), height = 10, width = 8)

# ---- Fig 3C: interferon NES across cell types -------------------------------
ifn <- do.call(rbind, lapply(names(enrich_h), function(x){
  subm <- enrich_h[[x]]; data.frame(celltype = x, subm[grep("INTERFERON", subm$pathway), ])
}))
lev <- rev(names(sort(vapply(split(ifn$NES, ifn$celltype), mean, numeric(1)))))
ifn$celltype <- factor(ifn$celltype, levels = lev)
ggplot(ifn, aes(celltype, NES, fill = NES)) +
  facet_wrap(~ pathway, ncol = 1) + geom_bar(stat = "identity") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig_dir, "Fig3C_interferon_nes.pdf"), height = 7, width = 7)

# ---- Fig 3D: PCA of per-cell-type fold-change profiles ----------------------
all_genes <- unique(unlist(lapply(des_deseq, function(x) rownames(x)[x$padj < 0.25])))
fc_mat <- do.call(cbind, lapply(des_deseq, function(x) x$log2FoldChange[match(all_genes, rownames(x))]))
rownames(fc_mat) <- all_genes; fc_mat[is.na(fc_mat)] <- 0
fc_mat <- fc_mat[apply(fc_mat, 1, var) > 0, ]
pca <- prcomp(t(fc_mat), scale. = TRUE)
ggplot(data.frame(cell = rownames(pca$x), pca$x), aes(PC1, PC2, label = cell)) +
  geom_point() + ggrepel::geom_text_repel() + theme_classic()
ggsave(file.path(fig_dir, "Fig3D_foldchange_pca.pdf"), height = 5, width = 5)
