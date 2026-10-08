# 05 - Hallmark gene-set enrichment and the fold-change PCA.
#
# Ported from code/figures_for_paper_dblt_removed.R (lines 433-495).
#
# Input   de_by_celltype.rds (from 04)
# Outputs Fig 3B    Hallmark NES heatmap across cell types
#         Fig 3C    interferon-alpha / -gamma NES across cell types
#         Fig 3D    PCA of per-cell-type log2 fold-change profiles
#         Table S4  full fgsea results per cell type
#
# NOTE The working script also ran GO:BP enrichment from a local GMT file.
# Those results do not appear in the manuscript, so that pass is omitted here
# and the GMT file is not a dependency of this repository.
#
# REPRODUCIBILITY Hallmark set membership changes between MSigDB releases, so
# the msigdbr version is recorded below and the NES values quoted in the
# Results are checked. These checks WARN rather than stop: a mismatch means
# the installed msigdbr differs from the one used for the paper, so restore it
# with renv::restore() before interpreting Fig 3B/3C or Table S4.

source(file.path("config", "paths.R"))
source(file.path("R", "figure_io.R"))
source(file.path("R", "plots.R"))
suppressMessages({
  library(msigdbr)
  library(fgsea)
  library(writexl)
  library(ggplot2)
})

des_deseq <- readRDS(deseq_rds)

# ---- Hallmark gene sets (mouse) ---------------------------------------------
msig_version <- as.character(packageVersion("msigdbr"))
message("  msigdbr version: ", msig_version)
hallmark <- msigdbr(species = "Mus musculus", collection = "H")
paths <- split(hallmark$gene_symbol, hallmark$gs_name)
expect_value("Hallmark gene sets", length(paths), 50, on_fail = "warn")

enrich <- run_enrich(des_deseq, paths)

# ---- Table S4 ---------------------------------------------------------------
# leadingEdge is a list column; flatten it so it survives the round trip.
flat <- lapply(enrich, function(x) {
  x$leadingEdge <- vapply(x$leadingEdge, paste, character(1), collapse = ";")
  as.data.frame(x)
})
write_xlsx(flat, path = gsea_xlsx)
message("  TableS4 -> ", basename(gsea_xlsx), " (",
        sum(vapply(flat, nrow, integer(1))), " rows across ",
        length(flat), " cell types)")
expect_value("Table S4 rows", sum(vapply(flat, nrow, integer(1))), 412,
             on_fail = "warn")

# ---- Fig 3B: Hallmark NES heatmap ------------------------------------------
m <- nes_matrix(enrich, padj = padj_de)
save_panel("Fig3B", plot_enrich(enrich, padj = padj_de),
  data = data.frame(pathway = rownames(m), m, check.names = FALSE),
  width = 8, height = 10, slug = "hallmark_nes_heatmap"
)

# ---- Fig 3C: the interferon dichotomy --------------------------------------
ifn <- do.call(rbind, lapply(names(enrich), function(ct) {
  x <- enrich[[ct]]
  data.frame(celltype = ct, x[grep("INTERFERON", x$pathway), ])
}))
lev <- rev(names(sort(vapply(split(ifn$NES, ifn$celltype), mean, numeric(1)))))
ifn$celltype <- factor(ifn$celltype, levels = lev)

# the contrast the Results are built on
get_nes <- function(ct, pathway) {
  ifn$NES[ifn$celltype == ct & ifn$pathway == pathway]
}
expect_value("interstitial macrophage IFN-alpha NES",
             get_nes("Interstitial macrophages",
                     "HALLMARK_INTERFERON_ALPHA_RESPONSE"), 2.14,
             tol = 0.05, on_fail = "warn")
expect_value("alveolar macrophage IFN-alpha NES",
             get_nes("Alveolar macrophages",
                     "HALLMARK_INTERFERON_ALPHA_RESPONSE"), -1.90,
             tol = 0.05, on_fail = "warn")

p3c <- ggplot(ifn, aes(.data$celltype, .data$NES, fill = .data$NES)) +
  facet_wrap(~pathway, ncol = 1) +
  geom_bar(stat = "identity") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_panel("Fig3C", p3c,
  data = ifn[, c("celltype", "pathway", "NES", "pval", "padj", "size")],
  width = 7, height = 7, slug = "interferon_nes"
)

# ---- Fig 3D: PCA of per-cell-type fold-change profiles ---------------------
# Restricted to genes differentially expressed in at least one cell type;
# missing values set to 0 so every cell type spans the same gene space.
all_genes <- unique(unlist(lapply(des_deseq, function(x) {
  rownames(x)[!is.na(x$padj) & x$padj < padj_de]
})))
fc_mat <- do.call(cbind, lapply(des_deseq, function(x) {
  x$log2FoldChange[match(all_genes, rownames(x))]
}))
rownames(fc_mat) <- all_genes
fc_mat[is.na(fc_mat)] <- 0
fc_mat <- fc_mat[apply(fc_mat, 1, var) > 0, ]

pca <- prcomp(t(fc_mat), scale. = TRUE)
ve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
message(sprintf("  fold-change PCA: PC1 %.1f%%, PC2 %.1f%%", ve[1], ve[2]))

pc <- data.frame(cell = rownames(pca$x), pca$x)
p3d <- ggplot(pc, aes(.data$PC1, .data$PC2, label = .data$cell)) +
  geom_point() +
  ggrepel::geom_text_repel() +
  theme_classic()
save_panel("Fig3D", p3d,
  data = data.frame(
    cell_type = pc$cell, PC1 = pc$PC1, PC2 = pc$PC2,
    PC1_percent = round(ve[1], 2), PC2_percent = round(ve[2], 2)
  ),
  width = 5, height = 5, slug = "foldchange_pca"
)

message("05_gsea.R complete")
