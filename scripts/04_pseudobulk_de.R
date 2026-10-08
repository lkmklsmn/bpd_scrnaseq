# 04 - Per-cell-type pseudo-bulk differential expression.
#
# Ported from code/figures_for_paper_dblt_removed.R (lines 386-418).
#
# Input   frozen Seurat object
# Outputs Fig 3A    number of DE genes per cell type
#         Table S3  full DESeq2 results for every cell type
#         de_by_celltype.rds  reused by 05_gsea, 06_macrophages, 08_drug_targets
#
# DESIGN NOTE: four pooled libraries, two per condition. Aggregating per
# library is replicate-aware but underpowered (2 vs 2), which is why the
# manuscript uses an adjusted-P threshold of 0.25 and leans on external
# replication rather than on these p-values alone.

source(file.path("config", "paths.R"))
source(file.path("R", "figure_io.R"))
source(file.path("R", "plots.R"))
suppressMessages(library(writexl))

seu <- load_annotated()

des_deseq <- run_de(seu)
saveRDS(des_deseq, deseq_rds)

# ---- Table S3 ---------------------------------------------------------------
write_xlsx(lapply(des_deseq, data.frame), path = deseq_xlsx)
message("  TableS3 -> ", basename(deseq_xlsx), " (",
        length(des_deseq), " cell types)")

# ---- the numbers quoted in the Results -------------------------------------
cnt <- de_counts(des_deseq, padj = padj_de)
am <- cnt[cnt$cell_type == "Alveolar macrophages", ]
expect_value("alveolar macrophage DE genes", am$n_total, 361)
expect_value("alveolar macrophage genes up", am$n_up, 262)
expect_value("alveolar macrophage genes down", am$n_down, 99)

im <- cnt[cnt$cell_type == "Interstitial macrophages", ]
expect_value("interstitial macrophage DE genes", im$n_total, 98)
expect_value("interstitial macrophage genes up", im$n_up, 96)

# ---- Fig 3A -----------------------------------------------------------------
save_panel("Fig3A", plot_de(des_deseq, padj = padj_de),
  data = cnt[order(-cnt$n_total), ], width = 6, height = 6,
  slug = "de_counts"
)

message("04_pseudobulk_de.R complete")
