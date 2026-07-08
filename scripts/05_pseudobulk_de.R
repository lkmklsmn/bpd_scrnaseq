# 05 - Per-cell-type pseudobulk differential expression (Figure 3A; Table S1).
# Source: figures_for_paper_dblt_removed.R (des_deseq + plot_de)
# Input : Seurat_object_scrublet.RData
# Output: condition_deg_by_celltype_deseq.xlsx  (-> used by 06, 07, 09)
#
# DESIGN NOTE: 4 pooled libraries (2 per condition). Pseudobulk aggregation per
# library is replicate-aware but underpowered (2 vs 2); read p-values with care.

source(file.path("config", "paths.R"))
source(file.path("scripts", "utils", "functions.R"))
source(file.path("scripts", "utils", "load_annotate.R"))   # -> seu
library(DESeq2); library(ggplot2)

des_deseq <- run_de(seu)   # per-cell-type DESeq2 (O2 vs Air), padj threshold 0.25

writexl::write_xlsx(lapply(des_deseq, data.frame), path = deseq_xlsx)  # Table S1

# ---- Fig 3A: DE-gene counts per cell type -----------------------------------
plot_de(des_deseq)
ggsave(file.path(fig_dir, "Fig3A_de_counts.pdf"), height = 6, width = 6)

saveRDS(des_deseq, file.path(tab_dir, "des_deseq.rds"))  # reused by 06/07
