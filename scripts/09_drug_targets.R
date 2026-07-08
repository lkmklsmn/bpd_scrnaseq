# 09 - Drug-target analysis of the alveolar-macrophage response (Figure 6).
# Source: drug_target_figures.R
# Input : Seurat_object_scrublet.RData, condition_deg_by_celltype_deseq.xlsx (from 05),
#         drug_target_candidates_curated.csv
# NOTE  : the curated CSV is produced by the babelgene -> Open Targets (GraphQL v4)
#         query step described in the Methods. Add that query script as
#         09a_opentargets_query.R (or .py) and cache its response in
#         data/external/opentargets/ for exact reproducibility (API version-dependent).

source(file.path("config", "paths.R"))
source(file.path("scripts", "utils", "load_annotate.R"))   # -> seu
library(ggplot2); library(ggrepel); library(readxl); library(dplyr); library(DESeq2); library(Seurat)

drug_file <- file.path(tab_dir, "drug_target_candidates_curated.csv")
drug <- read.csv(drug_file, stringsAsFactors = FALSE)
druggable_mouse <- unique(unlist(strsplit(drug$mouse_genes, "/")))
hasdrug_mouse   <- unique(unlist(strsplit(drug$mouse_genes[drug$tier == "A: existing drug"], "/")))

am <- as.data.frame(read_excel(deseq_xlsx, sheet = "Alveolar macrophages"))
am <- am[!is.na(am$padj), ]; am$pvalue[am$pvalue < 1e-50] <- 1e-50
am$drugclass <- ifelse(am$gene %in% hasdrug_mouse, "Druggable (existing drug)",
                ifelse(am$gene %in% druggable_mouse, "Druggable", "Not druggable"))
am$drugclass <- factor(am$drugclass, levels = c("Not druggable","Druggable","Druggable (existing drug)"))

# ---- Fig 6A: druggability-annotated volcano ---------------------------------
lab <- am[am$drugclass != "Not druggable" & am$log2FoldChange > 0.5 & am$pvalue < 1e-5, ]
lab <- lab[order(-lab$log2FoldChange), ]
ggplot(am[am$log2FoldChange > 0, ], aes(log2FoldChange, -log10(pvalue), color = drugclass)) +
  geom_point(alpha = 0.6, size = 1.4) +
  ggrepel::geom_text_repel(data = lab, aes(label = gene), color = "black", size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("Not druggable" = "grey80", "Druggable" = "#4393c3",
                                "Druggable (existing drug)" = "#d6604d"), name = NULL) +
  labs(x = "log2 fold change [O2/air]", y = "-log10 p-value") +
  theme_classic() + theme(legend.position = "top")
ggsave(file.path(fig_dir, "Fig6A_druggable_volcano.pdf"), width = 6, height = 6, device = cairo_pdf)

# ---- Fig 6B: ranked druggable targets ---------------------------------------
cand <- drug[grepl("Alveolar macrophages", drug$celltypes, ignore.case = TRUE) &
             drug$tier %in% c("A: existing drug", "B: clinically tractable"), ]
cand$am_fc   <- am$log2FoldChange[match(cand$mouse_genes, am$gene)]
cand$am_pval <- am$pvalue[match(cand$mouse_genes, am$gene)]
cand <- cand[order(-cand$am_fc), ]
cand$label_drug   <- sub(";.*$", "", cand$example_drugs)
cand$human_symbol <- factor(cand$human_symbol, levels = rev(cand$human_symbol))
ggplot(cand, aes(am_fc, human_symbol, color = tier)) +
  geom_segment(aes(x = 0, xend = am_fc, yend = human_symbol), color = "grey70", linewidth = 0.5) +
  geom_point(aes(size = -log10(am_pval))) +
  geom_text(aes(label = label_drug), hjust = -0.1, size = 3, color = "grey25") +
  scale_color_manual(values = c("A: existing drug" = "#d6604d",
                                "B: clinically tractable" = "#4393c3"), name = NULL) +
  xlim(0, max(cand$am_fc) * 1.45) + labs(x = "log2 fold change (O2/air)", y = NULL) +
  theme_classic() + theme(legend.position = "top")
ggsave(file.path(fig_dir, "Fig6B_ranked_targets.pdf"), width = 6.5, height = 6)

# ---- Fig 6C: Lgals3 specificity across CD45+ cell types ---------------------
VlnPlot(seu, features = "Lgals3", group.by = "celltype", pt.size = 0) +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Lgals3 across cell types", x = NULL)
ggsave(file.path(fig_dir, "Fig6C_lgals3_violin.pdf"), width = 7, height = 4)

# ---- Fig 6D: baseline-expression ECDF (Lgals3 highly expressed) -------------
am$rank <- rank(am$baseMean)
ggplot(am, aes(log(baseMean))) +
  labs(x = "Baseline expression [log]", y = "Percentile") + stat_ecdf() +
  ggrepel::geom_text_repel(data = am[am$gene %in% c("Lgals3","Spp1","Gpnmb","Txnrd1"), ],
                           aes(x = log(baseMean), y = rank / nrow(am), label = gene)) +
  theme_classic()
ggsave(file.path(fig_dir, "Fig6D_lgals3_ecdf.pdf"), width = 4, height = 4)
