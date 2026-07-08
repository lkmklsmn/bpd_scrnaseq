# aux - Hurskainen (GSE151974) cross-study validation.
# Source: frequency_analysis_gse151974.R + correlation_binoy_gse151974.R
# Produces:
#   Fig 2C  cross-study scatter of cell-type proportion changes (this study vs Hurskainen)
#   Fig S1  fold-change correlation heatmap (this study vs Hurskainen cell types)
#   Fig S9  Lgals3 in Hurskainen alveolar macrophages, hyperoxia vs normoxia per timepoint
# Inputs: data/raw/Seurat_object_scrublet.RData, data/external/GSE151974/Seurat_object_GSE151974.RData

source(file.path("config", "paths.R"))
suppressMessages({ library(Seurat); library(dplyr); library(ggplot2); library(pheatmap); library(patchwork) })

load(file.path(external_dir, "GSE151974", "Seurat_object_GSE151974.RData"))   # -> gse151974

# =============================================================================
# (A) Hurskainen cell-type proportion change (hyperoxia vs normoxia, per timepoint)
# =============================================================================
meta <- gse151974@meta.data
prop <- t(unclass(table(meta$Sample, meta$CellType)) /
          rowSums(unclass(table(meta$Sample, meta$CellType))))
offset <- 1e-6
logit  <- log((prop + offset) / (1 - prop + offset))
treat  <- do.call(rbind, lapply(colnames(prop), function(x) strsplit(x, "_", fixed = TRUE)[[1]]))
day <- treat[, 1]; condition <- treat[, 2]

# per-cell-type hyperoxia effect (logit-proportion LM, adjusting for day)
hur_effect <- t(apply(logit, 1, function(x){
  cf <- coefficients(summary(lm(x ~ day + condition)))
  cf[grep("^condition", rownames(cf)), c(1, 4)]
}))
colnames(hur_effect) <- c("coef", "pval")
hur_effect <- data.frame(hur_effect, cell = rownames(hur_effect))
hur_effect$coef <- -hur_effect$coef   # orient so positive = up in hyperoxia

# P14 hyperoxia/air percent change (timepoint matching this study)
percent_change <- prop[, c(1, 3, 5)] / prop[, c(2, 4, 6)] * 100
percent_change <- percent_change[, c(2, 3, 1)]   # -> P3, P7, P14 order

# =============================================================================
# (B) This study: cell-type proportion change (propeller output from script 04)
# =============================================================================
source(file.path("scripts", "utils", "load_annotate.R"))   # -> seu (+ celltype)
this_prop <- prop.table(table(seu$condition, seu$celltype), margin = 1)
this_change <- log2(this_prop["O2", ] / this_prop["Air", ])

# =============================================================================
# (C) Fig 2C - cross-study scatter (this study vs Hurskainen, P14)
#     NOTE: cell-type vocabularies differ between datasets; VERIFY this mapping.
# =============================================================================
map_ct <- c(
  "NK cells"                 = "NK cell",
  "Interstitial macrophages" = "Int Mf",
  "Alveolar macrophages"     = "Alv Mf",
  "T cells"                  = "gd T cell",   # adjust to the intended T-cell match
  "B cells"                  = "B cell",
  "Dendritic cells"          = "DC1",
  "Monocytes"                = "Mono",
  "Neutrophils"              = "Neutrophil")
hur_p14 <- log2(percent_change[, "P14_Hyperoxia"] / 100)   # log2 fold vs air at P14
comp <- data.frame(
  celltype   = names(map_ct),
  this_study = this_change[names(map_ct)],
  hurskainen = hur_p14[map_ct])
comp <- na.omit(comp)
rho <- cor(comp$this_study, comp$hurskainen, method = "spearman")
ggplot(comp, aes(this_study, hurskainen, label = celltype)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point() + ggrepel::geom_text_repel() +
  labs(x = "log2 proportion change (this study)",
       y = "log2 proportion change (Hurskainen, P14)",
       title = sprintf("Cross-study concordance (Spearman rho = %.2f)", rho)) +
  theme_classic()
ggsave(file.path(fig_dir, "Fig2C_crossstudy_scatter.pdf"), height = 5, width = 6)

# =============================================================================
# (D) Fig S1 - marker fold-change correlation heatmap (this study vs Hurskainen)
# =============================================================================
binoy_des      <- FindAllMarkers(seu,       group.by = "celltype", max.cells.per.ident = 100)
hurskainen_des <- FindAllMarkers(gse151974, group.by = "CellType", max.cells.per.ident = 100)
bt <- unique(binoy_des$cluster); ht <- unique(hurskainen_des$cluster)
cor_matrix <- matrix(NA, length(bt), length(ht), dimnames = list(bt, ht))
for (b in bt) for (h in ht) {
  shared <- inner_join(binoy_des[binoy_des$cluster == b, c("gene", "avg_log2FC")],
                       hurskainen_des[hurskainen_des$cluster == h, c("gene", "avg_log2FC")],
                       by = "gene", suffix = c("_b", "_h"))
  if (nrow(shared) >= 10) cor_matrix[b, h] <- cor(shared$avg_log2FC_b, shared$avg_log2FC_h)
}
p <- pheatmap(cor_matrix,
  color = colorRampPalette(c("#2166ac", "white", "#d6604d"))(100),
  breaks = seq(-1, 1, length.out = 101),
  display_numbers = TRUE, number_format = "%.2f", fontsize_number = 7,
  main = "Fold-change correlation (rows = this study, cols = Hurskainen)")
ggsave(p, filename = file.path(fig_dir, "S1_hurskainen_correlation.pdf"), height = 7.5, width = 10)

# =============================================================================
# (E) Fig S9 - Lgals3 in Hurskainen alveolar macrophages by timepoint
# =============================================================================
am <- gse151974[, gse151974$CellType == "Alv Mf"]
DefaultAssay(am) <- "RNA"; am <- NormalizeData(am, verbose = FALSE)
am$Oxygen <- factor(am$Oxygen, levels = c("Normoxia", "Hyperoxia"))
am$Age    <- factor(am$Age,    levels = c("P3", "P7", "P14"))
vln <- lapply(levels(am$Age), function(tp)
  VlnPlot(am[, am$Age == tp], features = "Lgals3", group.by = "Oxygen", pt.size = 1) +
    scale_fill_manual(values = c("Normoxia" = "grey", "Hyperoxia" = "red")) +
    labs(title = tp, x = NULL, y = "Lgals3 expression") +
    ggpubr::stat_compare_means() + theme_classic() + theme(legend.position = "none"))
wrap_plots(vln, nrow = 1)
ggsave(file.path(fig_dir, "S9_hurskainen_lgals3_am_vln.pdf"), height = 4, width = 9)
