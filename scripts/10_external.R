# 10 - Cross-study comparisons against Hurskainen et al. (GSE151974).
#
# Inputs  frozen Seurat object; data/external/GSE151974
# Outputs Fig 2C  cell-type frequency change, this study vs Hurskainen at P14
#         Fig 5E  gamma-delta T cells as a percentage of T cells, by timepoint
#         Fig S1  marker fold-change correlation across the two atlases
#
# PROVENANCE -----------------------------------------------------------------
# Fig 5E is ported from code/frequency_analysis_gse151974.R and Fig S1 from
# code/correlation_binoy_gse151974.R in the analysis repository.
#
# Fig 2C had NO source script - it was produced interactively and the code was
# not saved. It is reconstructed here from the description in the Results:
# per-cell-type proportion changes (hyperoxia vs air), both datasets restricted
# to CD45+ populations, Hurskainen taken at postnatal day 14 to match our own
# timepoint, compared by Spearman correlation. The reported value is rho =
# 0.64; the check below therefore warns rather than stops, because the
# cell-type grouping below is a reconstruction of an undocumented choice.
#
# REPRODUCIBILITY NOTE -------------------------------------------------------
# Fig S1 calls FindAllMarkers with max.cells.per.ident = 100, which subsamples
# at random. The original script set no seed, so the submitted panel cannot be
# reproduced exactly. The seed from config/paths.R is set here so that this
# repository is at least self-consistent run to run.

source(file.path("config", "paths.R"))
source(file.path("R", "figure_io.R"))
suppressMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
})

hur_file <- file.path(external_dir, "GSE151974",
                      "Seurat_object_GSE151974.RData")
if (!file.exists(hur_file)) {
  stop("Missing external data for GSE151974:\n  ", hur_file,
       "\nSee data/README.md.", call. = FALSE)
}

seu <- load_annotated()
env <- new.env(parent = emptyenv())
load(hur_file, envir = env)
gse151974 <- env[[ls(env)[1]]]
expect_value("Hurskainen cells loaded", ncol(gse151974), 61839)

# Hurskainen CD45+ populations grouped onto the nine used in this study.
# Non-immune clusters (fibroblasts, endothelium, epithelium, pericytes,
# mesothelium, lymphatics) are excluded so that both datasets are restricted
# to the CD45+ compartment. ILC2 and Mast Ba2 have no counterpart here, as do
# our cycling cells in Hurskainen.
hur_groups <- list(
  "B cells"                  = c("B cell 1", "B cell 2"),
  "T cells"                  = c("CD4 T cell 1", "CD4 T cell 2",
                                 "CD8 T cell 1", "CD8 T cell 2", "gd T cell"),
  "Monocytes"                = "Mono",
  "Alveolar macrophages"     = "Alv Mf",
  "Neutrophils"              = c("Neut 1", "Neut 2"),
  "Dendritic cells"          = c("DC1", "DC2"),
  "NK cells"                 = "NK cell",
  "Interstitial macrophages" = "Int Mf"
)

hm <- gse151974@meta.data
unmapped <- setdiff(unlist(hur_groups), unique(as.character(hm$CellType)))
if (length(unmapped)) {
  stop("Hurskainen cell types named in hur_groups are absent from the object: ",
       paste(unmapped, collapse = ", "), call. = FALSE)
}

# =============================================================================
# Fig 2C - frequency change in both studies, restricted to CD45+
# =============================================================================
hm$group <- NA_character_
for (g in names(hur_groups)) {
  hm$group[as.character(hm$CellType) %in% hur_groups[[g]]] <- g
}
hur_p14 <- hm[hm$Age == "P14" & !is.na(hm$group), ]

prop_by <- function(labels, grouping) {
  tab <- table(grouping, labels)
  sweep(tab, 1, rowSums(tab), "/")
}
hp <- prop_by(hur_p14$group, hur_p14$Oxygen)
hur_change <- log2(hp["Hyperoxia", ] / hp["Normoxia", ])

op <- prop_by(seu@meta.data$celltype, seu@meta.data$condition)
our_change <- log2(op["O2", ] / op["Air", ])

shared <- intersect(names(our_change), names(hur_change))
cmp <- data.frame(
  cell_type = shared,
  this_study = as.numeric(our_change[shared]),
  hurskainen_p14 = as.numeric(hur_change[shared])
)
rho <- cor(cmp$this_study, cmp$hurskainen_p14, method = "spearman")
message(sprintf("  Fig 2C: Spearman rho = %.2f over %d matched cell types",
                rho, nrow(cmp)))
expect_value("Fig 2C Spearman rho", rho, 0.64, tol = 0.05, on_fail = "warn")

p2c <- ggplot(cmp, aes(.data$this_study, .data$hurskainen_p14)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(size = 2.5, color = "#d6604d") +
  ggrepel::geom_text_repel(aes(label = .data$cell_type), size = 3,
                           color = "grey20") +
  labs(
    x = "log2 proportion change, this study (O2/air)",
    y = "log2 proportion change, Hurskainen P14 (hyperoxia/normoxia)",
    subtitle = sprintf("Spearman rho = %.2f", rho)
  ) +
  theme_classic()
save_panel("Fig2C", p2c, data = cmp, width = 6, height = 5,
           slug = "crossstudy_scatter")

# =============================================================================
# Fig 5E - gamma-delta T cells as a percentage of T cells
# =============================================================================
t_types <- hur_groups[["T cells"]]
tmat <- unclass(table(hm$Sample, hm$CellType))[, t_types, drop = FALSE]
gd <- data.frame(
  sample = rownames(tmat),
  percent = 100 * tmat[, "gd T cell"] / rowSums(tmat),
  condition = factor(sub(".*_", "", rownames(tmat)),
                     levels = c("Normoxia", "Hyperoxia")),
  day = factor(sub("_.*", "", rownames(tmat)), levels = c("P3", "P7", "P14"))
)

means <- tapply(gd$percent, gd$condition, mean)
expect_value("Fig 5E gd mean in normoxia (%)", means[["Normoxia"]], 5.4,
             tol = 0.15, on_fail = "warn")
expect_value("Fig 5E gd mean in hyperoxia (%)", means[["Hyperoxia"]], 10.2,
             tol = 0.15, on_fail = "warn")
higher <- all(vapply(split(gd, gd$day), function(d) {
  d$percent[d$condition == "Hyperoxia"] > d$percent[d$condition == "Normoxia"]
}, logical(1)))
expect_value("gd higher in hyperoxia at every timepoint", higher, TRUE)

p5e <- ggplot(gd, aes(.data$day, .data$percent, fill = .data$condition)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  scale_fill_manual(values = c("Normoxia" = "grey", "Hyperoxia" = "red")) +
  labs(
    title = "Hurskainen et al", x = "Postnatal day",
    y = "gd T cells frequency [% T cells]"
  ) +
  theme_classic()
save_panel("Fig5E", p5e, data = gd, width = 4, height = 4,
           slug = "hurskainen_gd_frequency")

# =============================================================================
# Fig S1 - marker fold-change correlation between the two atlases
# =============================================================================
set.seed(seed)
ours <- FindAllMarkers(seu, group.by = "celltype",
                       max.cells.per.ident = 100, verbose = FALSE)
set.seed(seed)
theirs <- FindAllMarkers(gse151974, group.by = "CellType",
                         max.cells.per.ident = 100, verbose = FALSE)

our_types <- unique(as.character(ours$cluster))
their_types <- unique(as.character(theirs$cluster))
cor_mat <- matrix(NA_real_, length(our_types), length(their_types),
                  dimnames = list(our_types, their_types))
n_mat <- cor_mat
for (b in our_types) {
  b_sub <- ours[ours$cluster == b, c("gene", "avg_log2FC")]
  for (h in their_types) {
    h_sub <- theirs[theirs$cluster == h, c("gene", "avg_log2FC")]
    shared <- dplyr::inner_join(b_sub, h_sub, by = "gene",
                                suffix = c("_b", "_h"))
    if (nrow(shared) >= 10) {
      cor_mat[b, h] <- cor(shared$avg_log2FC_b, shared$avg_log2FC_h)
      n_mat[b, h] <- nrow(shared)
    }
  }
}

ps1 <- pheatmap(cor_mat,
  color = colorRampPalette(c("#2166ac", "white", "#d6604d"))(100),
  breaks = seq(-1, 1, length.out = 101),
  cluster_rows = TRUE, cluster_cols = TRUE,
  display_numbers = TRUE, number_format = "%.2f",
  fontsize_number = 7, fontsize_row = 9, fontsize_col = 9,
  main = paste("Fold change correlation - This study (rows)",
               "and Hurskainen et al (columns)"),
  silent = TRUE
)
long <- data.frame(
  this_study = rep(rownames(cor_mat), times = ncol(cor_mat)),
  hurskainen = rep(colnames(cor_mat), each = nrow(cor_mat)),
  pearson_r = as.vector(cor_mat),
  n_shared_genes = as.vector(n_mat)
)
save_panel("FigS1", ps1, data = long[!is.na(long$pearson_r), ],
           width = 10, height = 7.5, slug = "hurskainen_correlation")

message("10_external.R complete")
