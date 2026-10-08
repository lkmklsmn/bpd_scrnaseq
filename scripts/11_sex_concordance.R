# 11 - Sex as a biological variable: concordance of the hyperoxia response.
#
# Ported from code/figure_S3_sex_concordance.R in the analysis repository.
#
# Input   frozen Seurat object; Table S3 (from 04)
# Outputs Fig S3A  PCA of all-cell pseudo-bulk profiles
#         Fig S3B  library-library correlation matrix
#         Fig S3C  alveolar macrophage fold-change concordance between sexes
#
# DESIGN ---------------------------------------------------------------------
# Four libraries, one male and one female per exposure:
#   library 1  hyperoxia, female      library 3  room air, female
#   library 2  room air,  male        library 4  hyperoxia, male
# Sex is therefore balanced across the exposure contrast and cannot confound
# it, but is confounded with library within each exposure, so a sex x exposure
# interaction cannot be tested. These panels ask only whether the main
# findings depend on sex.
#
# Nothing in the primary analysis is recomputed: both inputs are the same
# objects used for the main figures.

source(file.path("config", "paths.R"))
source(file.path("R", "figure_io.R"))
suppressMessages({
  library(Seurat)
  library(Matrix)
  library(readxl)
  library(ggplot2)
  library(reshape2)
})

n_var_heatmap <- 2000   # variable genes for the correlation matrix
n_var_pca     <- 5000   # variable genes for the PCA

seu <- load_annotated()

counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
lib <- as.character(seu$orig.ident)
lib_meta <- data.frame(
  lib      = c("1", "2", "3", "4"),
  exposure = c("Hyperoxia", "Room air", "Room air", "Hyperoxia"),
  sex      = c("Female", "Male", "Female", "Male"),
  stringsAsFactors = FALSE
)
lib_meta$label <- paste(lib_meta$exposure, tolower(lib_meta$sex), sep = ", ")

# Confirm the design encoded above matches the object, rather than assuming it
observed <- unique(seu@meta.data[, c("orig.ident", "condition", "sex")])
observed <- observed[order(observed$orig.ident), ]
expected_cond <- c("O2", "Air", "Air", "O2")
expected_sex <- c("female", "male", "female", "male")
if (!identical(as.character(observed$condition), expected_cond) ||
    !identical(as.character(observed$sex), expected_sex)) {
  stop("Library design in the object does not match this script:\n",
       paste(utils::capture.output(print(observed)), collapse = "\n"),
       call. = FALSE)
}

#' Pseudo-bulk log2 CPM over a set of cells, one column per library.
pseudobulk_logcpm <- function(cell_mask) {
  pb <- vapply(lib_meta$lib, function(l) {
    Matrix::rowSums(counts[, cell_mask & lib == l, drop = FALSE])
  }, numeric(nrow(counts)))
  colnames(pb) <- lib_meta$label
  keep <- rowSums(pb > 0) == 4 & rowSums(pb) >= 20   # expressed in all four
  cpm <- t(t(pb[keep, ]) / colSums(pb[keep, ]) * 1e6)
  log2(cpm + 1)
}

# ---- all-cell pseudo-bulk ---------------------------------------------------
lcpm <- pseudobulk_logcpm(rep(TRUE, ncol(seu)))
var_rank <- names(sort(apply(lcpm, 1, var), decreasing = TRUE))
message(sprintf("  all-cell pseudo-bulk: %d genes expressed in all four libraries",
                nrow(lcpm)))

# ---- Fig S3A: PCA -----------------------------------------------------------
hv <- var_rank[seq_len(min(n_var_pca, length(var_rank)))]
pca <- prcomp(t(lcpm[hv, ]), center = TRUE, scale. = TRUE)
ve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pc <- merge(data.frame(pca$x[, 1:2], label = rownames(pca$x)), lib_meta,
            by = "label")
message(sprintf("  PCA variance explained: PC1 %.1f%%, PC2 %.1f%%, PC3 %.1f%%",
                ve[1], ve[2], ve[3]))

# PC1 should separate exposure and PC2 sex: that ordering is the panel's claim
pc1_by_exposure <- abs(diff(tapply(pc$PC1, pc$exposure, mean)))
pc1_by_sex <- abs(diff(tapply(pc$PC1, pc$sex, mean)))
expect_value("PC1 separates exposure more than sex",
             pc1_by_exposure > pc1_by_sex, TRUE)
expect_value("Fig S3A PC1 variance explained (%)", ve[1], 66.3,
             tol = 1.0, on_fail = "warn")
expect_value("Fig S3A PC2 variance explained (%)", ve[2], 18.1,
             tol = 1.0, on_fail = "warn")

ps3a <- ggplot(pc, aes(.data$PC1, .data$PC2,
                       colour = .data$exposure, shape = .data$sex)) +
  geom_hline(yintercept = 0, colour = "grey88") +
  geom_vline(xintercept = 0, colour = "grey88") +
  geom_point(size = 4) +
  geom_text(aes(label = .data$label), vjust = -1.1, size = 3,
            colour = "grey20", show.legend = FALSE) +
  scale_colour_manual(values = c("Room air" = "grey40",
                                 "Hyperoxia" = "#d6604d")) +
  scale_shape_manual(values = c("Female" = 16, "Male" = 17)) +
  scale_x_continuous(expand = expansion(mult = 0.35)) +
  scale_y_continuous(expand = expansion(mult = 0.35)) +
  labs(
    title = "PCA of all-cell pseudo-bulks",
    subtitle = sprintf("top %d variable genes", length(hv)),
    x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]),
    colour = "Exposure", shape = "Sex"
  ) +
  theme_classic(base_size = 11)
save_panel("FigS3A", ps3a,
  data = data.frame(
    library = pc$label, exposure = pc$exposure, sex = pc$sex,
    PC1 = pc$PC1, PC2 = pc$PC2,
    PC1_percent = round(ve[1], 2), PC2_percent = round(ve[2], 2)
  ),
  width = 5, height = 4.8, slug = "pca"
)

# ---- Fig S3B: library-library correlation ----------------------------------
hv_b <- var_rank[seq_len(min(n_var_heatmap, length(var_rank)))]
cm <- cor(lcpm[hv_b, ])
ord <- lib_meta$label[c(3, 2, 1, 4)]   # air F, air M, O2 F, O2 M
m <- melt(cm[ord, ord])
colnames(m) <- c("x", "y", "r")
m$x <- factor(m$x, levels = ord)
m$y <- factor(m$y, levels = rev(ord))

ps3b <- ggplot(m, aes(.data$x, .data$y, fill = .data$r)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = sprintf("%.3f", .data$r)), size = 3.4,
            colour = ifelse(m$r > 0.96, "white", "grey15")) +
  scale_fill_gradient(low = "#f7f7f7", high = "#2166ac",
                      limits = c(0.90, 1), name = "Pearson r") +
  labs(title = "Library similarity", x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.line = element_blank(), axis.ticks = element_blank()
  )
save_panel("FigS3B", ps3b, data = m, width = 5, height = 4.8,
           slug = "library_correlation")

# ---- Fig S3C: alveolar macrophage fold-change concordance ------------------
lcpm_am <- pseudobulk_logcpm(seu@meta.data$celltype == "Alveolar macrophages")
fc <- data.frame(
  gene = rownames(lcpm_am),
  female = lcpm_am[, "Hyperoxia, female"] - lcpm_am[, "Room air, female"],
  male   = lcpm_am[, "Hyperoxia, male"]   - lcpm_am[, "Room air, male"]
)
am_de <- as.data.frame(read_excel(deseq_xlsx, sheet = "Alveolar macrophages"))
sig <- am_de$gene[!is.na(am_de$padj) & am_de$padj < padj_de]
fc$DE <- fc$gene %in% sig

r_de <- cor(fc$female[fc$DE], fc$male[fc$DE])
top_expressed <- names(sort(rowMeans(lcpm_am), decreasing = TRUE))[1:2000]
r_top <- cor(fc$female[fc$gene %in% top_expressed],
             fc$male[fc$gene %in% top_expressed])
pct_same <- 100 * mean(sign(fc$female[fc$DE]) == sign(fc$male[fc$DE]))
lg <- fc[fc$gene == "Lgals3", ]
message(sprintf(paste("  AM: r(DE)=%.3f n=%d | r(top 2000)=%.3f |",
                      "same direction %.1f%% | Lgals3 F %+.3f M %+.3f"),
                r_de, sum(fc$DE), r_top, pct_same, lg$female, lg$male))

expect_value("Fig S3C r across DE genes", r_de, 0.95, tol = 0.02)
expect_value("Fig S3C genes in the DE set", sum(fc$DE), 361)
expect_value("Fig S3C percent same direction", pct_same, 97.2, tol = 1.0)
expect_value("Lgals3 induced in both sexes",
             lg$female > 0 && lg$male > 0, TRUE)

ps3c <- ggplot(fc, aes(.data$female, .data$male)) +
  geom_hline(yintercept = 0, colour = "grey88") +
  geom_vline(xintercept = 0, colour = "grey88") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey65") +
  geom_point(data = subset(fc, !fc$DE), colour = "grey80", alpha = 0.25,
             size = 0.6) +
  geom_point(data = subset(fc, fc$DE), colour = "#2166ac", alpha = 0.75,
             size = 1.2) +
  geom_point(data = lg, colour = "#d6604d", size = 3) +
  annotate("text", x = lg$female, y = lg$male, label = "Lgals3",
           colour = "#d6604d", hjust = -0.25, vjust = -0.6,
           fontface = "italic", size = 3.6) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.6, size = 3.3,
           label = sprintf("DE genes: r = %.2f (n = %d)\ntop 2,000 expressed: r = %.2f",
                           r_de, sum(fc$DE), r_top)) +
  labs(
    title = "Alveolar macrophages",
    subtitle = sprintf("blue, differentially expressed (adj. P < %s); red, Lgals3",
                       padj_de),
    x = "log2 fold change, females (hyperoxia / air)",
    y = "log2 fold change, males (hyperoxia / air)"
  ) +
  theme_classic(base_size = 11)
save_panel("FigS3C", ps3c,
  data = fc[fc$DE, c("gene", "female", "male")],
  width = 5, height = 4.8, slug = "am_foldchange_concordance"
)

message("11_sex_concordance.R complete")
