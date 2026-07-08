# Shared helper functions for the discovery-data pipeline (scripts 03-08).
# Extracted from figures_for_paper_dblt_removed.R. `seu`, `paths` must be in scope.
# NOTE: the cross-dataset script (10) deliberately does NOT use these helpers, so
# per-dataset models/plots can be fine-tuned independently.

suppressMessages({
  library(DESeq2); library(dplyr); library(fgsea); library(ggplot2)
  library(pheatmap); library(Seurat); library(speckle)
})

# ---- volcano for a per-cell-type DESeq2 result table ------------------------
create_volcano2 <- function(subm){
  subm$gene <- rownames(subm)
  subm <- na.omit(subm)
  subm <- subm[sort.list(subm$padj), ]
  subm$pvalue[subm$pvalue < 1e-50] <- 1e-50
  subm$col <- "grey"
  subm$col[subm$padj < 0.25 & subm$log2FoldChange > 0] <- "red"
  subm$col[subm$padj < 0.25 & subm$log2FoldChange < 0] <- "blue"
  ggplot(subm, aes(log2FoldChange, -log10(pvalue))) +
    labs(x = "log2 fold change [O2/air]") +
    geom_point(aes(color = col)) + scale_color_identity() +
    ggrepel::geom_text_repel(data = subm[1:20, ], aes(label = gene), color = "black") +
    theme_classic()
}

# ---- pseudobulk boxplot of one gene in one cell type (DESeq2-normalized) -----
create_box_plot <- function(gene, cell, seu){
  sub <- seu[, seu$celltype == cell]
  bsplit <- split(1:ncol(sub), sub@meta.data$sample)
  sums <- do.call(cbind, lapply(bsplit, function(k) rowSums(sub@assays$SCT@counts[, k])))
  sums <- sums[rowSums(sums) > 0, ]
  detected <- apply(sums, 1, function(y) sum(y > 0))
  sums <- sums[which(detected >= 2), ]
  treat <- sub@meta.data$condition[match(colnames(sums), sub@meta.data$sample)]
  dds <- DESeqDataSetFromMatrix(sums, data.frame(treat), ~ treat); dds <- DESeq(dds)
  subm <- data.frame(gene = log(counts(dds, normalized = TRUE)[gene, ]),
                     sub@meta.data[match(colnames(sums), sub@meta.data$sample), ])
  ggplot(subm, aes(condition, gene, color = condition)) +
    labs(title = cell, y = paste(gene, "expression levels")) +
    geom_boxplot() + geom_point() +
    scale_color_manual(values = c("grey", "red")) + theme_classic()
}

# ---- differential abundance (propeller) -------------------------------------
run_differential_frequency <- function(sub)
  propeller(clusters = sub$celltype, sample = sub$orig.ident, group = sub$condition)

plot_differential_frequency <- function(res)
  ggplot(res, aes(BaselineProp.Freq, (PropMean.O2 / PropMean.Air) * 100,
                  size = -log10(P.Value), color = log2(PropMean.O2 / PropMean.Air))) +
    labs(y = "Percent [compared to air]", x = "Baseline frequency", color = "log2(O2/air)") +
    geom_hline(yintercept = 100, linetype = "dashed") + geom_point() +
    scale_color_gradient2(low = "blue", mid = "grey", high = "red") +
    ggrepel::geom_text_repel(aes(label = BaselineProp.clusters), color = "black") +
    theme_classic()

# ---- per-cell-type pseudobulk DE (DESeq2) -----------------------------------
run_de <- function(sub){
  asplit <- split(colnames(sub), sub$celltype)
  deseq <- lapply(asplit, function(x){
    xs <- sub[, x]
    bsplit <- split(1:ncol(xs), xs@meta.data$sample)
    sums <- do.call(cbind, lapply(bsplit, function(k) rowSums(xs@assays$SCT@counts[, k])))
    sums <- sums[rowSums(sums) > 0, ]
    detected <- apply(sums, 1, function(y) sum(y > 0))
    sums <- sums[which(detected >= 2), ]
    treat <- sub@meta.data$condition[match(colnames(sums), sub@meta.data$sample)]
    dds <- DESeqDataSetFromMatrix(sums, data.frame(treat), ~ treat); dds <- DESeq(dds)
    final <- results(dds); final$gene <- rownames(final); final
  })
  names(deseq) <- names(asplit); deseq
}

plot_de <- function(deseq){
  celltypes <- names(deseq)
  n_up <- n_down <- numeric(length(celltypes))
  for (i in seq_along(celltypes)) {
    res <- deseq[[i]]; sig <- !is.na(res$padj) & res$padj < 0.25
    n_up[i]   <- sum(sig & res$log2FoldChange > 0)
    n_down[i] <- sum(sig & res$log2FoldChange < 0)
  }
  total_de <- n_up + n_down; names(total_de) <- celltypes
  ord <- order(total_de, decreasing = TRUE)
  df <- rbind(
    data.frame(celltype = celltypes[ord], direction = "Up",   count =  n_up[ord]),
    data.frame(celltype = celltypes[ord], direction = "Down", count = -n_down[ord]))
  df$celltype <- factor(df$celltype, levels = rev(names(total_de)[ord]))
  ggplot(df, aes(y = celltype, x = count, fill = direction)) +
    geom_col() + geom_vline(xintercept = 0, color = "black") +
    scale_fill_manual(values = c("Up" = "red", "Down" = "blue")) +
    labs(title = "Significant DE genes per cell type",
         x = "Number of DE genes (padj < 0.25)", y = "Cell type") +
    theme_classic()
}

# ---- GSEA (Hallmark) over per-cell-type fold changes ------------------------
# requires `paths` (named list of gene sets) in scope
run_enrich <- function(deseq, paths){
  out <- lapply(deseq, function(x){
    x <- na.omit(x); stats <- x$log2FoldChange; names(stats) <- rownames(x)
    fgsea(pathways = paths, stats = stats, minSize = 5)
  })
  names(out) <- names(deseq); out
}

plot_enrich <- function(enrich){
  all_paths <- unique(unlist(lapply(enrich, function(x) x$pathway[x$padj < 0.25])))
  nes_mat <- do.call(cbind, lapply(enrich, function(x) x$NES[match(all_paths, x$pathway)]))
  rownames(nes_mat) <- gsub("HALLMARK_", "", fixed = TRUE, all_paths)
  nes_mat[is.na(nes_mat)] <- 0
  pheatmap(nes_mat, scale = "none", breaks = seq(-2, 2, length = 100),
           color = colorRampPalette(c("blue", "white", "red"))(100),
           fontsize_row = 7, main = "Hallmark NES across cell types")
}
