# Shared plotting and modelling helpers.
#
# Ported verbatim in behaviour from the preamble of
# code/figures_for_paper_dblt_removed.R in the analysis repository, with data
# columns referenced through the .data pronoun so the functions pass static
# checks. Any change here changes published figures: edit with care.

suppressMessages({
  library(DESeq2)
  library(dplyr)
  library(fgsea)
  library(ggplot2)
  library(pheatmap)
  library(Seurat)
  library(speckle)
})

# ---- differential cell-type abundance (propeller) ---------------------------
run_differential_frequency <- function(sub) {
  propeller(
    clusters = sub$celltype,
    sample = sub$orig.ident,
    group = sub$condition
  )
}

plot_differential_frequency <- function(res) {
  ggplot(res, aes(
    .data$BaselineProp.Freq,
    (.data$PropMean.O2 / .data$PropMean.Air) * 100,
    size = -log10(.data$P.Value),
    color = log2(.data$PropMean.O2 / .data$PropMean.Air)
  )) +
    labs(
      y = "Percent [compared to Air]",
      x = "Baseline frequency",
      color = "log2(O2/Air)"
    ) +
    geom_hline(yintercept = 100, linetype = "dashed") +
    geom_point() +
    scale_color_gradient2(low = "blue", mid = "grey", high = "red") +
    ggrepel::geom_text_repel(aes(label = .data$BaselineProp.clusters),
      color = "black"
    ) +
    theme_classic()
}

#' Proportion of one cell type per library, by exposure (Fig 2B, Fig 5D).
plot_freq <- function(cell, seu, ylab = NULL) {
  fr <- table(seu@meta.data$orig.ident[seu@meta.data$celltype == cell]) /
    table(seu@meta.data$orig.ident)
  meta <- seu@meta.data[match(names(fr), seu@meta.data$orig.ident), ]
  subm <- data.frame(fr, meta)
  ggplot(subm, aes(.data$condition, .data$Freq * 100, color = .data$condition)) +
    labs(
      y = if (is.null(ylab)) paste(cell, "frequency [%]") else ylab,
      x = "Condition"
    ) +
    geom_boxplot() +
    geom_point() +
    scale_color_manual(values = c("black", "red")) +
    theme_classic()
}

#' Tidy frame behind plot_freq(), for the panel _data.csv.
freq_table <- function(cell, seu) {
  fr <- table(seu@meta.data$orig.ident[seu@meta.data$celltype == cell]) /
    table(seu@meta.data$orig.ident)
  meta <- seu@meta.data[match(names(fr), seu@meta.data$orig.ident), ]
  data.frame(
    cell_type = cell, library = names(fr),
    percent = as.numeric(fr) * 100, condition = meta$condition
  )
}

# ---- per-cell-type pseudo-bulk differential expression (DESeq2) -------------
# Counts are aggregated per library (4 pooled libraries, 2 per condition);
# genes detected in fewer than 2 libraries are dropped.
run_de <- function(seu) {
  by_type <- split(colnames(seu), seu@meta.data$celltype)
  out <- lapply(by_type, function(cells) {
    x <- seu[, cells]
    by_sample <- split(seq_len(ncol(x)), x@meta.data$sample)
    sums <- do.call(cbind, lapply(by_sample, function(k) {
      rowSums(x@assays$SCT@counts[, k])
    }))
    sums <- sums[rowSums(sums) > 0, ]
    detected <- apply(sums, 1, function(y) sum(y > 0))
    sums <- sums[which(detected >= 2), ]
    treat <- seu@meta.data$condition[
      match(colnames(sums), seu@meta.data$sample)
    ]
    dds <- DESeqDataSetFromMatrix(
      countData = sums,
      colData = data.frame(treat),
      design = ~treat
    )
    dds <- DESeq(dds, quiet = TRUE)
    final <- results(dds)
    final$gene <- rownames(final)
    final
  })
  names(out) <- names(by_type)
  out
}

#' Counts of up- and down-regulated genes per cell type (Fig 3A).
de_counts <- function(deseq, padj = 0.25) {
  celltypes <- names(deseq)
  n_up <- n_down <- numeric(length(celltypes))
  for (i in seq_along(celltypes)) {
    res <- deseq[[i]]
    sig <- !is.na(res$padj) & res$padj < padj
    n_up[i] <- sum(sig & res$log2FoldChange > 0)
    n_down[i] <- sum(sig & res$log2FoldChange < 0)
  }
  data.frame(cell_type = celltypes, n_up = n_up, n_down = n_down,
             n_total = n_up + n_down)
}

plot_de <- function(deseq, padj = 0.25) {
  cnt <- de_counts(deseq, padj)
  ord <- order(cnt$n_total, decreasing = TRUE)
  df <- rbind(
    data.frame(cell_type = cnt$cell_type[ord], direction = "Up",
               count = cnt$n_up[ord]),
    data.frame(cell_type = cnt$cell_type[ord], direction = "Down",
               count = -cnt$n_down[ord])
  )
  df$cell_type <- factor(df$cell_type, levels = rev(cnt$cell_type[ord]))
  ggplot(df, aes(y = .data$cell_type, x = .data$count,
                 fill = .data$direction)) +
    geom_col() +
    geom_vline(xintercept = 0, color = "black") +
    scale_fill_manual(values = c("Up" = "red", "Down" = "blue")) +
    labs(
      title = "Significant DE genes per cell type",
      x = paste0("Number of DE genes (padj < ", padj, ")"), y = "Cell type"
    ) +
    theme_classic()
}

# ---- volcano for one cell type's DESeq2 table (Fig 4C/4D) -------------------
create_volcano <- function(subm, n_label = 20) {
  subm$gene <- rownames(subm)
  subm <- na.omit(subm)
  subm <- subm[sort.list(subm$padj), ]
  subm$pvalue[subm$pvalue < 1e-50] <- 1e-50
  subm$col <- "grey"
  subm$col[subm$padj < 0.25 & subm$log2FoldChange > 0] <- "red"
  subm$col[subm$padj < 0.25 & subm$log2FoldChange < 0] <- "blue"
  ggplot(subm, aes(.data$log2FoldChange, -log10(.data$pvalue))) +
    labs(x = "log2 fold change [O2/air]") +
    geom_point(aes(color = .data$col)) +
    scale_color_identity() +
    ggrepel::geom_text_repel(
      data = utils::head(subm, n_label),
      aes(label = .data$gene), color = "black"
    ) +
    theme_classic()
}

# ---- pseudo-bulk expression of one gene in one cell type (Fig 4E/4G) --------
create_box_plot <- function(gene, cell, seu) {
  sub <- seu[, seu@meta.data$celltype == cell]
  by_sample <- split(seq_len(ncol(sub)), sub@meta.data$sample)
  sums <- do.call(cbind, lapply(by_sample, function(k) {
    rowSums(sub@assays$SCT@counts[, k])
  }))
  sums <- sums[rowSums(sums) > 0, ]
  detected <- apply(sums, 1, function(y) sum(y > 0))
  sums <- sums[which(detected >= 2), ]
  treat <- sub@meta.data$condition[
    match(colnames(sums), sub@meta.data$sample)
  ]
  dds <- DESeqDataSetFromMatrix(sums, data.frame(treat), ~treat)
  dds <- DESeq(dds, quiet = TRUE)
  subm <- data.frame(
    expression = log(counts(dds, normalized = TRUE)[gene, ]),
    sub@meta.data[match(colnames(sums), sub@meta.data$sample), ]
  )
  p <- ggplot(subm, aes(.data$condition, .data$expression,
                        color = .data$condition)) +
    labs(title = cell, y = paste(gene, "expression levels")) +
    geom_boxplot() +
    geom_point() +
    scale_color_manual(values = c("grey", "red")) +
    theme_classic()
  list(plot = p, data = data.frame(
    gene = gene, cell_type = cell, library = rownames(subm),
    log_normalized = subm$expression, condition = subm$condition
  ))
}

# ---- Hallmark GSEA over per-cell-type fold changes (Fig 3B/3C) --------------
run_enrich <- function(deseq, paths) {
  out <- lapply(deseq, function(x) {
    x <- na.omit(x)
    stats <- x$log2FoldChange
    names(stats) <- rownames(x)
    fgsea(pathways = paths, stats = stats, minSize = 5)
  })
  names(out) <- names(deseq)
  out
}

nes_matrix <- function(enrich, padj = 0.25) {
  all_paths <- unique(unlist(lapply(enrich, function(x) {
    x$pathway[x$padj < padj]
  })))
  m <- do.call(cbind, lapply(enrich, function(x) {
    x$NES[match(all_paths, x$pathway)]
  }))
  rownames(m) <- gsub("HALLMARK_", "", all_paths, fixed = TRUE)
  m[is.na(m)] <- 0
  m
}

plot_enrich <- function(enrich, padj = 0.25) {
  m <- nes_matrix(enrich, padj)
  pheatmap(m,
    scale = "none", breaks = seq(-2, 2, length = 100),
    color = colorRampPalette(c("blue", "white", "red"))(100),
    fontsize_row = 7, main = "Hallmark NES Across Cell Types",
    silent = TRUE
  )
}
