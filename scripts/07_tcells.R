# 07 - T-cell subclustering and the gamma-delta expansion.
#
# Ported from code/figures_for_paper_dblt_removed.R (lines 616-785).
#
# Input   frozen Seurat object
# Outputs Fig 5A  T-cell subtype UMAP
#         Fig 5B  subtype marker dot plot
#         Fig 5C  differential abundance across subtypes
#         Fig 5D  gamma-delta proportion by exposure
#         Fig 5F  IL-17 module score on the T-cell UMAP
#         Fig 5G  IL-17 module within gamma-delta cells, by exposure
#         (Fig 5E is produced by 10_external.R from the Hurskainen data)
#
# REPRODUCIBILITY WARNING ----------------------------------------------------
# Unlike the macrophage UMAPs in script 06, the subtypes here are *derived* by
# re-clustering at analysis time, and the manuscript's central T-cell result
# (gamma-delta 7.9% -> 16.5%) depends on which cells land in which cluster.
# The working script attached labels to hardcoded cluster numbers
# ("0" = CD4, "1" = CD8, "2" = dg, ...), which silently mislabels every
# population if a different Seurat or igraph version renumbers the clusters.
#
# This port therefore assigns each label by the marker the cluster actually
# expresses most highly, then asserts the resulting subtype sizes against the
# manuscript. A renumbering is handled automatically; a genuinely different
# partition stops the script instead of producing relabelled results.

source(file.path("config", "paths.R"))
source(file.path("R", "figure_io.R"))
source(file.path("R", "plots.R"))
suppressMessages({
  library(Seurat)
  library(DESeq2)
  library(ggplot2)
})

seu <- load_annotated()

# ---- subcluster the T-cell compartment -------------------------------------
sub <- seu[, seu@meta.data$celltype == "T cells"]
expect_value("T cells extracted", ncol(sub), 6935)

DefaultAssay(sub) <- "RNA"
sub <- NormalizeData(sub, verbose = FALSE)
sub <- FindVariableFeatures(sub, nfeatures = 3000, verbose = FALSE)
sub <- ScaleData(sub,
  vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), verbose = FALSE
)
sub <- RunPCA(sub, npcs = 20, verbose = FALSE)
sub <- FindNeighbors(sub, dims = 1:20, verbose = FALSE)
sub <- RunUMAP(sub, dims = 1:20, verbose = FALSE)
sub <- FindClusters(sub, resolution = 0.1, verbose = FALSE)

expect_value("T-cell subclusters", length(unique(sub$seurat_clusters)), 5)

# ---- label clusters by their defining marker, not by cluster number --------
# Each marker is essentially absent from the other four subtypes (Fig 5B), so
# "the cluster with the highest mean expression" is an unambiguous assignment.
tcell_markers <- c(
  "CD4 T cells" = "Lef1",
  "CD8 T cells" = "Cd8b1",
  "dg T cells"  = "Trdc",
  "Th2 cells"   = "Il1rl1",
  "Tregs"       = "Foxp3"
)
expr <- GetAssayData(sub, assay = "RNA", layer = "data")
cl <- as.character(sub$seurat_clusters)
assignment <- vapply(tcell_markers, function(g) {
  means <- vapply(split(seq_along(cl), cl), function(idx) {
    mean(expr[g, idx])
  }, numeric(1))
  names(which.max(means))
}, character(1))

if (anyDuplicated(assignment)) {
  stop("Marker-based labelling is ambiguous: two subtypes claim the same ",
       "cluster.\n  ", paste(sprintf("%s -> cluster %s", names(assignment),
                                     assignment), collapse = "\n  "),
       "\nThe T-cell partition differs from the one used for the paper.")
}
message("  cluster assignment by marker:")
for (i in seq_along(assignment)) {
  message(sprintf("    %-12s <- cluster %s (%s)", names(assignment)[i],
                  assignment[i], tcell_markers[i]))
}

label_of <- setNames(names(assignment), assignment)
sub@meta.data$celltype <- unname(label_of[cl])

# ---- the subtype sizes reported in the Results -----------------------------
sizes <- table(sub@meta.data$celltype)
expect_value("naive CD4 T cells", as.integer(sizes[["CD4 T cells"]]), 3454)
expect_value("CD8 T cells",       as.integer(sizes[["CD8 T cells"]]), 1635)
expect_value("gamma-delta T cells", as.integer(sizes[["dg T cells"]]), 853)
expect_value("Th2 cells",         as.integer(sizes[["Th2 cells"]]), 638)
expect_value("regulatory T cells", as.integer(sizes[["Tregs"]]), 355)

# ---- Fig 5A / 5B ------------------------------------------------------------
umap_xy <- as.data.frame(Embeddings(sub, "umap"))
names(umap_xy) <- c("umap_1", "umap_2")
umap_xy$subtype <- sub@meta.data$celltype

save_panel("Fig5A", DimPlot(sub, group.by = "celltype", label = TRUE),
  data = umap_xy, width = 8, height = 7, slug = "umap_tcells"
)

p5b <- DotPlot(sub, features = unname(tcell_markers), group.by = "celltype") +
  labs(x = "Marker gene", y = "Cell type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_panel("Fig5B", p5b, data = p5b$data, width = 6, height = 5,
           slug = "dotplot_tcell_markers")

# ---- Fig 5C: differential abundance across subtypes ------------------------
res <- run_differential_frequency(sub)
res$subtype <- rownames(res)
save_panel("Fig5C", plot_differential_frequency(res),
  data = res[, c("subtype", "BaselineProp.Freq", "PropMean.Air",
                 "PropMean.O2", "P.Value", "FDR")],
  width = 6, height = 5, slug = "diff_freq_tcell"
)

gd <- res[res$subtype == "dg T cells", ]
expect_value("gamma-delta proportion in air (%)",
             100 * gd$PropMean.Air, 7.9, tol = 0.15)
expect_value("gamma-delta proportion in hyperoxia (%)",
             100 * gd$PropMean.O2, 16.5, tol = 0.15)
expect_value("gamma-delta abundance P value", gd$P.Value, 1.4e-3, tol = 2e-4)

# ---- Fig 5D: gamma-delta proportion by exposure ----------------------------
save_panel("Fig5D",
  plot_freq("dg T cells", sub, ylab = "dg T cell frequency [% T cells]"),
  data = freq_table("dg T cells", sub), width = 3, height = 5,
  slug = "gd_frequency"
)

# ---- Fig 5F: IL-17 module on the UMAP --------------------------------------
il17_genes <- c("Rorc", "Il23r", "Il1r1", "Il17a", "Il17f")
DefaultAssay(sub) <- "RNA"
sub <- AddModuleScore(sub,
  features = list(il17_genes), name = "IL17_module", assay = "RNA"
)
p5f <- FeaturePlot(sub,
  features = "IL17_module1", min.cutoff = "q05", max.cutoff = "q95"
) +
  labs(title = "IL-17 module score") +
  theme_classic()
save_panel("Fig5F", p5f,
  data = data.frame(umap_xy, il17_module = sub$IL17_module1),
  width = 8, height = 7, slug = "il17_module_umap"
)

# ---- Fig 5G: IL-17 module within gamma-delta cells -------------------------
# Pseudo-bulk per library, so the test is on four replicates rather than cells.
module_box <- function(genes, subtype) {
  ss <- sub[, sub@meta.data$celltype == subtype]
  by_sample <- split(seq_len(ncol(ss)), ss@meta.data$sample)
  sums <- do.call(cbind, lapply(by_sample, function(k) {
    rowSums(ss@assays$SCT@counts[, k])
  }))
  sums <- sums[rowSums(sums) > 0, ]
  treat <- ss@meta.data$condition[match(colnames(sums), ss@meta.data$sample)]
  dds <- DESeqDataSetFromMatrix(sums, data.frame(treat), ~treat)
  dds <- DESeq(dds, quiet = TRUE)
  norm <- counts(dds, normalized = TRUE)
  genes <- intersect(genes, rownames(norm))
  module <- colMeans(log1p(norm[genes, , drop = FALSE]))
  meta <- ss@meta.data[match(names(module), ss@meta.data$sample), ]
  subm <- data.frame(module = module, meta)
  p <- ggplot(subm, aes(.data$condition, .data$module,
                        color = .data$condition)) +
    labs(title = subtype, y = "IL-17 module score (pseudobulk)") +
    geom_boxplot() +
    geom_point() +
    ggpubr::stat_compare_means() +
    scale_color_manual(values = c("grey", "red")) +
    theme_classic()
  list(plot = p, data = data.frame(
    subtype = subtype, library = rownames(subm),
    il17_module_pseudobulk = subm$module, condition = subm$condition
  ))
}
gd_mod <- module_box(il17_genes, "dg T cells")
save_panel("Fig5G", gd_mod$plot, data = gd_mod$data,
           width = 3, height = 5, slug = "il17_module_gd")

message("07_tcells.R complete")
