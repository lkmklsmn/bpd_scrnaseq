# 06 - Alveolar vs interstitial macrophage programs.
#
# Ported from code/figures_for_paper_dblt_removed.R (lines 508-612).
#
# Input   frozen Seurat object, de_by_celltype.rds (from 04)
# Outputs Fig 4A  alveolar macrophage UMAP by exposure
#         Fig 4B  interstitial macrophage UMAP by exposure
#         Fig 4C  alveolar macrophage volcano
#         Fig 4D  interstitial macrophage volcano
#         Fig 4E  representative genes (Fabp1, Cd63 in AM; Il1b, Cdkn1a in IM)
#         Fig 4F  interferon running-enrichment, both populations
#         Fig 4G  Ifi30, oppositely regulated in the two populations
#
# NOTE Re-embedding each population (4A/4B) uses the RNA assay with cell-cycle
# and mitochondrial content regressed out, as in the working script. The UMAP
# is descriptive only; no result depends on its coordinates, and nothing is
# re-clustered, so version-dependent embedding drift cannot change a
# conclusion here.

source(file.path("config", "paths.R"))
source(file.path("R", "figure_io.R"))
source(file.path("R", "plots.R"))
suppressMessages({
  library(Seurat)
  library(msigdbr)
  library(fgsea)
  library(ggplot2)
})

seu <- load_annotated()
des_deseq <- readRDS(deseq_rds)

#' Re-embed one population on its own variable genes.
reembed <- function(cell, seu) {
  sub <- seu[, seu@meta.data$celltype == cell]
  DefaultAssay(sub) <- "RNA"
  sub <- NormalizeData(sub, verbose = FALSE)
  sub <- FindVariableFeatures(sub, nfeatures = 3000, verbose = FALSE)
  sub <- ScaleData(sub,
    vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"),
    verbose = FALSE
  )
  sub <- RunPCA(sub, npcs = 50, verbose = FALSE)
  sub <- FindNeighbors(sub, dims = 1:20, verbose = FALSE)
  RunUMAP(sub, dims = 1:20, verbose = FALSE)
}

umap_frame <- function(sub) {
  xy <- as.data.frame(Embeddings(sub, "umap"))
  names(xy) <- c("umap_1", "umap_2")
  xy$condition <- sub@meta.data$condition
  xy
}

# ---- Fig 4A / 4B: condition-coloured UMAPs ---------------------------------
am_embed <- reembed("Alveolar macrophages", seu)
save_panel("Fig4A", DimPlot(am_embed, group.by = "condition"),
  data = umap_frame(am_embed), width = 8, height = 7, slug = "umap_alv_mac"
)

im_embed <- reembed("Interstitial macrophages", seu)
save_panel("Fig4B", DimPlot(im_embed, group.by = "condition"),
  data = umap_frame(im_embed), width = 8, height = 7, slug = "umap_int_mac"
)

# ---- Fig 4C / 4D: volcanoes -------------------------------------------------
volcano_frame <- function(res) {
  d <- as.data.frame(res)
  d$gene <- rownames(d)
  d <- na.omit(d)
  d[order(d$padj), c("gene", "log2FoldChange", "pvalue", "padj", "baseMean")]
}

am_de <- des_deseq[["Alveolar macrophages"]]
im_de <- des_deseq[["Interstitial macrophages"]]

save_panel("Fig4C", create_volcano(am_de), data = volcano_frame(am_de),
           width = 4, height = 5, slug = "volcano_alv_mac")
save_panel("Fig4D", create_volcano(im_de), data = volcano_frame(im_de),
           width = 4, height = 5, slug = "volcano_int_mac")

# the two fold changes quoted in the Results
am_fc <- function(g) am_de$log2FoldChange[match(g, am_de$gene)]
expect_value("Fabp1 log2FC in alveolar macrophages", am_fc("Fabp1"),
             -3.52, tol = 0.02)
expect_value("Il1b log2FC in interstitial macrophages",
             im_de$log2FoldChange[match("Il1b", im_de$gene)], -1.62,
             tol = 0.02)

# ---- Fig 4E: representative genes ------------------------------------------
reps <- list(
  c("Fabp1", "Alveolar macrophages"),
  c("Cd63", "Alveolar macrophages"),
  c("Il1b", "Interstitial macrophages"),
  c("Cdkn1a", "Interstitial macrophages")
)
boxes <- lapply(reps, function(x) create_box_plot(x[1], x[2], seu))
p4e <- gridExtra::grid.arrange(grobs = lapply(boxes, `[[`, "plot"), nrow = 1)
save_panel("Fig4E", p4e,
  data = do.call(rbind, lapply(boxes, `[[`, "data")),
  width = 10, height = 5, slug = "representative_genes"
)

# ---- Fig 4F: interferon running enrichment ---------------------------------
hallmark <- msigdbr(species = "Mus musculus", collection = "H")
paths <- split(hallmark$gene_symbol, hallmark$gs_name)
ifn_genes <- union(
  paths[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]],
  paths[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]]
)

ranked <- function(res) {
  d <- na.omit(as.data.frame(res))
  stats <- d$log2FoldChange
  names(stats) <- d$gene
  stats
}
p_am <- fgsea::plotEnrichment(pathway = ifn_genes, stats = ranked(am_de)) +
  labs(title = "Alveolar macrophages") + theme_classic()
p_im <- fgsea::plotEnrichment(pathway = ifn_genes, stats = ranked(im_de)) +
  labs(title = "Interstitial macrophages") + theme_classic()

save_panel("Fig4F", gridExtra::grid.arrange(p_am, p_im, ncol = 1),
  data = data.frame(
    cell_type = c("Alveolar macrophages", "Interstitial macrophages"),
    n_interferon_genes_ranked = c(
      sum(names(ranked(am_de)) %in% ifn_genes),
      sum(names(ranked(im_de)) %in% ifn_genes)
    ),
    n_genes_ranked = c(length(ranked(am_de)), length(ranked(im_de)))
  ),
  width = 4, height = 6, slug = "interferon_enrichment"
)

# ---- Fig 4G: Ifi30, opposite direction in the two populations --------------
ifi30 <- lapply(c("Alveolar macrophages", "Interstitial macrophages"),
                function(ct) create_box_plot("Ifi30", ct, seu))
save_panel("Fig4G", gridExtra::grid.arrange(grobs = lapply(ifi30, `[[`, "plot"),
                                            nrow = 1),
  data = do.call(rbind, lapply(ifi30, `[[`, "data")),
  width = 5, height = 5, slug = "ifi30_boxplot"
)

expect_value("Ifi30 log2FC in alveolar macrophages", am_fc("Ifi30"),
             -0.70, tol = 0.02)

message("06_macrophages.R complete")
