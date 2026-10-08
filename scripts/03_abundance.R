# 03 - Differential cell-type abundance.
#
# Ported from code/figures_for_paper_dblt_removed.R (lines 344-381).
#
# Input   frozen Seurat object
# Outputs Fig 2A  propeller differential-abundance plot
#         Fig 2B  proportions of the three significantly expanded populations
#         Table S2  propeller results for all nine populations
#
# propeller fits a linear model with empirical Bayes moderation to
# logit-transformed proportions, treating the four libraries as the
# replication unit.

source(file.path("config", "paths.R"))
source(file.path("R", "figure_io.R"))
source(file.path("R", "plots.R"))

seu <- load_annotated()

# ---- Table S2 / Fig 2A ------------------------------------------------------
res <- run_differential_frequency(seu)
res$cell_type <- rownames(res)
save_table("TableS2_abundance", res, path = abundance_csv)

sig_up <- res$cell_type[res$P.Value < 0.05 &
  res$PropMean.O2 > res$PropMean.Air]
expect_value(
  "populations significantly increased in hyperoxia",
  sort(sig_up),
  sort(c("NK cells", "Dendritic cells", "Interstitial macrophages"))
)

save_panel("Fig2A", plot_differential_frequency(res),
  data = res[, c("cell_type", "BaselineProp.Freq", "PropMean.Air",
                 "PropMean.O2", "P.Value", "FDR")],
  width = 6, height = 5, slug = "abundance"
)

# ---- Fig 2B: the three expanded populations ---------------------------------
# One manuscript panel; the three boxplots are assembled side by side as in
# the submitted figure.
cells_2b <- c("NK cells", "Dendritic cells", "Interstitial macrophages")
p2b <- gridExtra::grid.arrange(
  grobs = lapply(cells_2b, plot_freq, seu = seu), nrow = 1
)
save_panel("Fig2B", p2b,
  data = do.call(rbind, lapply(cells_2b, freq_table, seu = seu)),
  width = 6, height = 4, slug = "frequency_boxplots"
)

message("03_abundance.R complete")
