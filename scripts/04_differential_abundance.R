# 04 - Differential cell-type abundance (Figure 2A/B; Fig S7 sex).
# Source: figures_for_paper_dblt_removed.R (propeller + frequency boxplots)
# Fig 2C (cross-study validation vs Hurskainen) is produced in
#   scripts/aux_hurskainen_frequency.R  (see data/README for GSE151974).
# Input : Seurat_object_scrublet.RData

source(file.path("config", "paths.R"))
source(file.path("scripts", "utils", "functions.R"))
source(file.path("scripts", "utils", "load_annotate.R"))   # -> seu
library(speckle); library(ggplot2)

# ---- Fig 2A: propeller differential abundance -------------------------------
res <- propeller(clusters = seu$celltype, sample = seu$orig.ident, group = seu$condition)
write.csv(res, file.path(tab_dir, "propeller_abundance.csv"), row.names = FALSE)
plot_differential_frequency(res)
ggsave(file.path(fig_dir, "Fig2A_frequency_volcano.pdf"), height = 5, width = 6)

# ---- Fig 2B: frequency boxplots for the significantly increased populations -
plot_freq <- function(cell){
  fr <- table(seu@meta.data$orig.ident[seu@meta.data$celltype == cell]) /
        table(seu@meta.data$orig.ident)
  subm <- data.frame(fr, seu@meta.data[match(names(fr), seu@meta.data$orig.ident), ])
  ggplot(subm, aes(condition, Freq * 100, color = condition)) +
    labs(y = paste(cell, "frequency [%]"), x = "Condition") +
    geom_boxplot() + geom_point() +
    scale_color_manual(values = c("black", "red")) + theme_classic()
}
p <- gridExtra::grid.arrange(
  plot_freq("NK cells"), plot_freq("Dendritic cells"),
  plot_freq("Interstitial macrophages"), nrow = 1)
ggsave(p, filename = file.path(fig_dir, "Fig2B_frequency_boxplots.pdf"), height = 4, width = 6)

# ---- Fig S7: abundance split by sex -----------------------------------------
res_sex <- propeller(clusters = seu$celltype, sample = seu$orig.ident, group = seu$sex)
write.csv(res_sex, file.path(tab_dir, "propeller_by_sex.csv"), row.names = FALSE)
