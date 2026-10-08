# 09 - Galectin-3 across mouse and human BPD datasets.
#
# Ported from code/lgals3_forest_figure_BACKUP_pre_sex.R in the analysis
# repository. That backup, NOT the current lgals3_forest_figure.R, is the code
# behind the submitted Figure 7: the current version adds a sex covariate and
# a per-infant random intercept to the GSE220135 model, which changes that
# estimate from +0.33 (P = 0.095) to +0.27 (P = 0.23). The published model is
# reproduced here deliberately.
#
# Inputs  frozen Seurat object; GSE151974; LungMAP LMEX0000004400; GSE32472;
#         GSE220135 (counts fetched from GEO, metadata frozen in data/frozen)
# Outputs Fig 7A   mouse alveolar macrophages, this study
#         Fig 7B   mouse alveolar macrophages, Hurskainen et al.
#         Fig 7C   human alveolar macrophages, LungMAP, by disease state
#         Fig 7D1  human blood, GSE32472 (microarray)
#         Fig 7D2  human blood, GSE220135 (bulk RNA-seq)
#         Fig 7E   forest plot of all five effect sizes
#         Table S6 harmonised effect sizes
#
# Each dataset is modelled INDEPENDENTLY, with no shared helper functions, so
# that per-dataset choices stay visible and tunable. Effect sizes are
# harmonised as a log2 fold change between the disease and control state.

source(file.path("config", "paths.R"))
source(file.path("R", "figure_io.R"))
suppressMessages({
  library(data.table)
  library(ggplot2)
  library(Matrix)
  library(glmmTMB)
  library(lme4)
  library(lmerTest)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(Seurat)
  library(DESeq2)
})

need <- function(path, accession) {
  if (!file.exists(path)) {
    stop("Missing external data for ", accession, ":\n  ", path,
         "\nRun scripts/00_download_external.R, or see data/README.md.",
         call. = FALSE)
  }
  path
}

res <- data.frame()
boxplot_theme <- function(p, title, ylab) {
  p + geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.12, size = 1.8) +
    scale_color_manual(values = c("grey40", "#d6604d")) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(title = title, x = NULL, y = ylab)
}

# =============================================================================
# (1) Mouse alveolar macrophages, this study - DESeq2 pseudo-bulk
# =============================================================================
seu <- load_annotated()
am <- seu[, seu@meta.data$celltype == "Alveolar macrophages"]
s <- as.character(am$orig.ident)
u <- sort(unique(s))
pb <- do.call(cbind, lapply(u, function(x) {
  Matrix::rowSums(am@assays$SCT@counts[, s == x, drop = FALSE])
}))
colnames(pb) <- u
rownames(pb) <- rownames(am@assays$SCT@counts)
grp <- factor(am$condition[match(u, am$orig.ident)], levels = c("Air", "O2"))

dds1 <- DESeqDataSetFromMatrix(
  round(pb[rowSums(pb) > 0, ]),
  colData = data.frame(treat = grp, row.names = colnames(pb)),
  design = ~treat
)
dds1 <- DESeq(dds1, quiet = TRUE)
e <- as.data.frame(results(dds1))["Lgals3", ]
res <- rbind(res, data.frame(
  dataset = "Mouse AM (this study)",
  log2FC = e$log2FoldChange,
  lo = e$log2FoldChange - 1.96 * e$lfcSE,
  hi = e$log2FoldChange + 1.96 * e$lfcSE,
  p = e$pvalue, model = "DESeq2 pseudo-bulk", note = "n=4 pooled samples"
))
expect_value("Fig 7A log2FC (this study)", e$log2FoldChange, 0.68, tol = 0.02)

df1 <- data.frame(counts = pb["Lgals3", ], total = colSums(pb), treat = grp)
df1$cpm <- log2(df1$counts / df1$total * 1e6 + 1)
save_panel("Fig7A",
  boxplot_theme(ggplot(df1, aes(.data$treat, .data$cpm, color = .data$treat)),
                "Mouse AM (this study)", "Lgals3 log2 CPM"),
  data = data.frame(library = rownames(df1), df1[, c("treat", "cpm")]),
  width = 3, height = 4, slug = "mouse_AM"
)
rm(seu, am, pb)
invisible(gc())

# =============================================================================
# (2) Hurskainen mouse alveolar macrophages (GSE151974) - NB pseudo-bulk
# =============================================================================
f <- need(file.path(external_dir, "GSE151974",
                    "Seurat_object_GSE151974.RData"), "GSE151974")
env <- new.env(parent = emptyenv())
load(f, envir = env)
gse151974 <- env[[ls(env)[1]]]

am <- gse151974[, gse151974$CellType == "Alv Mf"]
cc <- LayerData(am, assay = "RNA", layer = "counts")
s <- as.character(am$Sample)
u <- sort(unique(s))
pb <- do.call(cbind, lapply(u, function(x) {
  Matrix::rowSums(cc[, s == x, drop = FALSE])
}))
colnames(pb) <- u
rownames(pb) <- rownames(cc)
md <- am@meta.data[match(u, am$Sample), c("Age", "Oxygen")]
md$Oxygen <- relevel(factor(md$Oxygen), ref = "Normoxia")
df2 <- data.frame(counts = pb["Lgals3", ], total = colSums(pb), md)

fit2 <- glmmTMB(counts ~ offset(log(total)) + Age + Oxygen,
                family = nbinom2, data = df2)
e <- summary(fit2)$coefficients$cond["OxygenHyperoxia", ]
res <- rbind(res, data.frame(
  dataset = "Mouse AM (Hurskainen)",
  log2FC = e[1] / log(2),
  lo = (e[1] - 1.96 * e[2]) / log(2),
  hi = (e[1] + 1.96 * e[2]) / log(2),
  p = e[4], model = "NB pseudo-bulk",
  note = "n=1/group: CI/p overconfident"
))
expect_value("Fig 7B log2FC (Hurskainen)", e[1] / log(2), 0.72,
             tol = 0.02, on_fail = "warn")

df2$cpm <- log2(df2$counts / df2$total * 1e6 + 1)
save_panel("Fig7B",
  boxplot_theme(ggplot(df2, aes(.data$Oxygen, .data$cpm, color = .data$Oxygen)),
                "Hurskainen mouse AM", "Lgals3 log2 CPM"),
  data = data.frame(sample = rownames(df2),
                    df2[, c("Age", "Oxygen", "cpm")]),
  width = 3, height = 4, slug = "hurskainen_AM"
)
rm(gse151974, am, cc, pb)
invisible(gc())

# =============================================================================
# (3) Human alveolar macrophages, LungMAP - NB mixed, donor random intercept
# =============================================================================
lm_dir <- file.path(external_dir, "lungmap", "LMEX0000004400")
sce <- readH5AD(need(file.path(lm_dir, "BPD-adata_combined.h5ad"),
                     "LMEX0000004400"), use_hdf5 = TRUE)
counts <- round(as(assay(sce, assayNames(sce)[1]), "CsparseMatrix"))
cd <- as.data.frame(colData(sce))
cd$barcode <- colnames(sce)

ann <- read.delim(need(file.path(lm_dir, "BPD_RNA_author-clusters.txt"),
                       "LMEX0000004400"),
                  check.names = FALSE, stringsAsFactors = FALSE)
colnames(ann) <- c("barcode", "library", "cond", "celltype")
cd$celltype <- ann$celltype[match(cd$barcode, ann$barcode)]

am_idx <- which(cd$celltype == "AM")
cc <- counts[, am_idx]
mm <- cd[am_idx, ]
s <- as.character(mm$sample)
u <- sort(unique(s))
pb <- do.call(cbind, lapply(u, function(x) {
  Matrix::rowSums(cc[, s == x, drop = FALSE])
}))
colnames(pb) <- u
rownames(pb) <- rownames(cc)
mi <- mm[match(u, mm$sample), ]

# samples contributing at least 10 alveolar macrophages
keep <- as.integer(table(s)[u]) >= 10
pb <- pb[, keep]
mi <- mi[keep, ]

grp <- factor(
  ifelse(mi$disease %in% c("aeBPD", "eBPD", "BPD"),
         "BPD_active", "control_healed"),
  levels = c("control_healed", "BPD_active")
)
df3 <- data.frame(counts = pb["LGALS3", ], total = colSums(pb),
                  treat = grp, mi)
fit3 <- glmmTMB(counts ~ offset(log(total)) + sex + treat + (1 | donor),
                family = nbinom2, data = df3)
e <- summary(fit3)$coefficients$cond["treatBPD_active", ]
res <- rbind(res, data.frame(
  dataset = "Human AM (LungMAP)",
  log2FC = e[1] / log(2),
  lo = (e[1] - 1.96 * e[2]) / log(2),
  hi = (e[1] + 1.96 * e[2]) / log(2),
  p = e[4], model = "NB mixed, donor RE", note = "author AM annotation"
))
expect_value("Fig 7C log2FC (LungMAP)", e[1] / log(2), 0.52,
             tol = 0.03, on_fail = "warn")
expect_value("Fig 7C P value (LungMAP)", e[4], 0.085,
             tol = 0.01, on_fail = "warn")

df3$cpm <- log2(df3$counts / df3$total * 1e6 + 1)
df3$disease <- factor(df3$disease,
                      levels = c("control", "aeBPD", "eBPD", "hBPD"))
save_panel("Fig7C",
  boxplot_theme(ggplot(df3, aes(.data$disease, .data$cpm,
                                color = .data$treat)),
                "Human AM (LungMAP)", "LGALS3 log2 CPM"),
  data = data.frame(sample = rownames(df3),
                    df3[, c("disease", "treat", "cpm")]),
  width = 4, height = 4, slug = "lungmap_AM"
)
rm(sce, counts, cc, pb)
invisible(gc())

# =============================================================================
# (4) Human blood, GSE32472 (microarray) - LMM on log2 intensity
# =============================================================================
env <- new.env(parent = emptyenv())
load(need(file.path(external_dir, "GSE32472", "gse32472.RData"), "GSE32472"),
     envir = env)
ex <- env$ex
gpl <- env$gpl
meta <- env$meta

ok <- which(gpl$`Gene symbol` == "LGALS3")
if (length(ok) > 1) ok <- ok[which.max(rowMeans(ex[ok, , drop = FALSE]))]
df4 <- data.frame(lgals3 = ex[ok, ], meta)
na_label <- "bronchopulmonary dysplasia (bpd) group: #N/A!"
df4 <- df4[df4$characteristics_ch1.7 != na_label, ]
df4$treat <- df4$bronchopulmonary.dysplasia..bpd..group.ch1
df4 <- df4[df4$treat %in% c("0 (no BPD)", "3 (severe BPD)"), ]
df4$treat <- factor(df4$treat, levels = c("0 (no BPD)", "3 (severe BPD)"))

fit4 <- lmer(lgals3 ~ gender.ch1 + timepoint + treat + (1 | donor), data = df4)
ct <- summary(fit4)$coefficients
e <- ct[grep("severe", rownames(ct))[1], ]
res <- rbind(res, data.frame(
  dataset = "Human blood (GSE32472)",
  log2FC = e[1], lo = e[1] - 1.96 * e[2], hi = e[1] + 1.96 * e[2],
  p = e[grep("Pr", colnames(ct))],
  model = "LMM (microarray log2 intensity)",
  note = "microarray scale (compressed)"
))
expect_value("Fig 7D1 log2FC (GSE32472)", e[1], 0.073,
             tol = 0.01, on_fail = "warn")

df4$group <- factor(ifelse(grepl("^3", df4$treat), "severe BPD", "no BPD"),
                    levels = c("no BPD", "severe BPD"))
save_panel("Fig7D1",
  boxplot_theme(ggplot(df4, aes(.data$group, .data$lgals3,
                                color = .data$group)),
                "Human blood (GSE32472)", "LGALS3 (log2 intensity)"),
  data = df4[, c("group", "timepoint", "lgals3")],
  width = 3, height = 4, slug = "gse32472_blood"
)

# =============================================================================
# (5) Human blood, GSE220135 (bulk RNA-seq) - NB with library-size offset
# =============================================================================
# Counts are fetched from GEO. Sample metadata is read from a frozen file
# rather than parsed out of a directory listing, which is what the original
# script did and which breaks whenever the download folder gains a file.
urld <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/",
               "?format=file&type=rnaseq_counts")
tbl <- as.matrix(fread(
  paste(urld, "acc=GSE220135",
        "file=GSE220135_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep = "&"),
  header = TRUE, colClasses = "integer"
), rownames = 1)

tr <- read.csv(need(file.path(data_dir, "frozen",
                              "GSE220135_sample_metadata.csv"), "GSE220135"),
               stringsAsFactors = FALSE)
tr <- tr[match(colnames(tbl), tr$id), ]
stopifnot(!anyNA(tr$sample))

df5 <- data.frame(lgals3 = tbl["3958", ], total = colSums(tbl), tr)
df5 <- df5[df5$timepoint != "Cord", ]
df5$treat <- relevel(factor(df5$treat), ref = "nonBPD")

fit5 <- glmmTMB(lgals3 ~ offset(log(total)) + timepoint + treat,
                family = nbinom2, data = df5)
e <- summary(fit5)$coefficients$cond["treatBPD", ]
res <- rbind(res, data.frame(
  dataset = "Human blood (GSE220135)",
  log2FC = e[1] / log(2),
  lo = (e[1] - 1.96 * e[2]) / log(2),
  hi = (e[1] + 1.96 * e[2]) / log(2),
  p = e[4], model = "NB (offset)", note = ""
))
expect_value("Fig 7D2 log2FC (GSE220135)", e[1] / log(2), 0.33,
             tol = 0.02, on_fail = "warn")
expect_value("Fig 7D2 P value (GSE220135)", e[4], 0.095,
             tol = 0.01, on_fail = "warn")

df5$cpm <- log2(df5$lgals3 / df5$total * 1e6 + 1)
save_panel("Fig7D2",
  boxplot_theme(ggplot(df5, aes(.data$treat, .data$cpm, color = .data$treat)),
                "Human blood (GSE220135)", "LGALS3 log2 CPM"),
  data = df5[, c("sample", "timepoint", "treat", "cpm")],
  width = 3, height = 4, slug = "gse220135_blood"
)

# =============================================================================
# Table S6 and Fig 7E
# =============================================================================
expect_value("datasets in the forest plot", nrow(res), 5)
expect_value("all effects in the same direction", all(res$log2FC > 0), TRUE)
save_table("TableS6_crossdataset_lgals3", res, path = forest_csv)

res$sig <- ifelse(res$p < 0.05, "p < 0.05", "trend")
res$dataset <- factor(res$dataset, levels = rev(res$dataset))
res$label <- sprintf("%.2f (p=%s)", res$log2FC,
                     formatC(res$p, format = "g", digits = 2))

p7e <- ggplot(res, aes(.data$log2FC, .data$dataset, color = .data$sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_errorbarh(aes(xmin = .data$lo, xmax = .data$hi),
                 height = 0.25, linewidth = 0.7) +
  geom_point(size = 3) +
  geom_text(aes(label = .data$label), vjust = -1.1, size = 2.7,
            color = "grey20") +
  scale_color_manual(values = c("p < 0.05" = "#d6604d", "trend" = "grey50"),
                     name = NULL) +
  labs(x = "log2 fold change (BPD vs control)", y = NULL) +
  theme_classic() +
  theme(legend.position = "top")
save_panel("Fig7E", p7e,
  data = res[, c("dataset", "log2FC", "lo", "hi", "p", "model", "note")],
  width = 6, height = 4, slug = "forest"
)

print(res[, c("dataset", "log2FC", "lo", "hi", "p")], digits = 3)
message("09_crossdataset.R complete")
