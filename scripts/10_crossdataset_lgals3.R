# 10 - Galectin-3 conservation across mouse & human BPD datasets (Figure 7).
# Source: lgals3_forest_figure.R
# Each dataset is loaded and modeled INDEPENDENTLY (no shared helper functions) so
# per-dataset models/plots can be fine-tuned. Effect = log2 fold change (disease vs control).
# Inputs (see data/README for accessions):
#   this study : data/raw/Seurat_object_scrublet.RData
#   GSE151974  : data/external/GSE151974/Seurat_object_GSE151974.RData
#   LungMAP    : data/external/lungmap/LMEX0000004400/{BPD-adata_combined.h5ad, BPD_RNA_author-clusters.txt}
#   GSE32472   : data/external/GSE32472/gse32472.RData
#   GSE220135  : downloaded on the fly from GEO + data/external/GSE220135/ filenames

source(file.path("config", "paths.R"))
suppressMessages({
  library(data.table); library(ggplot2); library(Matrix); library(glmmTMB)
  library(lme4); library(lmerTest); library(zellkonverter); library(SingleCellExperiment)
  library(Seurat); library(DESeq2)
})
res <- data.frame()

# (1) Mouse alveolar macrophages (this study) — cluster 3 = AM, DESeq2 pseudo-bulk
load(seurat_obj_file)   # -> seu
am <- seu[, as.character(seu$seurat_clusters) == "3"]
s  <- as.character(am$orig.ident); u <- sort(unique(s))
pb <- do.call(cbind, lapply(u, function(x) Matrix::rowSums(am@assays$SCT@counts[, s == x, drop = FALSE])))
colnames(pb) <- u; rownames(pb) <- rownames(am@assays$SCT@counts)
grp <- factor(am$condition[match(u, am$orig.ident)], levels = c("Air", "O2"))
df1 <- data.frame(counts = pb["Lgals3", ], total = colSums(pb), treat = grp)
keepg <- rowSums(pb) > 0
dds1 <- DESeqDataSetFromMatrix(round(pb[keepg, ]),
          colData = data.frame(treat = grp, row.names = colnames(pb)), design = ~ treat)
dds1 <- DESeq(dds1); e <- as.data.frame(results(dds1))["Lgals3", ]
res <- rbind(res, data.frame(dataset = "Mouse AM (this study)",
  log2FC = e$log2FoldChange, lo = e$log2FoldChange - 1.96*e$lfcSE,
  hi = e$log2FoldChange + 1.96*e$lfcSE, p = e$pvalue,
  model = "DESeq2 pseudo-bulk", note = "n=4 pooled samples"))
df1$cpm <- log2(df1$counts / df1$total * 1e6 + 1)
p1 <- ggplot(df1, aes(treat, cpm, color = treat)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.12, size = 1.8) +
  scale_color_manual(values = c("grey40", "#d6604d")) +
  theme_classic() + theme(legend.position = "none") +
  labs(title = "Mouse AM (this study)", x = NULL, y = "Lgals3 log2 CPM")
rm(seu, am, pb); gc()

# (2) Hurskainen mouse AM (GSE151974) — CellType "Alv Mf", NB pseudo-bulk
load(file.path(external_dir, "GSE151974", "Seurat_object_GSE151974.RData"))  # -> gse151974
am <- gse151974[, gse151974$CellType == "Alv Mf"]
cc <- LayerData(am, assay = "RNA", layer = "counts")
s <- as.character(am$Sample); u <- sort(unique(s))
pb <- do.call(cbind, lapply(u, function(x) Matrix::rowSums(cc[, s == x, drop = FALSE])))
colnames(pb) <- u; rownames(pb) <- rownames(cc)
md <- am@meta.data[match(u, am$Sample), c("Age", "Oxygen")]
md$Oxygen <- relevel(factor(md$Oxygen), ref = "Normoxia")
df2 <- data.frame(counts = pb["Lgals3", ], total = colSums(pb), md)
fit2 <- glmmTMB(counts ~ offset(log(total)) + Age + Oxygen, family = nbinom2, data = df2)
e <- summary(fit2)$coefficients$cond["OxygenHyperoxia", ]
res <- rbind(res, data.frame(dataset = "Mouse AM (Hurskainen)",
  log2FC = e[1]/log(2), lo = (e[1]-1.96*e[2])/log(2), hi = (e[1]+1.96*e[2])/log(2),
  p = e[4], model = "NB pseudo-bulk", note = "n=1/group: CI/p overconfident"))
df2$cpm <- log2(df2$counts / df2$total * 1e6 + 1)
p2 <- ggplot(df2, aes(Oxygen, cpm, color = Oxygen)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.12, size = 1.8) +
  scale_color_manual(values = c("grey40", "#d6604d")) +
  theme_classic() + theme(legend.position = "none") +
  labs(title = "Hurskainen mouse AM", x = NULL, y = "Lgals3 log2 CPM")
rm(gse151974, am, cc, pb); gc()

# (3) Human AM — original LungMAP + author annotation, NB mixed (donor RE + sex)
sce <- readH5AD(file.path(external_dir, "lungmap", "LMEX0000004400", "BPD-adata_combined.h5ad"), use_hdf5 = TRUE)
counts <- round(as(assay(sce, assayNames(sce)[1]), "CsparseMatrix"))
cd <- as.data.frame(colData(sce)); cd$barcode <- colnames(sce)
ann <- read.delim(file.path(external_dir, "lungmap", "LMEX0000004400", "BPD_RNA_author-clusters.txt"),
                  check.names = FALSE, stringsAsFactors = FALSE)
colnames(ann) <- c("barcode", "library", "cond", "celltype")
cd$celltype <- ann$celltype[match(cd$barcode, ann$barcode)]
am <- which(cd$celltype == "AM"); cc <- counts[, am]; mm <- cd[am, ]
s <- as.character(mm$sample); u <- sort(unique(s))
pb <- do.call(cbind, lapply(u, function(x) Matrix::rowSums(cc[, s == x, drop = FALSE])))
colnames(pb) <- u; rownames(pb) <- rownames(cc)
mi <- mm[match(u, mm$sample), ]; ncells <- as.integer(table(s)[u]); keep <- ncells >= 10
pb <- pb[, keep]; mi <- mi[keep, ]
grp <- factor(ifelse(mi$disease %in% c("aeBPD", "eBPD", "BPD"), "BPD_active", "control_healed"),
              levels = c("control_healed", "BPD_active"))
df3 <- data.frame(counts = pb["LGALS3", ], total = colSums(pb), treat = grp, mi)
fit3 <- glmmTMB(counts ~ offset(log(total)) + sex + treat + (1 | donor), family = nbinom2, data = df3)
e <- summary(fit3)$coefficients$cond["treatBPD_active", ]
res <- rbind(res, data.frame(dataset = "Human AM (LungMAP)",
  log2FC = e[1]/log(2), lo = (e[1]-1.96*e[2])/log(2), hi = (e[1]+1.96*e[2])/log(2),
  p = e[4], model = "NB mixed, donor RE", note = "author AM annotation"))
df3$cpm <- log2(df3$counts / df3$total * 1e6 + 1)
df3$disease <- factor(df3$disease, levels = c("control", "aeBPD", "eBPD", "hBPD"))
p3 <- ggplot(df3, aes(disease, cpm, color = treat)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.12, size = 1.8) +
  scale_color_manual(values = c("grey40", "#d6604d")) +
  theme_classic() + theme(legend.position = "none") +
  labs(title = "Human AM (LungMAP)", x = NULL, y = "LGALS3 log2 CPM")
rm(sce, counts, cc, pb); gc()

# (4) Human blood, severe BPD (GSE32472, microarray) — LMM on log2 intensity
load(file.path(external_dir, "GSE32472", "gse32472.RData"))   # ex, gpl, meta
ok <- which(gpl$`Gene symbol` == "LGALS3")
if (length(ok) > 1) ok <- ok[which.max(rowMeans(ex[ok, , drop = FALSE]))]
df4 <- data.frame(lgals3 = ex[ok, ], meta)
df4 <- df4[df4$characteristics_ch1.7 != "bronchopulmonary dysplasia (bpd) group: #N/A!", ]
df4$treat <- df4$bronchopulmonary.dysplasia..bpd..group.ch1
df4 <- df4[df4$treat %in% c("0 (no BPD)", "3 (severe BPD)"), ]
df4$treat <- factor(df4$treat, levels = c("0 (no BPD)", "3 (severe BPD)"))
fit4 <- lmer(lgals3 ~ gender.ch1 + timepoint + treat + (1 | donor), data = df4)
ct <- summary(fit4)$coefficients; e <- ct[grep("severe", rownames(ct))[1], ]
res <- rbind(res, data.frame(dataset = "Human blood (GSE32472)",
  log2FC = e[1], lo = e[1]-1.96*e[2], hi = e[1]+1.96*e[2], p = e[grep("Pr", colnames(ct))],
  model = "LMM (microarray log2 intensity)", note = "microarray scale (compressed)"))
df4$group <- factor(ifelse(grepl("^3", df4$treat), "severe BPD", "no BPD"), levels = c("no BPD", "severe BPD"))
p4 <- ggplot(df4, aes(group, lgals3, color = group)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.12, size = 1.8) +
  scale_color_manual(values = c("grey40", "#d6604d")) +
  theme_classic() + theme(legend.position = "none") +
  labs(title = "Human blood (GSE32472)", x = NULL, y = "LGALS3 (log2 intensity)")

# (5) Human blood (GSE220135, bulk RNA-seq) — NB with library-size offset
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
tbl <- as.matrix(fread(paste(urld, "acc=GSE220135",
        "file=GSE220135_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep = "&"),
        header = TRUE, colClasses = "integer"), rownames = 1)
files <- list.files(file.path(external_dir, "GSE220135"))[-1]
tr <- do.call(rbind, lapply(files, function(x) strsplit(x, "_", fixed = TRUE)[[1]]))[, 1:4]
colnames(tr) <- c("id", "sample", "timepoint", "treat"); tr <- data.frame(tr)
tr <- tr[match(colnames(tbl), tr$id), ]
df5 <- data.frame(lgals3 = tbl["3958", ], total = colSums(tbl), tr)
df5 <- df5[df5$timepoint != "Cord", ]
df5$treat <- relevel(factor(df5$treat), ref = "nonBPD")
fit5 <- glmmTMB(lgals3 ~ offset(log(total)) + timepoint + treat, family = nbinom2, data = df5)
e <- summary(fit5)$coefficients$cond["treatBPD", ]
res <- rbind(res, data.frame(dataset = "Human blood (GSE220135)",
  log2FC = e[1]/log(2), lo = (e[1]-1.96*e[2])/log(2), hi = (e[1]+1.96*e[2])/log(2),
  p = e[4], model = "NB (offset)", note = ""))
df5$cpm <- log2(df5$lgals3 / df5$total * 1e6 + 1)
p5 <- ggplot(df5, aes(treat, cpm, color = treat)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.12, size = 1.8) +
  scale_color_manual(values = c("grey40", "#d6604d")) +
  theme_classic() + theme(legend.position = "none") +
  labs(title = "Human blood (GSE220135)", x = NULL, y = "LGALS3 log2 CPM")

# ---- Forest plot (Fig 7E) assembled from per-dataset coefficients -----------
write.csv(res, file.path(tab_dir, "lgals3_forest_res.csv"), row.names = FALSE)  # Table S4
res$sig <- ifelse(res$p < 0.05, "p < 0.05", "trend")
res$dataset <- factor(res$dataset, levels = rev(res$dataset))
res$label <- sprintf("%.2f (p=%s)", res$log2FC, formatC(res$p, format = "g", digits = 2))
pf <- ggplot(res, aes(log2FC, dataset, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 0.7) +
  geom_point(size = 3) + geom_text(aes(label = label), vjust = -1.1, size = 2.7, color = "grey20") +
  scale_color_manual(values = c("p < 0.05" = "#d6604d", "trend" = "grey50"), name = NULL) +
  labs(x = "log2 fold change (BPD vs control)", y = NULL) + theme_classic() +
  theme(legend.position = "top")

f <- gridExtra::grid.arrange(p1, p2, p3, p4, p5, pf, nrow = 2)
ggsave(file.path(fig_dir, "Fig7_lgals3_compilation.pdf"), f, width = 15, height = 10)
print(res, digits = 3)
