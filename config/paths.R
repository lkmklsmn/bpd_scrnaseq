# Path configuration. Sourced by every script.
#
# Reproduction starts from a frozen, checksummed Seurat object (see
# data/README.md), NOT from raw counts: SCTransform, UMAP and Louvain
# clustering are not stable across package versions, and every downstream
# script relies on the cluster numbering fixed in that object.
#
# Set BPD_DATA_DIR to wherever you downloaded the data. Everything else is
# derived from it, so no absolute paths appear in any script.

data_dir <- Sys.getenv("BPD_DATA_DIR", unset = file.path(getwd(), "data"))
if (!dir.exists(data_dir)) {
  stop("Data directory not found: ", data_dir, "\n",
       "Set BPD_DATA_DIR, or run scripts/00_download_external.R first.")
}

raw_dir      <- file.path(data_dir, "raw")
external_dir <- file.path(data_dir, "external")

results_dir  <- file.path(getwd(), "results")
fig_dir      <- file.path(results_dir, "figures")
tab_dir      <- file.path(results_dir, "tables")
int_dir      <- file.path(results_dir, "intermediate")
for (d in c(fig_dir, tab_dir, int_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ---- frozen analysis entry point -------------------------------------------
# 32,780 CD45+ cells; SCT + RNA assays; seurat_clusters at resolution 0.2
# (14 clusters); S.Score / G2M.Score; sample, sex and condition metadata.
seurat_obj_file <- file.path(raw_dir, "Seurat_object_scrublet.RData")

# ---- intermediates produced by the pipeline --------------------------------
# Supplemental tables are numbered in order of first citation in the
# manuscript, which is why abundance (cited in the Figure 2 paragraph) is S2
# and the cross-dataset effect sizes (cited at Figure 7) are S6.
markers_csv   <- file.path(tab_dir, "TableS1_celltype_markers.csv")
abundance_csv <- file.path(tab_dir, "TableS2_abundance.csv")
deseq_xlsx    <- file.path(tab_dir, "TableS3_de_by_celltype.xlsx")
gsea_xlsx     <- file.path(tab_dir, "TableS4_gsea_hallmark.xlsx")
drug_csv      <- file.path(tab_dir, "TableS5_drug_targets.csv")
forest_csv    <- file.path(tab_dir, "TableS6_crossdataset_lgals3.csv")
deseq_rds    <- file.path(int_dir, "de_by_celltype.rds")

# ---- cluster -> cell type annotation, defined in exactly one place ---------
# Fixed by the frozen object above; do not renumber.
cluster_annotations <- c(
  "0"  = "B cells",                  "1"  = "T cells",
  "2"  = "Monocytes",                "3"  = "Alveolar macrophages",
  "4"  = "Neutrophils",              "5"  = "Dendritic cells",
  "6"  = "NK cells",                 "7"  = "B cells",
  "8"  = "Cycling cells",            "9"  = "B cells",
  "10" = "Interstitial macrophages", "11" = "B cells",
  "12" = "B cells",                  "13" = "B cells"
)

# Canonical marker shown per cell type in Figure 1D. Each is drawn from the
# top 100 genes of its own cell type in Table S1.
fig1d_markers <- c(
  "B cells"                  = "Cd79a",
  "T cells"                  = "Cd3e",
  "Interstitial macrophages" = "C1qb",
  "Neutrophils"              = "S100a9",
  "Alveolar macrophages"     = "Atp6v0d2",
  "Monocytes"                = "Ace",
  "NK cells"                 = "Nkg7",
  "Dendritic cells"          = "Batf3",
  "Cycling cells"            = "Mki67"
)

# Thresholds used throughout, stated once.
padj_de      <- 0.25   # pseudo-bulk DE (4 pooled libraries; see Methods)
padj_marker  <- 0.05   # cell-type markers
min_pct      <- 0.1    # marker detection floor
n_heatmap    <- 100    # genes per cell type in Figure 1C
seed         <- 1234
set.seed(seed)
