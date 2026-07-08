# Load the discovery Seurat object and attach cell-type annotations.
# Sourced by scripts 03-10 so annotation is defined in exactly one place.
# Provides `seu` with `seu$celltype`.

load(seurat_obj_file)  # -> seu

cluster_annotations <- c(
  "0"  = "B cells",              "1"  = "T cells",
  "2"  = "Monocytes",            "3"  = "Alveolar macrophages",
  "4"  = "Neutrophils",          "5"  = "Dendritic cells",
  "6"  = "NK cells",             "7"  = "B cells",
  "8"  = "Cycling cells",        "9"  = "B cells",
  "10" = "Interstitial macrophages",
  "11" = "B cells",             "12" = "B cells",  "13" = "B cells")

seu@meta.data$celltype <- unname(cluster_annotations[as.character(seu$seurat_clusters)])

# Back-compat: objects built before the gender->sex rename carry `gender`.
if (is.null(seu@meta.data$sex) && !is.null(seu@meta.data$gender))
  seu@meta.data$sex <- seu@meta.data$gender
