# Data

Raw and external data are **not** committed (see `.gitignore`). Recreate the layout below.

## This study (`data/raw/`)
- `filtered_feature_bc_matrix.counts_full.h5ad` — scanpy/scrublet output (input to script 02)
- `Seurat_object_scrublet.RData` — produced by script 02 (input to 03–10)
- Raw counts deposited at GEO: **GSE######** *(add accession on deposit)*

## External datasets (`data/external/`)
| Dataset | Accession | Files | Used in |
|---|---|---|---|
| Hurskainen et al. mouse hyperoxia | **GSE151974** | `GSE151974/Seurat_object_GSE151974.RData` | Fig 2C, Fig 7B, Fig S1/S9 |
| LungMAP human BPD | **LMEX0000004400** | `lungmap/LMEX0000004400/BPD-adata_combined.h5ad`, `BPD_RNA_author-clusters.txt` | Fig 7C, Fig S9 |
| Human blood microarray | **GSE32472** | `GSE32472/gse32472.RData` (ex, gpl, meta) | Fig 7D |
| Human blood bulk RNA-seq | **GSE220135** | `GSE220135/` sample filenames (id_sample_timepoint_treat); counts fetched from GEO at runtime | Fig 7D |

## Reference resources
- MSigDB Hallmark (mouse) via `msigdbr` — **pin the msigdbr version** (Hallmark content changes between releases).
- Open Targets Platform GraphQL API **v4**, accessed **30 June 2026** — cache the query response under `data/external/opentargets/` for exact reproducibility (API is version-dependent).

## Excluded
GSE275938 and GSE225881 were evaluated during development but are **not** part of the final analysis.
