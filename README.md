# bpd_scrnaseq

Reproducibility code for *"Single-cell RNA sequencing of CD45+ lung cells reveals the immune landscape of hyperoxia-induced bronchopulmonary dysplasia."*

CD45+-enriched single-cell RNA-seq of a neonatal mouse hyperoxia model of BPD, with cross-dataset (mouse + human) validation, nominating galectin-3 as a candidate therapeutic target.

## Pipeline

Scripts run in order (`scripts/`). Python does raw processing/doublets; R does everything downstream. Paths are centralized in `config/paths.R` — edit the two roots there once.

| # | Script | Produces |
|---|--------|----------|
| 01 | `01_qc_preprocess.py` | scanpy + scrublet → counts `.h5ad` |
| 02 | `02_preprocess_cluster.R` | decontX, QC, SCTransform, clustering → `Seurat_object_scrublet.RData`; Fig S3 |
| 03 | `03_annotate_and_atlas.R` | annotation; **Fig 1**; Fig S2; marker tables |
| 04 | `04_differential_abundance.R` | propeller; **Fig 2A/B**; Fig S7 |
| 05 | `05_pseudobulk_de.R` | per-cell-type DESeq2; **Fig 3A**; Table S1 |
| 06 | `06_gsea_and_pca.R` | Hallmark GSEA, interferon, PCA; **Fig 3B/C/D**; Table S2 |
| 07 | `07_macrophages.R` | AM/IM programs; **Fig 4**; Fig S5, S8 |
| 08 | `08_tcells.R` | T-cell subclustering, γδ/IL-17; **Fig 5** |
| 09 | `09_drug_targets.R` | druggability analysis; **Fig 6**; Table S3 |
| 10 | `10_crossdataset_lgals3.R` | galectin-3 conservation forest; **Fig 7**; Fig S9; Table S4 |

Auxiliary scripts:
- `09a_opentargets_query.R` — babelgene → Open Targets (GraphQL v4) druggability query; writes the curated CSV consumed by `09`.
- `aux_hurskainen_frequency.R` — Hurskainen (GSE151974) cross-study validation: **Fig 2C** scatter, **Fig S1** correlation heatmap, **Fig S9** Lgals3 AM violins.

Shared helpers: `scripts/utils/functions.R` (discovery-data plotting/DE/GSEA) and `scripts/utils/load_annotate.R` (loads the object + cell-type labels). Script 10 is intentionally self-contained (no helpers) so each dataset's model can be tuned independently.

## Run

```r
# from repo root, in R
source("scripts/02_preprocess_cluster.R")
source("scripts/03_annotate_and_atlas.R")
# ... 04 -> 10
```
Data are not committed; see `data/README.md` for accessions and `scripts/00_download_data.sh` (to add) for retrieval.

## Environment
- R deps: restore with `renv::restore()` after you run `renv::snapshot()` (see `environment/` - `renv.lock` is a placeholder to be generated on your machine).
- Python deps: `conda env create -f environment/environment.yml`.
- Global seed and version pins (msigdbr Hallmark, Open Targets API date) are noted in `config/paths.R` and the relevant scripts.

## Notes / caveats
- Discovery cohort = 4 pooled libraries, 2 per condition - pseudobulk DE is replicate-aware but underpowered; treat extreme p-values with care.
- `09a_opentargets_query.R` is a **reconstruction** of the Open Targets query (API v4, accessed 30 June 2026); verify its output against the reported counts (257 orthologs mapped, 75 druggable) and cache responses under `data/external/opentargets/`, since the API is version-dependent.
- `aux_hurskainen_frequency.R` (Fig 2C) uses a cell-type name mapping between the two datasets — **verify the `map_ct` vector** before trusting the scatter.
- `00_download_data.sh` (GEO/LungMAP retrieval) and `environment/renv.lock` are still to be generated on your machine.

## Citation / archiving
Tag `v1.0-submission` and mint a Zenodo DOI from the GitHub release; cite that DOI in the manuscript's Code Availability statement.
