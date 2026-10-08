# Data

This pipeline needs 8.2 GB of input data. None of it is committed except the two
frozen files in `data/frozen/` (see below).

**Archive:** <https://doi.org/10.5281/zenodo.23249196>

```bash
Rscript scripts/00_download_external.R           # download, then verify
Rscript scripts/00_download_external.R --verify  # verify what is already on disk
```

The download script checks every file on both size and SHA-256 and refuses to proceed
on a mismatch, so a truncated or corrupted download fails immediately rather than
surfacing later as a wrong number in a figure.

To keep the data outside the repository — recommended, since the repository lives in a
synced folder for most of us — set `BPD_DATA_DIR`:

```bash
export BPD_DATA_DIR=/path/to/bpd_data      # Windows: setx BPD_DATA_DIR C:\bpd_data
```

Every path in every script derives from that one variable via `config/paths.R`. No
script contains an absolute path.

## Why the data are deposited rather than fetched from source

Several of these objects were processed before deposition — filtered, normalised,
annotated — and the exact intermediate files cannot be regenerated bit-for-bit from the
raw accessions, because the original processing scripts were not retained by the
depositing groups and the package versions are unknown. Reproducing the figures
therefore requires the processed objects themselves, not the raw accessions they came
from. The original accessions are listed below for provenance and should be cited
alongside this work.

## Layout

After extraction, `$BPD_DATA_DIR` looks like this:

```
raw/
  Seurat_object_scrublet.RData                       1.6 GB
external/
  GSE151974/Seurat_object_GSE151974.RData            2.6 GB
  lungmap/LMEX0000004400/BPD-adata_combined.h5ad     4.5 GB
  lungmap/LMEX0000004400/BPD_RNA_author-clusters.txt  16 MB
  GSE32472/gse32472.RData                             80 MB
frozen/                                      (version-controlled, not downloaded)
  drug_target_annotation_2026-06-30.csv
  GSE220135_sample_metadata.csv
```

## `data/raw/` — this study

### `Seurat_object_scrublet.RData` — the frozen entry point
32,780 CD45+ cells × 17,739 genes. SCT and RNA assays, `seurat_clusters` at resolution
0.2 (14 clusters), cell-cycle scores, and `sample` / `sex` / `condition` metadata.

Reproduction starts here, **not** from raw counts. SCTransform, PCA, UMAP and Louvain
clustering are not stable across package versions: a re-run can renumber the 14
clusters, which would silently invalidate the cluster → cell type map in
`config/paths.R` and every annotation downstream of it. `scripts/01_preprocess.R`
documents how this object was built and is gated behind `BPD_ALLOW_PREPROCESS=1` so it
cannot overwrite it by accident.

SHA-256 `a7b2aeae7c9bf628d74fb2386bad845093739e9f6bc96571ce443bc8daeaba71`

Required by: **every script.**

### `filtered_feature_bc_matrix.counts_full.h5ad` — not in the archive
The scanpy/scrublet output upstream of `01_preprocess.R`. Not distributed here; the raw
counts are deposited at GEO under **GSE346853**. You need this only if you intend to
rebuild the frozen object, which is not part of reproducing the manuscript.

## `data/external/` — reanalyzed datasets

| Dataset | Accession | Files | Required by |
|---|---|---|---|
| Hurskainen et al., mouse hyperoxia scRNA-seq | **GSE151974** | `GSE151974/Seurat_object_GSE151974.RData` | Fig 2C, Fig 5E, Fig 7B, Fig S1 |
| LungMAP, human BPD lung scRNA-seq | **LMEX0000004400** | `lungmap/LMEX0000004400/BPD-adata_combined.h5ad`, `BPD_RNA_author-clusters.txt` | Fig 7C |
| Human blood, microarray | **GSE32472** | `GSE32472/gse32472.RData` (`ex`, `gpl`, `meta`) | Fig 7D1 |
| Human blood, bulk RNA-seq | **GSE220135** | none — counts fetched from GEO at run time | Fig 7D2 |

SHA-256, in the order listed above:

```
b14fe3af379814ce4ba53c262d8344f44209dfbecaeddf01ff104b23287f9696
e1c21358c293f69ade480a038af214c8c1a511ad44af700d5a868549d997ae68
90f08675a1941492f375d97036b9970f545409efe2afcc645e98994358e380df
8cf306852320a77cecb799242f8266cab3183a3d00ed1b913a8fd48e5d1bf194
```

GSE220135 counts are fetched from GEO when `09_crossdataset.R` runs, so that dataset
needs no local file — but its **sample metadata is frozen** (next section), because
deriving it from directory listing order is unsafe.

Only Figures 2C, 5E, 7 and S1 need any of these files. Figures 1, 3, 4, 5A–D, 6, S2 and
S3 run from the frozen object alone, so you can reproduce most of the paper with the
1.6 GB download and skip the rest.

## `data/frozen/` — version-controlled, not downloadable

These two files are committed to the repository because **nothing in this pipeline can
regenerate them.** They are inputs, not outputs.

### `drug_target_annotation_2026-06-30.csv`
Produced by a query against the Open Targets Platform GraphQL API **v4** on
**30 June 2026**. The API responses were not retained and the schema has since changed
— `Target.knownDrugs` and `KnownDrug.drug.maxPhaseForIndication` no longer exist — so
the query cannot be replayed. This file is the annotation exactly as deposited with the
manuscript (Supplemental Table S5), and it is what `08_drug_targets.R` reads.

`scripts/aux_opentargets_requery.R` queries the current schema and caches its responses
under `external/opentargets_<date>/`. It returns a different target list (202 rows, a
different tier vocabulary) and therefore **does not reproduce Figure 6.** It sits
outside the numbered pipeline and exists only as a path to re-deriving the annotation in
future work.

### `GSE220135_sample_metadata.csv`
Subject, timepoint and BPD status for the 128 GSE220135 libraries (`id`, `sample`,
`timepoint`, `treat`). Frozen because the working script derived these labels from the
order of a directory listing, which scrambled 40 of the 128 sex assignments. Reading
them from a file instead removes a class of error that produces plausible-looking but
wrong results.

## Reference resources

- **MSigDB Hallmark (mouse)** via `msigdbr`. Hallmark set membership changes between
  MSigDB releases, so the version is pinned in `environment/renv.lock`. `05_gsea.R`
  checks the NES values quoted in the Results and **warns** on a mismatch, which means
  the installed `msigdbr` differs from the one used for the paper.
- **Open Targets Platform** — frozen as described above, not queried at run time.

## Evaluated but not used

GSE275938 and GSE225881 were examined during development and are **not** part of the
final analysis. They are not required, not downloaded, and no figure depends on them.
