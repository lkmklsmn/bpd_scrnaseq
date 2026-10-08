# bpd_scrnaseq

Reproducibility code for *"Single-cell RNA sequencing of CD45+ lung cells reveals the immune landscape of hyperoxia-induced bronchopulmonary dysplasia."*

CD45+-enriched single-cell RNA-seq of a neonatal mouse hyperoxia model of BPD (32,780 cells, 4 pooled libraries, 2 per exposure), with cross-dataset validation in mouse and human BPD, nominating galectin-3 (*Lgals3* / LGALS3) as a candidate therapeutic target.

**Data:** <https://doi.org/10.5281/zenodo.23249196>

Every panel in the paper is produced by exactly one script, and every script asserts the numbers quoted in the manuscript as it runs. A re-run that disagrees with the paper fails loudly rather than quietly producing different figures.

## Quick start

```bash
# 1. Get the data (~8.2 GB on disk). Either run the download script:
Rscript scripts/00_download_external.R

#    ...or download from Zenodo by hand and point the pipeline at it:
export BPD_DATA_DIR=/path/to/bpd_data     # Windows: setx BPD_DATA_DIR C:\bpd_data
Rscript scripts/00_download_external.R --verify

# 2. Run the pipeline, from the repository root
for s in 02 03 04 05 06 07 08 09 10 11 12; do Rscript scripts/${s}_*.R; done

# 3. Check the outputs against the committed reference numbers
Rscript scripts/99_verify.R
```

Scripts must be run **from the repository root**, in order, and `02` before anything else: later scripts consume intermediates written by earlier ones. No script contains an absolute path — everything derives from `BPD_DATA_DIR` (default: `./data`) and `config/paths.R`.

## Where the pipeline starts, and why

Reproduction starts from a **frozen, checksummed Seurat object** (`data/raw/Seurat_object_scrublet.RData`), not from raw counts.

`01_preprocess.R` documents how that object was built (decontX → QC → cell-cycle scoring → SCTransform → PCA/Louvain at resolution 0.2 → UMAP) but is **not on the critical path**. SCTransform, UMAP and Louvain clustering are not stable across package versions: a re-run can renumber the 14 clusters, which would silently invalidate the cluster → cell type map in `config/paths.R` and every annotation downstream of it. The script therefore refuses to run unless `BPD_ALLOW_PREPROCESS=1` is set.

`aux_qc_scrublet.py` documents the upstream scanpy/scrublet step for the same reason.

## Pipeline

| # | Script | Panels | Tables |
|---|--------|--------|--------|
| 00 | `00_download_external.R` | — | downloads and checksums all inputs |
| 01 | `01_preprocess.R` | *provenance only, gated* | — |
| 02 | `02_annotate.R` | Fig 1A–1D, Fig S2 | **S1** markers |
| 03 | `03_abundance.R` | Fig 2A, 2B | **S2** abundance |
| 04 | `04_pseudobulk_de.R` | Fig 3A | **S3** DE by cell type |
| 05 | `05_gsea.R` | Fig 3B, 3C, 3D | **S4** Hallmark GSEA |
| 06 | `06_macrophages.R` | Fig 4A–4G | — |
| 07 | `07_tcells.R` | Fig 5A–5D, 5F, 5G | — |
| 08 | `08_drug_targets.R` | Fig 6A–6D | **S5** druggability |
| 09 | `09_crossdataset.R` | Fig 7A–7E | **S6** effect sizes |
| 10 | `10_external.R` | Fig 2C, Fig 5E, Fig S1 | — |
| 11 | `11_sex_concordance.R` | Fig S3A–S3C | — |
| 12 | `12_supplement.R` | — | assembles the S1–S6 workbook |
| 99 | `99_verify.R` | — | regression-checks a re-run |

Dependencies between scripts: `04` writes `de_by_celltype.rds`, consumed by `05`, `06` and `08`; `04` also writes Table S3, read by `11`; `12` requires all six tables. Everything else reads only the frozen object and the external data.

Auxiliary, outside the numbered pipeline:
- `aux_opentargets_requery.R` — queries the *current* Open Targets schema. It does **not** reproduce Figure 6 (see caveats below); it exists as a path to re-deriving the annotation in future work.
- `aux_qc_scrublet.py` — the scanpy/scrublet doublet-detection step upstream of `01`.

## Outputs

Each panel writes two files to `results/figures/`:

```
Fig1C_marker_heatmap.pdf        the panel
Fig1C_marker_heatmap_data.csv   the numbers plotted in it
```

The CSV exists because re-running never produces byte-identical PDFs — fonts, embedded timestamps and graphics-device versions all differ — so "reproducible" cannot be checked on the files themselves. `99_verify.R` diffs the CSVs and Tables S1–S6 against committed reference copies in `reference/` at a tolerance of 1e-6, reports any panel whose numbers have moved, and exits non-zero so it can be wired into CI.

```bash
Rscript scripts/99_verify.R --snapshot   # populate reference/ from a run you trust
Rscript scripts/99_verify.R              # compare a later run against it
```

## Repository layout

```
config/paths.R      data roots, table filenames, cluster -> cell type map,
                    Figure 1D markers, every statistical threshold
R/figure_io.R       save_panel() / save_table() / load_annotated() /
                    expect_value() -- the assertion helper
R/plots.R           shared DE, abundance and enrichment helpers
scripts/            the numbered pipeline (above)
data/frozen/        version-controlled inputs that CANNOT be regenerated
data/raw/           }  downloaded, gitignored
data/external/      }  see data/README.md
results/            figures, tables, intermediates (gitignored)
reference/          committed reference outputs for 99_verify.R
environment/        conda environment + renv lockfile
```

`config/paths.R` is the single place to change thresholds. The values used in the paper: pseudo-bulk DE at adjusted *P* < 0.25, markers at Bonferroni-adjusted *P* < 0.05 with `min.pct = 0.1`, top 100 markers per cell type in Figure 1C, global seed 1234.

## Caveats — read before interpreting a re-run

These are the places where the pipeline is not a pure function of its inputs. Each is handled explicitly in the relevant script; none is hidden.

**The discovery cohort is 2 vs 2.** Four pooled libraries, two per exposure. Pseudo-bulk aggregation per library is replicate-aware but underpowered, which is why the manuscript uses an adjusted-*P* threshold of 0.25 and leans on external replication rather than on these p-values alone.

**T-cell subtypes are derived at run time** (`07`). Unlike the macrophage UMAPs in `06`, which are descriptive only, the central T-cell result (γδ cells 7.9% → 16.5% of T cells) depends on which cells land in which subcluster. The original script attached labels to hardcoded cluster numbers, which mislabels every population if a package version renumbers them. This port instead assigns each label by the marker the cluster actually expresses most highly, then asserts the resulting subtype sizes. A renumbering is absorbed automatically; a genuinely different partition stops the script.

**Figure 2C is a reconstruction** (`10`). It was produced interactively and the code was not saved. It is rebuilt from the description in the Results, including an explicit mapping of Hurskainen's immune clusters onto our nine populations. The reported ρ = 0.64 is checked as a **warning**, not an error, because that grouping is a reconstruction of an undocumented choice.

**Figure S1 cannot be reproduced exactly** (`10`). It calls `FindAllMarkers` with `max.cells.per.ident = 100`, which subsamples at random, and the original script set no seed. The seed from `config/paths.R` is set here so the repository is at least self-consistent run to run.

**Figure 6 depends on a frozen annotation** (`08`). `data/frozen/drug_target_annotation_2026-06-30.csv` came from an Open Targets GraphQL **v4** query on 30 June 2026. The responses were not retained and the schema has since changed, so the query cannot be replayed — hence the file is version-controlled rather than regenerated. The funnel it reproduces: 262 genes upregulated in alveolar macrophages → 257 human orthologs → 92 carrying a druggability annotation → 11 candidate targets. `aux_opentargets_requery.R` against the current schema returns a different list (202 rows, different tier vocabulary) and so does not reproduce the published figure.

**Figure 7 uses the published model, deliberately** (`09`). For GSE220135 the manuscript's model is `lgals3 ~ offset(log(total)) + timepoint + treat`. A later revision adding a sex covariate and a per-infant random intercept moves that estimate from +0.33 (*P* = 0.095) to +0.27 (*P* = 0.23). The published model is the one reproduced here; the difference is noted in the script header.

**Hallmark gene sets drift between MSigDB releases** (`05`). The NES values quoted in the Results are checked, but as warnings: a mismatch means the installed `msigdbr` differs from the one used for the paper. Restore the pinned version with `renv::restore()` before interpreting Fig 3B/3C or Table S4.

**Sex is balanced but confounded with library** (`11`). One male and one female library per exposure, so sex cannot confound the exposure contrast, but a sex × exposure interaction cannot be tested. Figure S3 asks only whether the main findings depend on sex; it does not test for one.

## Data

See [`data/README.md`](data/README.md) for accessions, the directory layout and provenance of each file.

In brief: single-cell data generated in this study are at GEO **GSE346853**; reanalyzed datasets are GEO **GSE151974** (Hurskainen et al., mouse hyperoxia), **GSE32472** and **GSE220135** (human blood), and LungMAP **LMEX0000004400** (human BPD). The exact processed objects the pipeline consumes are deposited on Zenodo at the DOI above, because several were processed before deposition and cannot be regenerated bit-for-bit from the raw accessions.

Only Figures 2C, 5E, 7 and S1 need the external datasets. Figures 1, 3, 4, 5A–D, 6, S2 and S3 run from the frozen object alone — a 1.6 GB download rather than 8.2 GB.

Two inputs are version-controlled in `data/frozen/` rather than downloaded, because nothing in this pipeline can regenerate them: the Open Targets annotation described above, and the GSE220135 sample metadata.

## Environment

- **R:** `renv::restore()` from `environment/renv.lock`. *(Currently `renv.lock.PLACEHOLDER` — generate the real lockfile with `renv::snapshot()`.)* Key packages: Seurat, DESeq2, glmmTMB, lme4/lmerTest, fgsea, msigdbr, speckle, celda, zellkonverter, babelgene, readxl/writexl, pheatmap, ggplot2.
- **Python** (provenance scripts only): `conda env create -f environment/environment.yml`.
- `scripts/00_download_external.R` needs the `digest` package for checksum verification; without it, size checks still run.

## Citation

Please cite the paper, the data DOI above, and — if you use the code — the archived release DOI.

## License

Code is released under the MIT License (see [`LICENSE`](LICENSE)).
