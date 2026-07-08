#!/usr/bin/env python
"""
01 - Raw processing and doublet flagging (Python / scanpy).

Upstream of the R pipeline. Produces the counts h5ad that script 02 reads with
decontX. This documents the scanpy + scrublet step described in the Methods;
fill in with your notebook if it differs.

Inputs : CellRanger filtered_feature_bc_matrix (10x, 5' GEX) per library
Output : data/raw/filtered_feature_bc_matrix.counts_full.h5ad  (-> script 02)

Methods reference: scanpy for processing; scrublet for doublet flagging;
ambient-RNA removal (decontX) is done downstream in R (script 02).
"""
import scanpy as sc
import scrublet as scr
import numpy as np

# --- load (edit path/loader to your CellRanger output) -----------------------
adata = sc.read_10x_mtx("data/raw/filtered_feature_bc_matrix",
                        var_names="gene_symbols", cache=True)

# --- doublet detection (per library) ----------------------------------------
scrub = scr.Scrublet(adata.X)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
adata.obs["scrublet_score"]   = doublet_scores
adata.obs["predicted_doublet"] = predicted_doublets

# Keep counts intact (decontX in R needs raw counts); export.
adata.write("data/raw/filtered_feature_bc_matrix.counts_full.h5ad")
