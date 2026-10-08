# 12 - Assemble the supplemental tables workbook.
#
# Inputs  Tables S1-S6 written by scripts 02-09
# Output  results/tables/Supplemental_Tables_S1_S6.xlsx
#
# Supplemental tables are numbered in ORDER OF FIRST CITATION in the
# manuscript, which is why differential abundance (cited in the Figure 2
# paragraph) is S2 and the cross-dataset effect sizes (cited at Figure 7) are
# S6. Do not renumber without renumbering the citations.
#
# The workbook that accompanied the submission was assembled by hand; this
# script replaces that step so the deposited file can be regenerated.

source(file.path("config", "paths.R"))
source(file.path("R", "figure_io.R"))
suppressMessages({
  library(readxl)
  library(writexl)
})

out_file <- file.path(tab_dir, "Supplemental_Tables_S1_S6.xlsx")

spec <- list(
  list(id = "S1", sheet = "S1_Markers", path = markers_csv, kind = "csv",
       title = "Cell type marker genes",
       script = "02_annotate.R",
       desc = paste(
         "All significantly enriched marker genes per annotated immune",
         "population (Seurat FindAllMarkers, Wilcoxon rank-sum test, positive",
         "markers only, detected in at least 10% of cells in either group,",
         "Bonferroni-adjusted P < 0.05), ranked within each population by the",
         "difference in the proportion of expressing cells (pct.1 - pct.2).",
         "The top 100 genes per population are shown in Figure 1C.")),
  list(id = "S2", sheet = "S2_Abundance", path = abundance_csv, kind = "csv",
       title = "Differential cell type abundance, hyperoxia vs room air",
       script = "03_abundance.R",
       desc = paste(
         "propeller differential abundance test for each immune population",
         "between hyperoxia and room air, treating the four pooled libraries",
         "as the replication unit. Supports Figure 2, A and B.")),
  list(id = "S3", sheet = "S3_DE", path = deseq_xlsx, kind = "xlsx",
       title = "Pseudo-bulk differential expression by cell type",
       script = "04_pseudobulk_de.R",
       desc = paste(
         "DESeq2 pseudo-bulk differential expression between hyperoxia and",
         "room air within each cell type, counts aggregated per library",
         "(4 libraries, 2 per condition). Genes are called differentially",
         "expressed at Benjamini-Hochberg adjusted P < 0.25.")),
  list(id = "S4", sheet = "S4_GSEA", path = gsea_xlsx, kind = "xlsx",
       title = "Hallmark gene set enrichment by cell type",
       script = "05_gsea.R",
       desc = paste(
         "fgsea results for MSigDB Hallmark gene sets per cell type, ranked",
         "on the pseudo-bulk log2 fold change (hyperoxia versus room air).")),
  list(id = "S5", sheet = "S5_DrugTargets", path = drug_csv, kind = "csv",
       title = "Druggability annotation of alveolar macrophage-induced genes",
       script = "08_drug_targets.R",
       desc = paste(
         "Human orthologs of the 262 genes upregulated in alveolar",
         "macrophages, annotated against the Open Targets Platform (query of",
         "30 June 2026). Orthologs with no Open Targets record are retained",
         "and marked, so the 262 -> 257 -> annotated funnel can be checked",
         "directly.")),
  list(id = "S6", sheet = "S6_CrossDataset", path = forest_csv, kind = "csv",
       title = "Cross-dataset galectin-3 (LGALS3) effect sizes",
       script = "09_crossdataset.R",
       desc = paste(
         "Harmonised galectin-3 effect size (log2 fold change, 95% CI, P) for",
         "each mouse and human dataset in Figure 7, with the model used for",
         "each."))
)

missing <- vapply(spec, function(s) !file.exists(s$path), logical(1))
if (any(missing)) {
  stop("Cannot assemble the workbook; these tables have not been generated:\n",
       paste(sprintf("  %s  (%s, from %s)",
                     vapply(spec[missing], `[[`, character(1), "id"),
                     basename(vapply(spec[missing], `[[`, character(1), "path")),
                     vapply(spec[missing], `[[`, character(1), "script")),
             collapse = "\n"),
       call. = FALSE)
}

sheets <- list()
contents <- data.frame()

for (s in spec) {
  if (identical(s$kind, "csv")) {
    sheets[[s$sheet]] <- read.csv(s$path, stringsAsFactors = FALSE,
                                  check.names = FALSE)
    n_rows <- nrow(sheets[[s$sheet]])
  } else {
    # multi-sheet inputs (one sheet per cell type) are stacked into one table
    tabs <- excel_sheets(s$path)
    parts <- lapply(tabs, function(tb) {
      d <- as.data.frame(read_excel(s$path, sheet = tb))
      cbind(cell_type = tb, d)
    })
    sheets[[s$sheet]] <- do.call(rbind, parts)
    n_rows <- nrow(sheets[[s$sheet]])
  }
  contents <- rbind(contents, data.frame(
    Table = s$id, Title = s$title, Description = s$desc,
    N_rows = n_rows, Source_script = s$script,
    stringsAsFactors = FALSE
  ))
  message(sprintf("  %-18s %7d rows  <- %s", s$sheet, n_rows,
                  basename(s$path)))
}

write_xlsx(c(list(Contents = contents), sheets), path = out_file)
message("  workbook -> ", basename(out_file), " (",
        length(sheets) + 1, " sheets)")

# The row counts reported in the manuscript and in data/README.md
expect_value("Table S1 rows", contents$N_rows[contents$Table == "S1"], 8388)
expect_value("Table S2 rows", contents$N_rows[contents$Table == "S2"], 9)
expect_value("Table S6 rows", contents$N_rows[contents$Table == "S6"], 5)

message("12_supplement.R complete")
