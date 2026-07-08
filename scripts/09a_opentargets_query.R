# 09a - Druggability annotation via babelgene + Open Targets Platform (GraphQL v4).
# Produces: results/tables/drug_target_candidates_curated.csv  (consumed by 09)
#
# Reconstruction of the query run during analysis (Open Targets GraphQL API v4,
# accessed 30 June 2026). The API is version-dependent: cache responses under
# data/external/opentargets/ and record the access date for exact reproducibility.
# VERIFY against the paper: 262 AM-upregulated genes -> 257 human orthologs mapped,
# 75 with an approved/clinical-stage drug.
#
# Output columns (as expected by 09_drug_targets.R):
#   mouse_genes | human_symbol | ensembl_id | celltypes | tier | example_drugs

source(file.path("config", "paths.R"))
suppressMessages({ library(readxl); library(babelgene); library(httr); library(jsonlite) })

cache_dir <- file.path(external_dir, "opentargets"); dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
OT <- "https://api.platform.opentargets.org/api/v4/graphql"

# ---- 1. AM-upregulated genes (padj < 0.25, log2FC > 0) ----------------------
am <- as.data.frame(read_excel(deseq_xlsx, sheet = "Alveolar macrophages"))
up_mouse <- am$gene[!is.na(am$padj) & am$padj < 0.25 & am$log2FoldChange > 0]

# ---- 2. mouse -> human orthologs --------------------------------------------
orth <- orthologs(genes = up_mouse, species = "mouse", human = FALSE)  # human_symbol, symbol(mouse)
map <- data.frame(mouse = orth$symbol, human = orth$human_symbol, stringsAsFactors = FALSE)
map <- map[!is.na(map$human) & map$human != "", ]

# ---- 3. GraphQL helpers (cached) --------------------------------------------
gql <- function(query, variables){
  res <- POST(OT, body = list(query = query, variables = variables), encode = "json")
  content(res, as = "parsed", type = "application/json")
}
resolve_ensembl <- function(symbol){
  q <- 'query($q:String!){ search(queryString:$q, entityNames:["target"]){ hits{ id name } } }'
  hit <- tryCatch(gql(q, list(q = symbol))$data$search$hits, error = function(e) NULL)
  if (length(hit)) hit[[1]]$id else NA_character_
}
target_info <- function(ensembl_id){
  f <- file.path(cache_dir, paste0(ensembl_id, ".json"))
  if (file.exists(f)) return(fromJSON(f, simplifyVector = FALSE))
  q <- 'query($id:String!){ target(ensemblId:$id){
          approvedSymbol
          tractability{ modality value label }
          knownDrugs{ rows{ drug{ name maxPhaseForIndication } } } } }'
  out <- gql(q, list(id = ensembl_id))
  write_json(out, f, auto_unbox = TRUE); out
}

# ---- 4. classify each ortholog ----------------------------------------------
classify <- function(symbol){
  eid <- resolve_ensembl(symbol)
  if (is.na(eid)) return(data.frame(human_symbol = symbol, ensembl_id = NA,
                                    tier = "Not druggable", example_drugs = ""))
  ti <- tryCatch(target_info(eid)$data$target, error = function(e) NULL)
  drugs <- character(0)
  if (!is.null(ti$knownDrugs$rows))
    drugs <- unique(vapply(ti$knownDrugs$rows, function(r) r$drug$name %||% "", character(1)))
  drugs <- drugs[drugs != ""]
  tract <- if (!is.null(ti$tractability))
    vapply(ti$tractability, function(t) isTRUE(t$value), logical(1)) else logical(0)
  tier <- if (length(drugs)) "A: existing drug" else
          if (any(tract))   "B: clinically tractable" else
                            "Druggable"                # tractability record exists but low
  data.frame(human_symbol = symbol, ensembl_id = eid, tier = tier,
             example_drugs = paste(head(drugs, 3), collapse = ";"), stringsAsFactors = FALSE)
}
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

ann <- do.call(rbind, lapply(unique(map$human), classify))

# ---- 5. assemble curated table ----------------------------------------------
ann$mouse_genes <- vapply(ann$human_symbol, function(h)
  paste(unique(map$mouse[map$human == h]), collapse = "/"), character(1))
ann$celltypes <- "Alveolar macrophages"
# drop non-specific housekeeping (ribosomal / mitochondrial) as in Methods
ann <- ann[!grepl("^RP[LS]|^MRP|^MT-", ann$human_symbol), ]

write.csv(ann[, c("mouse_genes","human_symbol","ensembl_id","celltypes","tier","example_drugs")],
          file.path(tab_dir, "drug_target_candidates_curated.csv"), row.names = FALSE)
message(sprintf("orthologs mapped: %d | with existing drug: %d",
                length(unique(map$human)), sum(ann$tier == "A: existing drug")))
