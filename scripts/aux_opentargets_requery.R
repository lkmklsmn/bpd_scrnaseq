# 09b - Druggability annotation via babelgene + Open Targets Platform
#       (GraphQL v4).
# Updated re-implementation of 09a for the CURRENT Open Targets schema.
#
# Why this script exists
# ----------------------
# 09a_opentargets_query.R was written against the Open Targets GraphQL schema
# as of 30 June 2026. That schema has since changed and 09a no longer runs:
#   * Target.knownDrugs -> removed; replaced by
#     Target.drugAndClinicalCandidates
#   * KnownDrug.drug.maxPhaseForIndication -> removed; use
#     ClinicalTargetFromTarget.maxClinicalStage
# This script queries the current fields, keeps 09a's tier logic unchanged,
# and additionally records subcellular location and target class (used to
# argue accessibility of the target).
#
# API VERSION / ACCESS DATE
# -------------------------
# Endpoint : https://api.platform.opentargets.org/api/v4/graphql
# Accessed : 2026-08-29
# The API is version-dependent. Every response is cached as JSON under
# <external_dir>/opentargets_2026-08-29/ so that the annotation can be
# reproduced exactly
# even if the live schema changes again. Delete the cache to re-query.
#
# Output: results/tables/drug_target_candidates_curated_2026-08-29.csv
#         (written alongside, NOT over,
#         09a's drug_target_candidates_curated.csv)
#
# VERIFY against the paper: 262 AM-upregulated genes -> 257 human orthologs
# mapped, 75 with an approved/clinical-stage drug.

source(file.path("config", "paths.R"))
suppressMessages({
  library(readxl)
  library(babelgene)
  library(httr)
  library(jsonlite)
})

ot_url <- "https://api.platform.opentargets.org/api/v4/graphql"
access_date <- "2026-08-29"
cache_dir <- file.path(external_dir, paste0("opentargets_", access_date))
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ---- 1. AM-upregulated genes (padj < 0.25, log2FC > 0) ----------------------
am <- as.data.frame(read_excel(deseq_xlsx, sheet = "Alveolar macrophages"))
up_mouse <- am$gene[!is.na(am$padj) & am$padj < 0.25 & am$log2FoldChange > 0]
message(sprintf("AM-upregulated mouse genes: %d", length(up_mouse)))

# ---- 2. mouse -> human orthologs --------------------------------------------
orth <- orthologs(genes = up_mouse, species = "mouse", human = FALSE)
map <- data.frame(
  mouse = orth$symbol,
  human = orth$human_symbol,
  stringsAsFactors = FALSE
)
map <- map[!is.na(map$human) & map$human != "", ]
message(sprintf("human orthologs mapped: %d", length(unique(map$human))))

# ---- 3. GraphQL helpers (cache-first) ---------------------------------------
gql <- function(query, variables) {
  res <- POST(
    ot_url,
    body = list(query = query, variables = variables), encode = "json"
  )
  content(res, as = "parsed", type = "application/json")
}

q_search <- paste(
  "query($q:String!){",
  '  search(queryString:$q, entityNames:["target"]){ hits{ id name } }',
  "}"
)

# CURRENT schema: drugAndClinicalCandidates replaces knownDrugs
q_target <- "query($id:String!){ target(ensemblId:$id){
               approvedSymbol
               targetClass{ label }
               tractability{ modality value label }
               subcellularLocations{ location }
               drugAndClinicalCandidates{
                 count rows{ maxClinicalStage drug{ name drugType } }
               } } }"

resolve_ensembl <- function(symbol) {
  f <- file.path(cache_dir, paste0("search_", symbol, ".json"))
  if (file.exists(f)) {
    return(fromJSON(f, simplifyVector = FALSE)$id)
  }
  hits <- tryCatch(
    gql(q_search, list(q = symbol))$data$search$hits,
    error = function(e) NULL
  )
  # prefer an exact symbol match over the top-ranked hit
  eid <- NA_character_
  if (length(hits)) {
    exact <- Filter(
      function(h) toupper(h$name %||% "") == toupper(symbol),
      hits
    )
    eid <- if (length(exact)) exact[[1]]$id else hits[[1]]$id
  }
  write_json(list(id = eid), f, auto_unbox = TRUE)
  eid
}

target_info <- function(ensembl_id) {
  f <- file.path(cache_dir, paste0(ensembl_id, ".json"))
  if (file.exists(f)) {
    return(fromJSON(f, simplifyVector = FALSE))
  }
  out <- gql(q_target, list(id = ensembl_id))
  write_json(out, f, auto_unbox = TRUE)
  Sys.sleep(0.34) # be polite to the public API
  out
}

# ---- 4. classify each ortholog (tier logic identical to 09a) -----------------
classify <- function(symbol) {
  eid <- resolve_ensembl(symbol)
  if (is.na(eid)) {
    return(data.frame(
      human_symbol = symbol, ensembl_id = NA, tier = "Not druggable",
      example_drugs = "", max_clinical_stage = "", secreted = NA,
      subcellular = "", target_class = "", stringsAsFactors = FALSE
    ))
  }

  ti <- tryCatch(target_info(eid)$data$target, error = function(e) NULL)

  drugs <- character(0)
  stages <- character(0)
  rows <- ti$drugAndClinicalCandidates$rows
  if (!is.null(rows)) {
    drugs <- unique(vapply(rows, function(r) r$drug$name %||% "", character(1)))
    drugs <- drugs[drugs != ""]
    stages <- unlist(lapply(rows, function(r) r$maxClinicalStage %||% NULL))
  }
  tract <- if (!is.null(ti$tractability)) {
    vapply(ti$tractability, function(t) isTRUE(t$value), logical(1))
  } else {
    logical(0)
  }

  # unchanged from 09a
  tier <- if (length(drugs)) {
    "A: existing drug"
  } else if (any(tract)) {
    "B: clinically tractable"
  } else {
    "Druggable"
  }

  locs <- unlist(lapply(ti$subcellularLocations, function(l) {
    l$location %||% NULL
  }))
  locs <- unique(locs[!is.na(locs)])
  sec <- any(grepl("secret|extracellular", locs, ignore.case = TRUE))
  tcl <- unlist(lapply(ti$targetClass, function(t) t$label %||% NULL))

  data.frame(
    human_symbol = symbol,
    ensembl_id = eid,
    tier = tier,
    example_drugs = paste(head(drugs, 3), collapse = ";"),
    max_clinical_stage = if (length(stages)) {
      sort(stages, decreasing = TRUE)[1]
    } else {
      ""
    },
    secreted = sec,
    subcellular = paste(head(locs, 4), collapse = ";"),
    target_class = paste(head(tcl, 2), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

ann <- do.call(rbind, lapply(unique(map$human), classify))

# ---- 5. assemble curated table ----------------------------------------------
ann$mouse_genes <- vapply(ann$human_symbol, function(h) {
  paste(unique(map$mouse[map$human == h]), collapse = "/")
}, character(1))
ann$celltypes <- "Alveolar macrophages"

# Counts are reported at TWO stages, matching the manuscript:
#   (i)  classification stage - all mapped orthologs. The paper's "75 were
#        targeted by an approved drug or clinical-stage compound" is this
#        number.
#   (ii) ranking stage - after excluding non-specific housekeeping genes
#        (ribosomal / mitochondrial), which is what Fig 6B ranks. Most
#        tier-A genes removed here are ribosomal proteins that carry
#        incidental drug annotations (e.g. ataluren).
n_tier_a_all <- sum(ann$tier == "A: existing drug")
n_tier_b_all <- sum(ann$tier == "B: clinically tractable")
housekeeping <- grepl("^RP[LS]|^MRP|^MT-", ann$human_symbol)
message(sprintf(
  "classification stage: %d orthologs | tier A: %d | tier B: %d",
  nrow(ann), n_tier_a_all, n_tier_b_all
))
message(sprintf(
  paste(
    "housekeeping (ribosomal/mitochondrial) removed for ranking:",
    "%d (of which tier A: %d)"
  ),
  sum(housekeeping), sum(housekeeping & ann$tier == "A: existing drug")
))

# drop non-specific housekeeping (ribosomal / mitochondrial) as in Methods
ann <- ann[!housekeeping, ]

out_file <- file.path(
  tab_dir,
  paste0("drug_target_candidates_curated_", access_date, ".csv")
)
write.csv(
  ann[, c(
    "mouse_genes", "human_symbol", "ensembl_id", "celltypes", "tier",
    "example_drugs", "max_clinical_stage", "secreted", "subcellular",
    "target_class"
  )],
  out_file,
  row.names = FALSE
)

message(sprintf(
  "ranking stage (housekeeping removed): %d targets | tier A: %d | tier B: %d",
  nrow(ann),
  sum(ann$tier == "A: existing drug"),
  sum(ann$tier == "B: clinically tractable")
))
message(paste(
  "EXPECTED (paper): 262 AM-up genes -> 257 orthologs -> 75 with an",
  "existing drug at the classification stage"
))
message(sprintf(
  "Open Targets endpoint %s, accessed %s; responses cached in %s",
  ot_url, access_date, cache_dir
))
message(sprintf("written: %s", out_file))
