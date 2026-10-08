# lintr configuration for the bpd_scrnaseq reproducibility pipeline.
#
# Default tidyverse linters apply (80-column lines, snake_case objects,
# no semicolons, 2-space indentation), with one exception:
#
# commented_code_linter is disabled. Every instance it flags in this repo is a
# provenance or input-path comment in a script header that happens to parse as
# R code, for example:
#
#   #   GSE32472   : data/external/GSE32472/gse32472.RData   -> parses as `:`
#   # Source: figures_for_paper_dblt_removed.R (des_deseq + plot_de)  -> a call
#
# These comments record where every input file and every figure comes from,
# which is the point of this repository, so they are kept deliberately.
#
# Check the repo with:  lintr::lint_dir(".")
linters <- linters_with_defaults(
  commented_code_linter = NULL
)
