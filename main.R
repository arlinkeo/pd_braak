setwd("C:/Users/dkeo/surfdrive/pd_braak")

# Prepare data
source("PD/roiSamples.R")
source("PD/braakStageSamples.R")

# Select progression genes
source("PD/braaklabel_cor.R")
source("PD/diff_expr.R")
source("PD/select_genes.R")
source("PD/diff_expr_hemispheres.R")

# Visualize data
source("PD/volcanoplot_summary_diff.R")
source("PD/boxplot_braak.R")

# Gene conversion table
source("ahba_entrez2ensembl.R")

# Validate progression genes
source("PD/diff_expr_GTEX.R")
source("PD/diff_expr_UKBEC.R")
source("PD/boxplot_UKBEC.R")

# Braak co-expression modules
source("PD/coexpr_modules.R")
source("PD/eigengene_cor.R")

source("PD/funcEnrich.R")
source("PD/module_enrichment.R")
# source("PD/")
# source("PD/")
# source("PD/")
# source("PD/")