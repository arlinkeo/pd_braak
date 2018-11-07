setwd("C:/Users/dkeo/surfdrive/pd_braak")

# Prepare Braak data
source("PD/braakStageSamples.R")

# Single gene analysis
source("PD/diff_expr.R")
source("PD/braaklabel_cor.R")
source("PD/braakgenes.R")
source("PD/diff_expr_hemispheres.R")

# Visualize data
source("PD/boxplot_braak.R")

# Braak genes corrected for cell-type expression
# source("PD/lm_celltypes.R")
source("PD/diff_expr_lm.R")

# Gene conversion table
source("ahba_entrez2ensembl.R")

# Validate Braak stage related genes
source("PD/diff_expr_GTEX.R")
# source("PD/ukbecGeneID.R")
source("PD/diff_expr_UKBEC.R")

# Braak co-expression modules
source("PD/gene_coexpr.R")
source("PD/coexpr_modules.R")
source("PD/eigengene_cor.R")
source("PD/module_enrichment.R")
