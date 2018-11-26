setwd("C:/Users/dkeo/surfdrive/pd_braak")

# Data preprocessing
source("../AHBA_Arlin/probe2gene.R")

# Prepare Braak data
source("PD/braakStageSamples.R")

# Braak stage related genes (not corrected for cell-types)
source("PD/diff_expr.R")
source("PD/braaklabel_cor.R")
source("PD/braakgenes.R")
source("PD/diff_expr_hemispheres.R")

#  Braak stage related genes (corrected for cell-types)
# source("PD/lm_celltypes.R")
# source("PD/diff_expr_compare_references.R")
source("PD/diff_expr_lm.R")
# source("PD/diff_expr_lm_eigengene.R")

# Visualize data
source("PD/boxplot_braak.R")

# Gene conversion table
source("ahba_entrez2ensembl.R")

# Validate Braak stage related genes
source("PD/diff_expr_GTEX.R")
source("PD/diff_expr_UKBEC.R")

# Braak co-expression modules
source("PD/gene_coexpr.R")
source("PD/coexpr_modules.R")
source("PD/eigengene_cor.R")
source("PD/module_enrichment.R")
