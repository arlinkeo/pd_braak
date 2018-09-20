setwd("C:/Users/dkeo/surfdrive/pd_braak")

# Prepare Braak data
source("PD/braakStageSamples.R")

# # Cell-type markers
# source("PD/lm_celltypes.R")
# source("PD/diff_expr_lm.R")

# Single gene analysis
source("PD/diff_expr.R")
source("PD/braaklabel_cor.R")
source("PD/select_genes.R")


source("PD/funcEnrich.R")
source("PD/diff_expr_hemispheres.R")

# Visualize data
source("PD/boxplot_braak.R")

# Gene conversion table
source("ahba_entrez2ensembl.R")

# Validate progression genes
source("PD/diff_expr_GTEX.R")
source("PD/ukbecGeneID.R")
source("PD/diff_expr_UKBEC.R")
source("PD/boxplot_UKBEC.R")

# Braak co-expression modules
source("PD/gene_coexpr.R")
source("PD/coexpr_modules.R")
source("PD/eigengene_cor.R")
source("PD/module_enrichment.R")
# source("PD/")

# source("PD/")