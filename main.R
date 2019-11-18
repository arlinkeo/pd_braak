# Main script to run all the analyses
setwd("C:/Users/dkeo/surfdrive/pd_braak/pd_braak")
options(stringsAsFactors = FALSE)

# Useful variables
donorNames <- c("donor9861", "donor10021", "donor12876", "donor14380", "donor15496", "donor15697")
names(donorNames) <- donorNames
braakRoi <- c("R1","R2","R3","R4","R5","R6")
names(braakRoi) <- braakRoi

# Source directory with functions
fun_dir <- paste0(getwd(), "/functions/")
R.utils::sourceDirectory(fun_dir, modifiedOnly = FALSE)

# Make output folder
dir.create("output")

# AHBA Data preprocessing
source("PD/probe2gene_AHBA.R")

# Assign AHBA samples to Braak stage-related regions
source("PD/braakStageSamples.R")

# Select Braak stage-related genes (BRGs)
source("PD/diff_expr_ahba.R")
source("PD/braaklabel_cor.R")
source("PD/braakgenes.R")

# Visualize BRGs in boxplots
source("PD/boxplot_braak.R")

# Gene conversion table for analysis of RNA-seq datasets
source("ahba_entrez2ensembl.R")

# Validate Braak stage related genes in datasets of non-neurological controls
source("PD/diff_expr_UKBEC.R")
source("PD/diff_expr_GTEX.R")

# Validate Braak stage related genes in datasets of PD patients
source("PD/probe2genes_PD_microarray.R")
source("PD/diffexpr_PD_microarray.R")
source("PD/diffexpr_PD_RNAseq.R")

# Braak co-expression modules
source("PD/gene_coexpr.R") # Ran on HP server
source("PD/coexpr_modules.R")
source("PD/eigengene_cor.R")
source("PD/module_enrichment.R")

# Presence of PD-implicated genes
source("PD/pd_implicated_genes.R")

# Other  scripts 
# forestplot.R
# lineplot_dopaminergic_genes.R
# psea.R