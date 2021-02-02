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
source("pd_braak/probe2gene_AHBA.R")

# Assign AHBA samples to Braak stage-related regions
source("pd_braak/braakStageSamples.R")

# Select Braak stage-related genes (BRGs)
source("pd_braak/diff_expr_ahba.R")
source("pd_braak/braaklabel_cor.R")
source("pd_braak/braakgenes.R")

# Visualize BRGs in boxplots
source("pd_braak/boxplot_braak.R")

# Gene conversion table for analysis of RNA-seq datasets
source("pd_braak/ahba_entrez2ensembl.R")

# Validate Braak stage related genes in datasets of non-neurological controls
source("pd_braak/diff_expr_UKBEC.R")
source("pd_braak/diff_expr_GTEX.R")

# Validate Braak stage related genes in datasets of PD patients
source("pd_braak/probe2genes_PD_microarray.R")
source("pd_braak/diffexpr_PD_microarray.R")
source("pd_braak/diffexpr_PD_RNAseq.R")

# Braak co-expression modules
source("pd_braak/gene_coexpr.R") # Ran on HP server
source("pd_braak/coexpr_modules.R")
source("pd_braak/eigengene_cor.R")
source("pd_braak/module_enrichment.R")

# Presence of PD-implicated genes
source("pd_braak/pd_implicated_genes.R")

# Other  scripts 
source("pd_braak/lineplot_braak.R")
# forestplot.R

# psea.R