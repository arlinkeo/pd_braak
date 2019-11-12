# Main script to run all the analyses

# Useful variables
donorNames <- c("donor9861", "donor10021", "donor12876", "donor14380", "donor15496", "donor15697")
names(donorNames) <- donorNames
braakRoi <- c("R1","R2","R3","R4","R5","R6")
names(braakRoi) <- braakRoi

# AHBA data directory and data
ahba_dir <-"C:/Users/dkeo/surfdrive/AHBA_Arlin"
probeInfo <- read.csv(paste0(ahba_dir, "/probe_info_2018-11-18.csv"))
brainExpr <- readRDS(paste0(ahba_dir, "/gene_expr.RDS"))
ontology <- read.csv(paste0(ahba_dir, "/Ontology.csv"))
sample_annot <- lapply(donorNames, function(d){ # Sample info per donor
  t <- read.csv(paste0(ahba_dir, "/sample_info_", d, "_2018-11-18.csv"))
  info <- ontology[match(t$structure_id, ontology$id), ]
  info$color_hex_triplet <- paste0("#", info$color_hex_triplet)
  info
})

# Source directory with functions
fun_dir <- paste0(getwd(), "/functions/")
R.utils::sourceDirectory(fun_dir, modifiedOnly = FALSE)

# Make output folder
dir.create("output")

# Data preprocessing
# source(paste0(ahba_dir, "/probe2gene.R"))

# Run scripts
source("PD/braakStageSamples.R")
source("PD/diff_expr_ahba.R")
source("PD/braaklabel_cor.R")

source("PD/braakgenes.R")

# Visualize data
source("PD/boxplot_braak.R")

# Gene conversion table for analysis of RNA-seq datasets
source("ahba_entrez2ensembl.R")

# Validate Braak stage related genes
source("PD/diff_expr_UKBEC.R")
source("PD/diff_expr_GTEX.R")

# Braak co-expression modules
source("PD/gene_coexpr.R") # Ran on HP server
source("PD/coexpr_modules.R")
source("PD/eigengene_cor.R")
source("PD/module_enrichment.R")

# Presence of PD-implicated genes
source("PD/pd_implicated_genes.R")