# Gene co-expression in whole Braak region
setwd("/tudelft.net/staff-bulk/ewi/insy/DBL/Arlin/pd_braak")

source("PD/base_script.R")

# Random subselection of genes
genes <- ahba.genes()

# Labels for whole Braak region
load("resources/braakLabels.RData") # Braak stage label vectors
wholeBraak <- lapply(braakLabels, function(d){ d != 0 })# Merge labels 1-6

# Labels for non-Braak region
load("resources/nonBraak.RData") # Braak stage label vectors
nonBraak <- nonBraak$nonBraakA2 # all except braak and cerebellum samples

# Labels for each Braak region 
load("resources/braakStages.RData")
braakStages <- lapply(braakNames, function(b){
  lapply(braakStages, function(d) unlist(d[, b]))
})

# Combine all labels of regions in list
sampleList <- c(braakStages, 'braak1-6' = list(wholeBraak), nonbraak = list(nonBraak))

# Average co-expression for all regions in sampleList
lapply(names(sampleList), function(n){
  r <- sampleList[[n]]
  expr <- select.expr(genes, r)
  gene_coexpr <- gene.coexpr(expr)
  
  avgCoexpr <- apply(simplify2array(gene_coexpr), 1:2, mean)
  save(avgCoexpr, file = paste0("resources/avgCoexpr_", n, ".RData"))
})