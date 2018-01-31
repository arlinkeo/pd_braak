# Gene co-expression in whole Braak region
setwd("/tudelft.net/staff-bulk/ewi/insy/DBL/Arlin/pd_braak")

library("Hmisc")

source("PD/base_script.R")
load("../ABA_Rdata/BrainExprNorm.RData")
load("resources/nonBraak.RData") # Braak stage label vectors

# Select non-Braak region
nonBraak <- nonBraak$nonBraakA2 # all except braak and cerebellum samples

# Select all Braak samples in expression data
braakExpr <- lapply(donorNames, function(d){
  # Subselect brain samples
  labels <- as.logical(nonBraak[[d]])
  brainExprNorm[[d]][, labels]
})
sapply(braakExpr, dim)
remove(brainExprNorm)

# p-values of correlations
gene_coexpr_nonbraak <- lapply(donorNames, function(d){
  expr <- braakExpr[[d]]
  m <- rcorr(t(expr), type = "pearson") #List with correlation-, size-, and p-value-matrix
  m$n <- m$n[1, 1]
  m
})
save(gene_coexpr_nonbraak, file = "resources/gene_coexpr_nonbraak.Rdata")