# Gene co-expression in whole Braak region
setwd("C:/Users/dkeo/surfdrive/Parkinson")

library("psych")

source("PD/base_script.R")
load("../ABA_Rdata/BrainExprNorm.RData")
load("resources/braakLabels.RData") # Braak stage label vectors

# Merge labels 1-6
wholeBraak <- lapply(braakLabels, function(d){
  d != 0
})
sapply(wholeBraak, sum)

# Select all Braak samples in expression data
braakExpr <- lapply(donorNames, function(d){
  # Subselect expression matrices
  labels <- wholeBraak[[d]]
  brainExprNorm[[d]][, labels]
})
sapply(braakExpr, dim)

# Gene co-expression in Braak samples, for each donor
gene_coexpr <- lapply(donorNames, function(d){
  expr <- braakExpr[[d]]
  x=corr.test(t(expr), adjust = "BH")
})
save(gene_coexpr, file = "gene_coexpr.Rdata")

