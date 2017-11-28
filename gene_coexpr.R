# Gene co-expression in whole Braak region
setwd("C:/Users/dkeo/surfdrive/Parkinson")

library("Hmisc")

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
remove(brainExprNorm)

# p-values of correlations
gene_coexpr <- lapply(donorNames, function(d){
  expr <- braakExpr[[d]]
  m <- rcorr(t(expr), type = "pearson") #List with correlation-, size-, and p-value-matrix
  m$n <- m$n[1, 1]
  m
})
save(gene_coexpr, file = "gene_coexpr.Rdata")