# select significant genes based on significant summary estimate

setwd("C:/Users/dkeo/surfdrive/Parkinson")

source("PD/base_script.R")
load("resources/sumEffectSize.RData") # For each gene, a list of tables for each braak stage (donor x statistics)

#Filter for summary effect
bsperGene <- lapply(sumEffectSize, function(g){
  # summary statistics for each braak stage
  t <- t(sapply(g, function(bs){
    bs["SummaryEff", ]
  }))
  t <- as.data.frame(t)
  t
})

#Table of all genes per braak stage
genesPerBs <- lapply(braakNames, function(bs){
  bsTab <- t(sapply(bsperGene, function(g){
    gene <- g[bs, ]
  }))
  bsTab <- as.data.frame(bsTab)
  bsTab
})

#Number of diff. expressed genes based on summary effect p-value
sapply(genesPerBs, function(t){
  sum(unlist(t$pvalue) < 0.05)
})

#Number of diff. expressed genes based on summary effect p-value corrected
sapply(genesPerBs, function(t){
  p <- unlist(t$pvalue)
  correctedP <- p.adjust(p, method = "bonferroni", n = length(p))
  sum(correctedP < 0.05)
})