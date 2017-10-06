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
save(genesPerBs, file = "resources/diffExprBraak.RData")
load("resources/diffExprBraak.RData")

#Number of diff. expressed genes based on summary effect p-value (uncorrected)
sapply(genesPerBs, function(t){
  sum(unlist(t$pvalue) < 0.05)
})

#Average correlation across brains
load("resources/geneLabelCor.RData")
corr <- apply(geneLabelCor, 1, mean)

#Bonferroni-corrected p-values and effect size
pvalEffsize <- lapply(genesPerBs, function(t){
  p <- unlist(t$pvalue)
  names(p) <- rownames(t)
  p <- p.adjust(p, method = "bonferroni", n = length(p))
  data.frame(p, effectSize = unlist(t$meanDiff), r = corr)
})
sapply(pvalEffsize, function(x) sum(x$p<.05 & abs(x$effectSize)>1))
sapply(pvalEffsize, function(x) sum(x$p<.05 & x$effectSize>1))
sapply(pvalEffsize, function(x) sum(x$p<.05 & x$effectSize< -1))
sapply(pvalEffsize, function(x) sum(x$p<.05 & abs(x$effectSize)>1 & abs(x$r)>0.5))
sapply(pvalEffsize, function(x) sum(x$p<.05 & x$effectSize>1 & abs(x$r)>0.5))
sapply(pvalEffsize, function(x) sum(x$p<.05 & x$effectSize< -1 & abs(x$r)>0.5))

