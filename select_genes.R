# select significant genes based on significant summary estimate

setwd("C:/Users/dkeo/surfdrive/Parkinson")

source("PD/base_script.R")
# load("resources/sumEffectSize.RData") # For each gene, a list of tables for each braak stage (donor x statistics)
load("resources/summaryEffect.RData") # new
load("resources/summaryCorr.RData")

#Filter for summary effect
summaryMeanDiff <- lapply(summaryEffect, function(ref){
  lapply(ref, function(g){
    as.data.frame(t(sapply(g, function(bs){
      bs["summary", ]
    })))
  })
})
summaryCorr2 <- as.data.frame(t(sapply(summaryCorr, function(g){
  unlist(g["summary", ])
})))


# #Table of all genes per braak stage
# genesPerBs <- lapply(braakNames, function(bs){
#   bsTab <- t(sapply(bsperGene, function(g){
#     gene <- g[bs, ]
#   }))
#   bsTab <- as.data.frame(bsTab)
#   bsTab
# })
# save(genesPerBs, file = "resources/diffExprBraak.RData")
# load("resources/diffExprBraak.RData")
# 
# #Number of diff. expressed genes based on summary effect p-value (uncorrected)
# sapply(genesPerBs, function(t){
#   sum(unlist(t$pvalue) < 0.05)
# })
# 
# #Average correlation across brains
# load("resources/geneLabelCor.RData")
# corr <- apply(geneLabelCor, 1, mean)
# 
# #Bonferroni-corrected p-values and effect size
# pvalEffsize <- lapply(genesPerBs, function(t){
#   p <- unlist(t$pvalue)
#   names(p) <- rownames(t)
#   p <- p.adjust(p, method = "bonferroni", n = length(p))
#   data.frame(p, effectSize = unlist(t$meanDiff), r = corr)
# })
# # sapply(pvalEffsize, function(x) sum(x$p<.05 & abs(x$effectSize)>1))
# # sapply(pvalEffsize, function(x) sum(x$p<.05 & x$effectSize>1))
# # sapply(pvalEffsize, function(x) sum(x$p<.05 & x$effectSize< -1))
# sapply(pvalEffsize, function(x) sum(x$p<.05 & abs(x$effectSize)>1 & abs(x$r)>0.5))
# sapply(pvalEffsize, function(x) sum(x$p<.05 & x$effectSize>1 & abs(x$r)>0.5))
# sapply(pvalEffsize, function(x) sum(x$p<.05 & x$effectSize< -1 & abs(x$r)>0.5))
# 
# #Diff. expr. in 1-3 OR 4-6 OR both
# genes <- sapply(pvalEffsize, function(x) rownames(x)[which(x$p<.05 & abs(x$effectSize)>1 & abs(x$r)>0.5)])
# genes1to3 <- setdiff(genes$`braak1-3`, genes$`braak4-6`) #significant in 1-3 and not 4-6
# genes4to6 <- setdiff(genes$`braak4-6`, genes$`braak1-3`)
# genes1to6 <- intersect(genes$`braak4-6`, genes$`braak1-3`)
# upGenes <- sapply(pvalEffsize, function(x) rownames(x)[which(x$p<.05 & x$effectSize>1 & abs(x$r)>0.5)])
# upGenes1to3 <- setdiff(upGenes$`braak1-3`, upGenes$`braak4-6`)
# upGenes4to6 <- setdiff(upGenes$`braak4-6`, upGenes$`braak1-3`)
# upGenes1to6 <- intersect(upGenes$`braak4-6`, upGenes$`braak1-3`)
# totalUpGenes1to6 <- union(upGenes$`braak4-6`, upGenes$`braak1-3`)
# downGenes <- sapply(pvalEffsize, function(x) rownames(x)[which(x$p<.05 & x$effectSize< -1 & abs(x$r)>0.5)])
# downGenes1to3 <- setdiff(downGenes$`braak1-3`, downGenes$`braak4-6`)
# downGenes4to6 <- setdiff(downGenes$`braak4-6`, downGenes$`braak1-3`)
# downGenes1to6 <- intersect(downGenes$`braak4-6`, downGenes$`braak1-3`)
# totalDownGenes1to6 <- union(downGenes$`braak4-6`, downGenes$`braak1-3`)
# relatedGenes1to3 <- c(upGenes1to3, downGenes1to3)
# relatedGenes4to6 <- c(upGenes4to6, downGenes4to6)
# relatedGenes1to6 <- c(upGenes1to6, downGenes1to6)
# relatedGenes <- list(braak1to3 = relatedGenes1to3, braak4to6 = relatedGenes4to6, braak1to6 = relatedGenes1to6)
# save(relatedGenes, file = "braakRelatedGenes.RData")
# 
# #Presence of PD-implicated genes
# lapply(relatedGenes, function(b){
#   overlap <- lapply(pdGenesID, function(pd){intersect(pd, b)})
#   sapply(overlap, entrezId2Name)
# })
#   
# Get correlated genes |r>0.5|
corrGenes <- rownames(summaryCorr2)[abs(as.numeric(summaryCorr2$z)) > 0.5]

# save(corrGenes, file = "resources/correlated_genes1.RData")
# load("resources/correlated_genes1.RData")

# Correlated genes with diff. expression in one of the braak regions
braakGenes <- lapply(summaryMeanDiff, function(ref){
  tabll <- ref[corrGenes]
  genes <- sapply(tabll, function(g){
    any(g$meanDiff > 1 & g$pvalue <0.05)
  })
  names(tabll)[genes]
})
overlapBraakGenes <- Reduce(intersect, braakGenes)
sapply(braakGenes, function(x) any(x== "6622"))
 