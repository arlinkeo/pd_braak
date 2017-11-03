# select significant genes based on significant summary estimate

setwd("C:/Users/dkeo/surfdrive/Parkinson")

library("metap")

source("PD/base_script.R")
load("resources/summaryMeanDiff.RData")
load("resources/summaryCorr.RData")

#Filter for summary effect
summaryMeanDiff <- lapply(summaryMeanDiff, function(ref){
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

# Diff. expressed genes in one of the merged braak regions for each non-Braak region
diffGenes <- lapply(summaryMeanDiff, function(ref){
  genes <- sapply(ref, function(g){
    g <- g[braakNamesMerged1, ] # check only in merged braak
    any(abs(unlist(g$meanDiff)) > 1 & g$benjamini_hochberg <0.05)
  })
  names(genes)[genes]
})

similarity <- sapply(diffGenes, function(a){ #jaccard
  sapply(diffGenes, function(b){
    shared <- length(intersect(a, b))
    unshared <- length(union(a, b))
    shared/unshared*100
  })
})#[c(1,3,5,2,4,6), c(1,3,5,2,4,6)]

overlapRef <- Reduce(intersect, diffGenes)
overlapRefA <- Reduce(intersect, diffGenes[c(1,3,5)])
overlapRefB <- Reduce(intersect, diffGenes[c(2,4,6)])

# Correlated genes with diff. expression in one of the braak regions
braakGenesPerRef <- lapply(diffGenes, function(ref){
  intersect(ref, corrGenes)
})
braakGenes <- braakGenesPerRef$nonBraakA2

save(braakGenes, file = "resources/braakGenes.RData")
any(braakGenes== "6622")

braakGenesList <- summaryMeanDiff$nonBraakA2[braakGenes]
#Check which Braak regions show diff. expression in 1-2, 3-4, and/or 5-6
braakGenesList2 <- lapply(braakGenesList, function(x)x[braakNamesMerged1, ])
braakProfile <- sapply(braakGenesList2, function(t) {
  diff <- unlist(t$meanDiff)
  pval <- unlist(t$benjamini_hochberg)
  profile <- as.numeric(abs(diff) > 1 & pval < 0.05)
  neg <- ifelse(diff <0, -1, 1)
  profile <- neg*profile
  paste(as.character(profile), collapse = "")
})
occurences <- as.data.frame(table(braakProfile))
profiles <- as.character(occurences$braakProfile)
names(profiles) <- profiles

profileBraakGenes <- lapply(profiles, function(p){
  idx <- braakProfile %in% p
  names(braakProfile)[idx]
})
save(profileBraakGenes, file = "resources/profileBraakGenes.RData")

# Correlation of Braak profile genes
profileCorr <- lapply(profileBraakGenes, function(p){
  summaryCorr2[p, ]
})

#Presence of PD-implicated genes
lapply(profileBraakGenes, function(p){
  overlap <- lapply(pdGenesID, function(pd){intersect(pd, p)})
  sapply(overlap, entrezId2Name)
})

