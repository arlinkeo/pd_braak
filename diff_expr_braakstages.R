# Differential expression in braak regions of PD

setwd("C:/Users/dkeo/surfdrive/Parkinson")
library("metap")

source("PD/base_script.R")
load("../ABA_Rdata/BrainExprNorm.RData")
nDonors <- length(donorNames)
load("resources/braakStages.RData")

nGenes <- nrow(brainExprNorm$donor9861)

# T-test for each gene
ttestGene <- function(b, nB) {# Braak vs. non-braak
  test2tail <- t.test(b, nB, alternative = "two.sided")
  estimate <- test2tail$estimate
  confidence95 <- test2tail$conf.int
  # fChange <- log2(2^test2tail$estimate[1] / 2^test2tail$estimate[2]) # estimate is mean
  # names(fChange) <- 'fold-change'
  c('pvalue' = test2tail$p.value, 
    'meanB' = estimate[1], 'varB' = var(unlist(b)), 
    'meanNB' = estimate[2], 'varNB' = var(unlist(nB)), 
    'lower95' = confidence95[1], 'upper95' = confidence95[2])
}

# T-test for all genes in matrix
ttestTab <- function(df, nV, nBV){ # expr. data and binary vectors to select columns/samples
  tab <- as.data.frame(t(apply(df, 1, function(g){
    b <- g[nV]
    nB <- g[nBV]
    ttestGene(b, nB)
  })))
  tab$'benjamini_hochberg' <- p.adjust(tab$'pvalue' , method = "BH", n = nGenes) #corrected p
  tab$'bonferroni' <- p.adjust(tab$'pvalue' , method = "bonferroni", n = nGenes) #corrected p
  tab
}

#Combine table of braak and merged braak 
braakStages2 <- lapply(donorNames, function(d){
  braak <- braakStages[[d]]
  braakM1 <- braakMerged1[[d]]
  braakM2 <- braakMerged2[[d]]
  cbind(braak, braakM1, braakM2)
})
braakNames2 <- c(braakNames, braakNamesMerged1, braakNamesMerged2)

diffExprRef <- lapply(names(nonBraak), function(r){# Diff. expr. for each non-Braak reference
  ref <- nonBraak[[r]]
  diffExprPerBrain <- lapply(donorNames, function(d){# Diff. expr. for each donor
    expr <- brainExprNorm[[d]]
    braak <- braakStages2[[d]]
    nonBraakVec <- as.logical(ref[[d]])
    lapply(braakNames2, function(bs){ # Diff. expr. for each braak stage
      print(paste0(r, ", ", d, ", ", bs)) # print progress
      braakVec <- as.logical(braak[, bs])
      pvalues <- ttestTab(expr, braakVec, nonBraakVec)
    })
  })
})
names(diffExprRef) <- names(nonBraak)

save(diffExprRef, file = "resources/diffExprBraak.RData")

lapply(diffExprRef, function(ref){
  t(sapply(ref, function(d){
    sapply(d, function(b){
      sum(b[, "bonferroni"] < 0.05)
    })
  }))
})

##################################################################################
# #Differential expression for each braak stage compared to non-braak stage
# diffGenesList <- lapply(braakNames, function(bs){
#   braakInfo <- braakStages[[bs]]
#   # Differential expression for each donor
#   diffGenesPerBrain <- lapply(donorNames, function(d){
#     expr <- brainExpr[[d]]
#     nGenes <- nrow(expr)
#     braakVec <- as.logical(braakInfo[[d]])
#     nonBraakVec <- as.logical(nonBraak[[d]])
#     # T-test for each gene
#     pvalues <- as.data.frame(t(apply(expr, 1, function(g){
#       exprB <- g[braakVec]
#       exprNB <- g[nonBraakVec]
#       test2tail <- t.test(exprB, exprNB, alternative = "two.sided")
#       pval2tail <- test2tail$p.value
#       meanB <- test2tail$estimate[1]
#       meanNB <- test2tail$estimate[2]
#       names(meanB) <- 'meanB'
#       names(meanNB) <- 'meanNB'
#       varB <- var(unlist(exprB))
#       varNB <- var(unlist(exprNB))
#       confidence95 <- test2tail$conf.int
#       # fChange <- log2(2^test2tail$estimate[1] / 2^test2tail$estimate[2]) # estimate is mean
#       # names(fChange) <- 'fold-change'
#       c('pvalue' = pval2tail, meanB, 'varB' = varB, meanNB, 'varNB' = varNB, 'lower95' = confidence95[1], 'upper95' = confidence95[2])
#     })))
#     pvalues$'benjamini_hochberg' <- p.adjust(pvalues$'pvalue' , method = "BH", n = nGenes) #corrected p
#     pvalues$'bonferroni' <- p.adjust(pvalues$'pvalue' , method = "bonferroni", n = nGenes) #corrected p
#     pvalues
#   })
# })
# save(diffGenesList, file = "resources/diffGenesBraak.RData")
# 
# # Number of diff. expr. genes in each braak stage en donor after multiple testing correction
# nDiffgenes <- function(x){
#   sapply(diffGenesList, function(bs){
#     sapply(bs, function(d){
#       sum(d[[x]] < 0.05)
#     })
#   })
# }
# 
# nDiffgenes("bonferroni")
# nDiffgenes("benjamini_hochberg")
# 
# nDonors <- 6
# 
# # convert lists to gene x donor table for each braak stage
# dg.table <- function(x){
#   lapply(braakNames, function(bs){
#     brainList <- diffGenesList[[bs]]
#     genes <- rownames(brainList[[1]])
#     tab <- sapply(brainList, function(t){v <- t[ , x]})
#     rownames(tab) <- genes
#     tab
#   })
# }
# 
# dgTable1 <- dg.table("benjamini_hochberg")
# dgTable2 <- dg.table("bonferroni")
##################################################################################