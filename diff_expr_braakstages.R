# Differential expression in braak regions of PD

setwd("C:/Users/dkeo/surfdrive/Parkinson")
library("metap")

source("PD/base_script.R")
load("../ABA_Rdata/BrainExpr.RData")
nDonors <- length(donorNames)
load("resources/braakStages.RData")
braakNames <- names(braakStages)
names(braakNames) <- braakNames

#Differential expression for each braak stage compared to non-braak stage
diffGenesList <- lapply(braakNames, function(bs){
  braakInfo <- braakStages[[bs]]
  # Differential expression for each donor
  diffGenesPerBrain <- lapply(donorNames, function(d){
    expr <- brainExpr[[d]]
    nGenes <- nrow(expr)
    braakVec <- as.logical(braakInfo[[d]])
    nonBraakVec <- as.logical(nonBraak[[d]])
    # T-test for each gene
    pvalues <- as.data.frame(t(apply(expr, 1, function(g){
      exprB <- g[braakVec]
      exprNB <- g[nonBraakVec]
      test2tail <- t.test(exprB, exprNB, alternative = "two.sided")
      pval2tail <- test2tail$p.value
      meanB <- test2tail$estimate[1]
      meanNB <- test2tail$estimate[2]
      names(meanB) <- 'meanB'
      names(meanNB) <- 'meanNB'
      varB <- var(unlist(exprB))
      varNB <- var(unlist(exprNB))
      confidence95 <- test2tail$conf.int
      # fChange <- log2(2^test2tail$estimate[1] / 2^test2tail$estimate[2]) # estimate is mean
      # names(fChange) <- 'fold-change'
      c('pvalue' = pval2tail, meanB, 'varB' = varB, meanNB, 'varNB' = varNB, 'lower95' = confidence95[1], 'upper95' = confidence95[2])
    })))
    pvalues$'benjamini_hochberg' <- p.adjust(pvalues$'pvalue' , method = "BH", n = nGenes) #corrected p
    pvalues$'bonferroni' <- p.adjust(pvalues$'pvalue' , method = "bonferroni", n = nGenes) #corrected p
    pvalues
  })
})

save(diffGenesList, file = "resources/diffGenesBraak.RData")

# Number of diff. expr. genes in each braak stage en donor after multiple testing correction
nDiffgenes <- function(x){
  sapply(diffGenesList, function(bs){
    sapply(bs, function(d){
      sum(d[[x]] < 0.05)
    })
  })
}

nDiffgenes("bonferroni")
nDiffgenes("benjamini_hochberg")

# #############################################
# # Examine differential expressed genes
# 
# load("resources/diffGenesBraak.RData")
# 
# #Genes differentially expressed in all regions
# allUp <- lapply(diffGenesBraak2, function(bs){bs[["upregulated"]][ , "gene_symbol"]})
# allDown <- lapply(diffGenesBraak2, function(bs){bs[["downregulated"]][ , "gene_symbol"]})
# allStagesUp <- Reduce(intersect, allUp)
# allStagesDown <- Reduce(intersect, allDown)
# allStagesUp
# allStagesDown
# # 
# # # Comparing top 10 genes for each region to other regions
# # top10up <- lapply(diffGenesBraak2, function(bs){
# #   bs[["upregulated"]][1:10, "gene_symbol"]
# # })
# # top10down <- lapply(diffGenesBraak2, function(bs){
# #   bs[["downregulated"]][1:10, "gene_symbol"]
# # })
# # lapply(top10up, is.present.up)
# # lapply(top10down, is.present.down)