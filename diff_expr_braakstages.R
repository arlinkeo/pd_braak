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
      confidence95 <- test2tail$conf.int
      # fChange <- log2(2^test2tail$estimate[1] / 2^test2tail$estimate[2]) # estimate is mean
      # names(fChange) <- 'fold-change'
      c('pvalue' = pval2tail, meanB, meanNB, 'lower95' = confidence95[1], 'upper95' = confidence95[2])
    })))
    pvalues$'pvalue' <- p.adjust(pvalues$'pvalue' , method = "BH", n = nGenes) #corrected p
    pvalues
  })
})

# # convert lists to gene x donor tablefor each braak stage
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
# dgTable2tail <- dg.table("p-value")

save(diffGenesList, file = "resources/diffGenesBraak.RData")
# 
# # Diff. expressed in all donors
# diff.all <- function(list, th) {
#   braakTab <- sapply(list, function(tab){
#     apply(diffUp < th, 1, sum) == nDonors
#   })
#   rownames(braakTab) <- rownames(list[[1]])
#   braakTab
# }
# 
# dgAll2tail <- diff.all(dgTable2tail, 0.05)
# 
# #Combine p-values
# combine.p <- function(list){
#   sapply(list, function(tab){
#     apply(tab, 1, function(p){
#       pchisq(-2*sum(log(p)), df = length(p)*2, lower.tail = FALSE)
#     })
#   })
# }
# 
# fisherP2tail <- combine.p(dgTable2tail)
# 
# #Get list of sorted significant genes
# signif.genes <- function(tab){
#   lapply(braakNames, function(bs){
#     (head(tab[, bs])) < 0.05
#   })
# }
# 
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