# Differential expression in braak regions of PD

options(stringsAsFactors = FALSE)
setwd("C:/Users/dkeo/surfdrive/Parkinson")
load("../ABA_Rdata/BrainExpr.RData")
donorNames <- names(brainExpr)
names(donorNames) <- donorNames
nDonors <- length(donorNames)
load("resources/braakStages.RData")
braakNames <- names(braakStages)
names(braakNames) <- braakNames
probeInfo <- read.csv("../ABA_human_processed/probe_info_2014-11-11.csv")
entrezId2Name <- function (x) {probeInfo$gene_symbol[match(x, probeInfo$entrez_id)]} #Input is vector

#Differential expression for each braak stage compare to non-braak stage
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
      testUp <- t.test(exprB, exprNB, alternative = "greater")
      testDown <- t.test(exprB, exprNB, alternative = "less")
      test2tail <- t.test(exprB, exprNB, alternative = "two.sided")
      pvalUp <- testUp$p.value
      pvalDown <- testDown$p.value
      pval2tail <- test2tail$p.value
      fChange <- log2(test2tail$estimate[1] / test2tail$estimate[2])
      names(fChange) <- NULL
      c(pval_up = pvalUp, pval_down = pvalDown, pval_2tail = pval2tail, 'fold-change' = fChange)
    })))
    pvalues$pval_up <- p.adjust(pvalues$pval_up, method = "BH", n = nGenes)
    pvalues$pval_down <- p.adjust(pvalues$pval_down, method = "BH", n = nGenes)
    pvalues$pval_2tail <- p.adjust(pvalues$pval_2tail , method = "BH", n = nGenes)
    pvalues
  })
})



diffGenesTable <- lapply(diffGenesList, function(tab){
  diffGenesPerBrain <- tab
  # Convert list to two tables with p-values in all donors
  diffUp <- sapply(diffGenesPerBrain, function(t){t[ , "pval_up"]})
  diffDown <- sapply(diffGenesPerBrain, function(t){t[ , "pval_down"]})
  diff2tail <- sapply(diffGenesPerBrain, function(t){t[ , "pval_2tail"]})
  # Select genes signficant in all donors
  sigDiffUp <-apply(diffUp < 0.025, 1, sum) == nDonors
  sigDiffDown <- apply(diffDown < 0.025, 1, sum) == nDonors
  sigDiff2tail <- apply(diff2tail < 0.05, 1, sum) == nDonors
  # Reduce diffUp and diffDown table to only significant genes
  diffUp <- diffUp[sigDiffUp, ]
  diffDown <- diffDown[sigDiffDown, ]
  diff2tail <- diff2tail[sigDiff2tail, ]
  # Sort genes based on lowest p-value found across donors
  sortedUp <- names(sort(apply(diffUp, 1, min)))
  sortedDown <- names(sort(apply(diffDown, 1, min)))
  sorted2tail <- names(sort(apply(diff2tail, 1, min)))
  diffUp <- diffUp[sortedUp, ]
  diffDown <- diffDown[sortedDown, ]
  diff2tail <- diff2tail[sorted2tail, ]
  #Format and add genes in table
  diffUp <- format(diffUp, digits = 4)
  diffDown <- format(diffDown, digits = 4)
  diff2tail <- format(diff2tail, digits = 4)
  geneIdUp <- rownames(diffUp)
  geneIdDown <- rownames(diffDown)
  geneId2tail <- rownames(diff2tail)
  diffUp <- cbind(entrez_id = geneIdUp, gene_symbol = entrezId2Name(geneIdUp), diffUp)
  diffDown <- cbind(entrez_id = geneIdDown, gene_symbol = entrezId2Name(geneIdDown), diffDown)
  diff2tail <- cbind(entrez_id = geneId2tail, gene_symbol = entrezId2Name(geneId2tail), diff2tail)
  #save tables
  write.table(diffUp, file = paste("DiffExpr_braak/diffUp_", "bs", ".txt", sep = ""), quote = FALSE, row.names = FALSE)
  write.table(diffDown, file = paste("DiffExpr_braak/diffDown_", "bs", ".txt", sep = ""), quote = FALSE, row.names = FALSE)
  write.table(diff2tail, file = paste("DiffExpr_braak/diff2tail_", "bs", ".txt", sep = ""), quote = FALSE, row.names = FALSE)
  # Output
  list(upregulated = diffUp, downregulated = diffDown, two_tailed = diff2tail)
})
save(diffGenesList, diffGenesTable, file = "resources/diffGenesBraak.RData")
#############################################
# Examine differential expressed genes

load("resources/diffGenesBraak.RData")

#Genes differentially expressed in all regions
allUp <- lapply(diffGenesBraak2, function(bs){bs[["upregulated"]][ , "gene_symbol"]})
allDown <- lapply(diffGenesBraak2, function(bs){bs[["downregulated"]][ , "gene_symbol"]})
allStagesUp <- Reduce(intersect, allUp)
allStagesDown <- Reduce(intersect, allDown)
allStagesUp
allStagesDown
# 
# # Comparing top 10 genes for each region to other regions
# top10up <- lapply(diffGenesBraak2, function(bs){
#   bs[["upregulated"]][1:10, "gene_symbol"]
# })
# top10down <- lapply(diffGenesBraak2, function(bs){
#   bs[["downregulated"]][1:10, "gene_symbol"]
# })
# lapply(top10up, is.present.up)
# lapply(top10down, is.present.down)