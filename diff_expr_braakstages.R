# Differential expression in braak regions of PD

options(stringsAsFactors = FALSE)
setwd("C:/Users/dkeo/surfdrive/Parkinson")
load("../polyQ_coexpression/resources/BrainExpr.RData")
donorNames <- names(brainExpr)
names(donorNames) <- donorNames
nDonors <- length(donorNames)
load("resources/braakStages.RData")
braakNames <- names(braakStages)
names(braakNames) <- braakNames
probeInfo <- read.csv("../ABA_human_processed/probe_info_2014-11-11.csv")
entrezId2Name<-function(x){
  sapply(x, function(id){
    row <- which(probeInfo$entrez_id == id)
    probeInfo[row, 4]
  })
}

#Differential expression for each braak stage compare to non-braak stage
diffGenesPerBraakStage <- lapply(braakNames, function(bs){
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
      pvalUp <- testUp[["p.value"]]
      pvalDown <- testDown[["p.value"]]
      c(pval_up = pvalUp, pval_down = pvalDown)
    })))
    corr_pvalUp <- p.adjust(pvalues$pval_up, method = "BH", n = nGenes)
    corr_pvalDown <- p.adjust(pvalues$pval_down, method = "BH", n = nGenes)
    tab <- cbind(corr_pvalUp, corr_pvalDown)
    rownames(tab) <- rownames(pvalues)
    tab
  })
  # Convert list to two tables with p-values in all donors
  diffUp <- sapply(diffGenesPerBrain, function(t){
    t[ , "corr_pvalUp"]
  })
  diffDown <- sapply(diffGenesPerBrain, function(t){
    t[ , "corr_pvalDown"]
  })
  # Select genes signficant in all donors
  sigDiffUp <-apply(diffUp < 0.025, 1, sum) == nDonors
  sigDiffDown <- apply(diffDown < 0.025, 1, sum) == nDonors
  # Reduce diffUp and diffDown table to only significant genes
  diffUp <- diffUp[sigDiffUp, ]
  diffDown <- diffDown[sigDiffDown, ]
  # Sort genes based on lowest p-value found across donors
  sortedUp <- names(sort(apply(diffUp, 1, min)))
  sortedDown <- names(sort(apply(diffDown, 1, min)))
  diffUp <- diffUp[sortedUp, ]
  diffDown <- diffDown[sortedDown, ]
  #Format and add genes in table
  diffUp <- format(diffUp, digits = 4)
  diffDown <- format(diffDown, digits = 4)
  geneIdUp <- rownames(diffUp)
  geneIdDown <- rownames(diffDown)
  diffUp <- cbind(entrez_id = geneIdUp, gene_symbol = entrezId2Name(geneIdUp), diffUp)
  diffDown <- cbind(entrez_id = geneIdDown, gene_symbol = entrezId2Name(geneIdDown), diffDown)
  #save tables
  write.table(diffUp, file = paste("DiffExpr_braakStages/diffUp_", "bs", ".txt", sep = ""), quote = FALSE, row.names = FALSE)
  write.table(diffDown, file = paste("DiffExpr_braakStages/diffDown_", "bs", ".txt", sep = ""), quote = FALSE, row.names = FALSE)
  # Output
  list(upregulated = diffUp, downregulated = diffDown)
})
save(diffGenesPerBraakStage, file = "resources/diffGenesPerBraakStage.RData")
#############################################
# Examine differential expressed genes

load("resources/diffGenesPerBraakStage.RData")
is.present.up <- function(x){
  res <- sapply(diffGenesPerBraakStage, function(bs){
    upGenes <- bs[["upregulated"]]
    as.numeric(x %in% upGenes[ , "gene_symbol"])
  })
  rownames(res) <- x
  res
}

is.present.down <- function(x){
  res <- sapply(diffGenesPerBraakStage, function(bs){
    downGenes <- bs[["downregulated"]]
    as.numeric(x %in% downGenes[ , "gene_symbol"])
  })
  rownames(res) <- x
  res
}

# Presence of high impact genes
hiGenes <- c("GBA", "LRRK2", "PINK1", "PARK7", "SCNA", "VPS35", "DNAJC13", "CHCHD2") # High impact genes
presence_hiGenes <- is.present(hiGenes)

#Genes differentially expressed in all regions
allUp <- lapply(diffGenesPerBraakStage, function(bs){bs[["upregulated"]][ , "gene_symbol"]})
allDown <- lapply(diffGenesPerBraakStage, function(bs){bs[["downregulated"]][ , "gene_symbol"]})
allStagesUp <- Reduce(intersect, allUp)
allStagesDown <- Reduce(intersect, allDown)
allStagesUp
allStagesDown

# Comparing top 10 genes for each region to other regions
top10up <- lapply(diffGenesPerBraakStage, function(bs){
  bs[["upregulated"]][1:10, "gene_symbol"]
})
top10down <- lapply(diffGenesPerBraakStage, function(bs){
  bs[["downregulated"]][1:10, "gene_symbol"]
})
lapply(top10up, is.present.up)
lapply(top10down, is.present.down)