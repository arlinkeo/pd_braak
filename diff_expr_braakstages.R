# Differential expression in braak regions of PD

options(stringsAsFactors = FALSE)
setwd("C:/Users/dkeo/surfdrive/Parkinson")
library("metap")

load("../ABA_Rdata/BrainExpr.RData")
donorNames <- names(brainExpr)
names(donorNames) <- donorNames
nDonors <- length(donorNames)
load("resources/braakStages.RData")
braakNames <- names(braakStages)
names(braakNames) <- braakNames
probeInfo <- read.csv("../ABA_human_processed/probe_info_2014-11-11.csv")
entrezId2Name <- function (x) {probeInfo$gene_symbol[match(x, probeInfo$entrez_id)]} #Input is vector

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
      fChange <- log2(2^test2tail$estimate[1] / 2^test2tail$estimate[2]) # estimate is mean
      names(fChange) <- NULL
      c(pval_2tail = pval2tail, 'fold-change' = fChange)
    })))
    pvalues$pval_2tail <- p.adjust(pvalues$pval_2tail , method = "BH", n = nGenes)
    pvalues
  })
})

# convert lists to gene x donor table for each braak stage
dg.table <- function(x){
  lapply(braakNames, function(bs){
    brainList <- diffGenesList[[bs]]
    genes <- rownames(brainList[[1]])
    tab <- sapply(brainList, function(t){v <- t[ , x]})
    rownames(tab) <- genes
    tab
  })
}

dgTable2tail <- dg.table("pval_2tail")
dgTableUp <- dg.table("pval_up")
dgTableDown <- dg.table("pval_down")

save(diffGenesList, dgTable2tail, dgTableUp, dgTableDown, file = "resources/diffGenesBraak.RData")

# Diff. expressed in all donors
diff.all <- function(list, th) {
  braakTab <- sapply(list, function(tab){
    apply(diffUp < th, 1, sum) == nDonors
  })
  rownames(braakTab) <- rownames(list[[1]])
  braakTab
}

dgAll2tail <- diff.all(dgTable2tail, 0.05)
dgAllUp <- diff.all(dgTableUp, 0.025)
dgAllDown <- diff.all(dgTableDown, 0.025)

#Combine p-values
combine.p <- function(list){
  sapply(list, function(tab){
    apply(tab, 1, function(p){
      pchisq(-2*sum(log(p)), df = length(p)*2, lower.tail = FALSE)
    })
  })
}

fisherP2tail <- combine.p(dgTable2tail)
fisherPUp <- combine.p(dgTableUp)
fisherPDown <- combine.p(dgTableDown)

#Get list of sorted significant genes
signif.genes <- function(tab){
  lapply(braakNames, function(bs){
    (head(tab[, bs])) < 0.05
  })
}

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