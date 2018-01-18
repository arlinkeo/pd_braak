# select significant genes based on significant summary estimate

setwd("C:/Users/dkeo/surfdrive/pd_braak")

library("metap")

source("PD/base_script.R")
load("resources/summaryDiffExpr.RData")
load("resources/summaryLabelCorr.RData")

# Select region pair
diffExpr <- summaryDiffExpr$`braak1-braak6`

#Filter for summary effect
diffExpr <- do.call(rbind.data.frame, lapply(diffExpr, function(g) g["summary",]))
labelCor <- do.call(rbind.data.frame, lapply(summaryLabelCorr, function(g) g["summary",]))

# Get correlated genes |r>0.6|
corrGenes <- rownames(labelCor)[abs(labelCor$r) > 0.6]
sum(labelCor[corrGenes, "r"] < 0)
sum(labelCor[corrGenes, "r"] > 0)

# Diff. expressed genes braak 1 vs. braak 6
diffExpr$benjamini_hochberg <- p.adjust(diffExpr$pvalue, method = "BH")
# diffExpr$bonferroni <- p.adjust(diffExpr$pvalue, method = "bonferroni")
diffGenes <- rownames(diffExpr)[diffExpr$benjamini_hochberg < 0.001 & abs(diffExpr$meanDiff) > 1]

# Correlated and diff. expressed
braakGenes <- intersect(corrGenes, diffGenes)
positive_r <- braakGenes[labelCor[braakGenes, "r"] > 0]
negative_r <- braakGenes[labelCor[braakGenes, "r"] < 0]
braakGenes <- list(positive_r = positive_r, negative_r = negative_r)

# # Sort by correlation
# braakGenes <- lapply(braakGenes, function(ll){
#   r <- labelCor[ll, "r", drop = FALSE]
#   r <- r[order(r), ,drop = FALSE]
#   rownames(r)
# })


save(braakGenes, file = "resources/braakGenes.RData")

#Presence of PD-implicated genes
pd.genes <- function(x){
  lapply(pdGenesID, function(l){
   res <- intersect(l, x)
   m <- labelCor[res, ]
   rownames(m) <- entrezId2Name(rownames(m))
   m
  })
}

pd.genes(unlist(braakGenes))
