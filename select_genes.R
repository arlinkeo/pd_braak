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

# Get top 10% correlated genes |r>0.6|
order <- rev(order(abs(labelCor$r))) # order absolute corr.
corrGenes <- rownames(labelCor)[order[1:1992]] # top 10% genes
sum(labelCor[corrGenes, "r"] < 0)
sum(labelCor[corrGenes, "r"] > 0)

# Top 10% significant Diff. expressed genes braak 1 vs. braak 6
diffExpr$benjamini_hochberg <- p.adjust(diffExpr$pvalue, method = "BH")
diffExpr$bonferroni <- p.adjust(diffExpr$pvalue, method = "bonferroni")
order <- order(diffExpr$benjamini_hochberg) # order absolute corr.
diffGenes1 <- rownames(diffExpr)[order[1:1992]] # top 10% genes
sum(diffExpr[diffGenes, "meanDiff"] < 0)
sum(diffExpr[diffGenes, "meanDiff"] > 0)

# Top 10% mean Diff. expressed genes braak 1 vs. braak 6
order <- rev(order(abs(diffExpr$meanDiff))) # order absolute corr.
diffGenes2 <- rownames(diffExpr)[order[1:1992]] # top 10% genes
sum(diffExpr[diffGenes2, "meanDiff"] < 0)
sum(diffExpr[diffGenes2, "meanDiff"] > 0)

# Correlated and diff. expressed
braakGenes <- Reduce(intersect, list(corrGenes, diffGenes1, diffGenes2))
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
   intersect(l, x)
  })
}

pd <- pd.genes(unlist(braakGenes))
lapply(pd, entrezId2Name)

# Print braak genes with correlations, mean diff, and p-values
gene.stats <- function(g){
  r <- labelCor[g, "r", drop = FALSE] 
  geneorder <- rownames(r)[order(r)]
  r <- r[geneorder, ]
  de <- diffExpr[geneorder, c("meanDiff", "benjamini_hochberg")]
  id <- geneorder
  name <- entrezId2Name(geneorder)
  df <- data.frame(id, name, r, de)
  df
}

tab <- Reduce(rbind, lapply(pd , gene.stats))
write.table(tab, file = "pdgenes_stats.txt", sep ="\t", quote = FALSE, row.names = FALSE)