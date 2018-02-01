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
corrGenes <- rownames(labelCor)[abs(labelCor$r) > 0.65]
sum(labelCor[corrGenes, "r"] < 0)
sum(labelCor[corrGenes, "r"] > 0)

# Diff. expressed genes braak 1 vs. braak 6
diffExpr$benjamini_hochberg <- p.adjust(diffExpr$pvalue, method = "BH")
diffExpr$bonferroni <- p.adjust(diffExpr$pvalue, method = "bonferroni")
diffGenes <- rownames(diffExpr)[diffExpr$benjamini_hochberg < 0.00125 & abs(diffExpr$meanDiff) > 1.5]

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
   intersect(l, x)
  })
}

pd <- pd.genes(unlist(braakGenes))

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
