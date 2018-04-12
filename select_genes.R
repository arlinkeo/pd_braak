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

# number of genes in top 10
top10 <- floor(nrow(labelCor)*0.1)

# Get top 10% correlated genes
order <- rev(order(abs(labelCor$r))) # order absolute corr.
corrGenes <- rownames(labelCor)[order[1:top10]] # top 10% genes
corrGenes_neg <- corrGenes[labelCor[corrGenes, "r"] < 0]
corrGenes_pos <- corrGenes[labelCor[corrGenes, "r"] > 0]
min(abs(labelCor[corrGenes, "r"]))

# Top 10% significant Diff. expressed genes braak 1 vs. braak 6
diffExpr$benjamini_hochberg <- p.adjust(diffExpr$pvalue, method = "BH")
diffExpr$bonferroni <- p.adjust(diffExpr$pvalue, method = "bonferroni")
order <- order(diffExpr$benjamini_hochberg) # order absolute corr.
diffGenes1 <- rownames(diffExpr)[order[1:top10]] # top 10% genes
diffGenes1_up <- diffGenes1[diffExpr[diffGenes1, "meanDiff"] > 0]
diffGenes1_down <- diffGenes1[diffExpr[diffGenes1, "meanDiff"] < 0]
max(diffExpr[diffGenes1, "benjamini_hochberg"])

# Top 10% mean Diff. expressed genes braak 1 vs. braak 6
order <- rev(order(abs(diffExpr$meanDiff))) # order absolute corr.
diffGenes2 <- rownames(diffExpr)[order[1:top10]] # top 10% genes
diffGenes2_up <- diffGenes2[diffExpr[diffGenes2, "meanDiff"] > 0]
diffGenes2_down <- diffGenes2[diffExpr[diffGenes2, "meanDiff"] < 0]
min(abs(diffExpr[diffGenes2, "meanDiff"]))

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

# Print list of progression genes
progressionList <- unlist(braakGenes)
names <- entrezId2Name(progressionList)
tab <- cbind(entrez_id = progressionList, gene_symbol = names)
write.table(tab, file = "progression_genes.txt", row.names = FALSE, quote = FALSE, sep = "\t")

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
  r <- round(r, digits = 2)
  de <- diffExpr[geneorder, c("meanDiff", "benjamini_hochberg")]
  de$meanDiff <- round(de$meanDiff, digits = 2)
  de$benjamini_hochberg <- format(de$benjamini_hochberg, digits = 3, scientific = TRUE)
  id <- geneorder
  name <- entrezId2Name(geneorder)
  df <- data.frame(id, name, r, de)
  df
}

tab <- Reduce(rbind, lapply(pd , gene.stats))
write.table(tab, file = "pdgenes_stats.txt", sep ="\t", quote = FALSE, row.names = FALSE)
