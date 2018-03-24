# Differential expression MEDU (1) vs. FCTX (6)

setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")

load("resources/braakGenes.RData")
load("../UKBEC/regionExpr.RData")
load("../UKBEC/probe2gene.map.RData")

load("resources/summaryDiffExpr.RData")
diffExpr <- summaryDiffExpr$`braak1-braak6`
diffExpr <- do.call(rbind.data.frame, lapply(diffExpr, function(g) g["summary",]))
diffExpr$benjamini_hochberg <- p.adjust(diffExpr$pvalue, method = "BH")

# Gene probe mapping functions
AHBA2exprID<- function(x){probe2gene.map$exprID[match(x, probe2gene.map$AHBA)]}
entrezId2exprID <- function(x) AHBA2exprID(entrezId2Name(x))
exprID2AHBA <- function(x){probe2gene.map$AHBA[match(x, probe2gene.map$exprID)]}
exprID2entrezID <- function(x){name2EntrezId(exprID2AHBA(x))}

# T-test in UKBEC
ttest_ukbec <- as.data.frame(t(sapply(probe2gene.map$exprID, function(p){
  medu <- unlist(regionExpr$MEDU[p, ])
  fctx <- unlist(regionExpr$FCTX[p, ])
  t <- t.test(medu, fctx)
  meanDiff <- t$estimate[1] - t$estimate[2]
  names(meanDiff) <- NULL
  pvalue <- t$p.value
  c(meanDiff = meanDiff, pvalue = pvalue)
})))
ttest_ukbec$benjamini_hochberg <- p.adjust(ttest$pvalue, method = "BH")
save(ttest_ukbec, file = "resources/ttest_ukbec.RData")
# 
# # Differential expressed genes top 10% based on corrected p-value and mean difference
# n <- floor(nrow(ttest)*0.1)
# deg1 <- rownames(ttest)[order(ttest$benjamini_hochberg)[1:n]]
# deg2 <- rownames(ttest)[order(ttest$meanDiff)[1:n]]
# deg <- intersect(deg1, deg2)
# deg_up <- deg[ttest[deg, "meanDiff"] > 0]
# deg_down <- deg[ttest[deg, "meanDiff"] < 0]
# 
# # Overlap with Braak progression genes
# prog_up <- entrezId2exprID(braakGenes$negative_r)
# prog_up <- prog_up[!is.na(prog_up)]
# overlap_up <- intersect(prog_up, deg_up)
# 
# prog_down <- entrezId2exprID(braakGenes$positive_r)
# prog_down<- prog_down[!is.na(prog_down)]
# overlap_down <- intersect(prog_down, deg_down)
# 
# prog_genes <- entrezId2exprID(unlist(braakGenes))
# prog_genes <- prog_genes[!is.na(prog_genes)]
# overlap <- intersect(prog_genes, deg)
# 
# # Differential expression of progression genes
# tab_ukbec <- ttest[prog_genes, c("meanDiff", "benjamini_hochberg")]
# tab_ahba <- diffExpr[exprID2entrezID(prog_genes), c("meanDiff", "benjamini_hochberg")]
# tab <- cbind(ukbec = tab_ukbec, ahba = tab_ahba, gene = exprID2AHBA(prog_genes))
# tab <- tab[rev(order(abs(tab$ukbec.meanDiff))), ]
# 
# # Braak genes ordered by p-value
# bgOrder <- intersect(order, bgAll)
# write.table(ttest[bgOrder,], file = "diffexpr_braakgenes_UKBEC.txt", quote = FALSE, row.names = FALSE)
# 
# 
# ttest[intersect(order, entrezId2exprID(pdGenesID$susceptible)),]
# ttest[entrezId2exprID(pdGenesID$lysosome),]
# 
# ttest[entrezId2exprID(pdGenesID$susceptible),]
