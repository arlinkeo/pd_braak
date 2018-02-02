# Differential expression MEDU (1) vs. FCTX (6)

setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")

# Load datasets
load("../UKBEC_RData/regionExpr.RData")
load("../UKBEC/expr.maps2.RData")

# All genes
probeNames <- rownames(regionExpr$CRBL)

# Gene probe mapping functions
geneProbeMap <- expr.map2[probeNames, ]
AHBA2exprID<- function(x){geneProbeMap$exprID[match(x, geneProbeMap$AHBA)]}
entrezId2exprID <- function(x) AHBA2exprID(entrezId2Name(x))
exprID2AHBA <- function(x){geneProbeMap$AHBA[match(x, geneProbeMap$exprID)]}

# T-test
ttest <- as.data.frame(t(sapply(probeNames, function(p){
  medu <- unlist(regionExpr$MEDU[p, ])
  fctx <- unlist(regionExpr$FCTX[p, ])
  t <- t.test(medu, fctx)
  meanDiff <- t$estimate[1] - t$estimate[2]
  names(meanDiff) <- NULL
  pvalue <- t$p.value
  c(meanDiff = meanDiff, pvalue = pvalue)
})))
ttest$benjamini_hochberg <- p.adjust(ttest$pvalue, method = "BH")
ttest$bonferroni <- p.adjust(ttest$pvalue, method = "bonferroni")

ttest$id <- probeNames
ttest$gene <- exprID2AHBA(probeNames)

order <- ttest$id[order(ttest$benjamini_hochberg)]
ttest <- ttest[order, ]

bgOrder <- intersect(order, bgAll)

# Diff. expr. braak genes
write.table(ttest[bgOrder,], file = "diffexpr_braakgenes_UKBEC.txt", quote = FALSE, row.names = FALSE)


ttest[intersect(order, entrezId2exprID(pdGenesID$susceptible)),]
ttest[entrezId2exprID(pdGenesID$lysosome),]

ttest[entrezId2exprID(pdGenesID$susceptible),]
