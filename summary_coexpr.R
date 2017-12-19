#Summary co-expression for all gene pairs in whole Braak region
setwd("/tudelft.net/staff-bulk/ewi/insy/DBL/Arlin/pd_braak")

library(metafor)

source("PD/base_script.R")
# load("resources/gene_coexpr.Rdata")
load("resources/gene_coexpr_subset.RData") # Subset genes
load("resources/braakLabels.RData") # Braak stage label vectors
back.transform <- dget("PD/back.transform.R")

#Correlations across donors per gene pair
llr <- lapply(gene_coexpr, function(x) x$r)
llp <- lapply(gene_coexpr, function(x) x$P)

genes <- rownames(gene_coexpr$donor9861$r)
geneNames <- entrezId2Name(genes)
geneNamePairs <- as.data.frame(t(combn(geneNames, 2)))
colnames(geneNamePairs) <- c("gene_A", "gene_B")
rownames(geneNamePairs) <- apply(geneNamePairs, 1, function(x) paste0(x, collapse = "-"))
genepairs <- as.data.frame(t(combn(genes, 2)))
colnames(genepairs) <- c("gene_A", "gene_B")
rownames(genepairs) <- rownames(geneNamePairs)

# Values for each genepair
genepair.mat <- function(ll, gp) {
  arr <- apply(simplify2array(ll), 1:2, c) # 3D data array genes*genes*donors
  t(apply(gp, 1, function(gene){
    A <- gene[1]
    B <- gene[2]
    arr[ , A, B]
  }))
}

genepairCor <- genepair.mat(llr, genepairs)
genepairPval <- genepair.mat(llp, genepairs)

# Sampling size
braakSize <- sapply(braakLabels, function(d){
  sum(d != 0)
})

# Transform correlations to Fisher's z-scale and get corresponding sampling variances and confidence intervals
genepairZscore <- apply(genepairCor, 1, function(r){
  t <- escalc(measure = "ZCOR", ri = r, ni = braakSize)
  t <- summary(t)
  t <- t[, c(1,2,5,6)]
  colnames(t) <- c("r", "variance", "lower95", "upper95")
  t
})

# Summary correlation (and uncorrected p-values)
summaryCor <- lapply(rownames(genepairs), function(gp){
  t <- genepairZscore[[gp]]
  t$pvalue <- genepairPval[gp, ]
  donors <- sapply(rownames(t), function(n){ gsub("donor", "Donor ", n)})
  summary <- rma(t$r, t$variance, method = "DL", test = "t")
  weight <- round(weights(summary), digits = 2)
  t <- cbind(donors, t, braakSize, weight)
  rbind(t, 'summary' = list("Summary", summary$beta, summary$se^2 , summary$ci.lb, summary$ci.ub, 
                            summary$pval, sum(braakSize), sum(weight)))
})

summaryGeneCoexpr <- lapply(summaryCor, function(t){
  t[, c("r", "variance", "lower95", "upper95")] <- back.transform(t[, c("r", "variance", "lower95", "upper95")])
  t
})

save(summaryGeneCoexpr, file = "resources/summaryGeneCoexpr_subset.RData")

#Print memory usage
print(object.size(x=lapply(ls(), get)), units="Mb")
