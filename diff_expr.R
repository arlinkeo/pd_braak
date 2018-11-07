# Differential expression based on mean expression
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(metafor)
library(reshape2)
library(ggplot2)
library(plyr)
load("resources/braakInfo.RData")
load("../ABA_Rdata/BrainExpr.RData")

# Load functions
source("PD/plot.deg.numbers.R")
source("PD/t.test.table.R")

# Pairwise combinations of Braak regions 1-6
braakPairs <- t(combn(braakNames, 2))
rownames(braakPairs) <- apply(braakPairs, 1, paste, collapse = "-")
colnames(braakPairs) <- c("region_A", "region_B")

# T-test to get p-values and CI's (only needed for forest plot of meta-analysis)
ttest <- lapply(donorNames, function(d){
  print(d)
  expr <- brainExpr[[d]]
  exprll <- lapply(braak_idx[[d]], function(b){ # expr. in Braak regions 1-6
    expr[, b]
  })
  tab <- alply(braakPairs, 1, function(r){
    print(unname(r))
    a <- exprll[[r[1]]] # matrix group 1
    b <- exprll[[r[2]]] # matrix group 2
    t.test.table(a,b) # Test for all genes
  }, .dims = TRUE) # keep names
  simplify2array(tab) # 3D array: genes x measures x region pairs
})
ttest <- simplify2array(ttest) # 4D array: genes x measures x region pairs x donors
save(ttest, file = "resources/ttest.RData")

# Print number of diff. expr. genes
apply(ttest, c(3,4), function(x){
  meanDiff <- x[, "meanA"] - x[, "meanB"]
  sum(x[, "BH"] < 0.05 & abs(meanDiff) > 1)
})

# Meta-analysis across donors
summaryDiffExpr <- aaply(ttest, c(1,3), function(g){ # For each Braak region pair and gene
  gene <- as.data.frame(t(g))
  t <- escalc(measure = "MD",  # Get estimates, variance (needed for meta-analysis) and confidence intervals (region B vs. A)
              m1i = gene[, "meanB"], m2i = gene[, "meanA"], # estimate
              n1i = gene[, "sizeB"], n2i = gene[, "sizeA"], 
              sd1i = sqrt(gene[, "varB"]), sd2i = sqrt(gene[, "varA"]))
  t <- summary(t)[, -c(3,4)]
  colnames(t) <- c("Estimate", "Var", "lower95", "upper95")
  summary <- rma(t$Estimate, t$Var, method = "DL", test = "t")
  gene$weight <- weights(summary)
  t <- cbind(t, pvalue = gene[, "pvalue"], weight = weights(summary))
  t <- rbind(t, 'summary' = list(summary$beta, summary$se^2 , summary$ci.lb, summary$ci.ub,
                            summary$pval, sum(weights(summary))))
  as.matrix(t)
}) # 4D-array: genes x regions x measures x donors
save(summaryDiffExpr, file = "resources/summaryDiffExpr.RData")

# Filter summary estimates, and correct P-values
summaryDiffExpr2 <- summaryDiffExpr[,,"summary",] # genes x region x measures
summaryDiffExpr2 <- aaply(summaryDiffExpr2, 2, function(t){ # P-value corrected for genes
  b <- p.adjust(t[, "pvalue"], method = "BH")
  cbind(t, BH = b)
}) # region x genes x measures

diffGenes <- apply(summaryDiffExpr2, 1, function(t){
  down <-  rownames(t)[which(t[, "Estimate"] < -1 & t[,"BH"] < 0.05)]
  up <- rownames(t)[which(t[, "Estimate"] > 1 & t[,"BH"] < 0.05)]
  list(down = down, up = up)
})
numbers <- t(sapply(diffGenes, function(x){
  deg <- sapply(x, length)
  c(deg, sum = sum(deg))
}))

########## Bar plot of differentially expressed genes between all Braak regions ##########

p <- plot.deg.numbers(numbers[, -3])
pdf("diff_expr_barplot.pdf", 6, 4)
print(p)
dev.off()