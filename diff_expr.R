# Differential expression based on mean expression
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(metafor)
library(reshape2)
library(ggplot2)
library(plyr)
load("resources/braakInfo.RData")
load("../ABA_Rdata/BrainExpr.RData")
# expr1 <- readRDS("resources/expr_corrected_lm_mean_allsamples.rds")
# expr2 <- readRDS("resources/expr_corrected_lm_eg_allsamples.rds")

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
  expr <- expr2[[d]]
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
save(ttest, file = "resources/ttest_corrected_lm_eg_allsamples.RData")

# Print number of diff. expr. genes
apply(ttest, c(3,4), function(x){
  meanDiff <- x[, "meanA"] - x[, "meanB"]
  sum(x[, "BH"] < 0.05 & abs(meanDiff) > 1)
})

summaryDiffExpr <- sapply(dimnames(ttest)[[3]], function(r){
  print(r)
  sapply(dimnames(ttest)[[1]], function(gene){
    # print(paste0(r, ";", gene))
    gene <- ttest[gene,,r,]
    gene <- t(gene)
    
    # Get effect sizes and its variance (needed for meta-analysis) and confidence intervals (region B vs. A)
    t <- escalc(measure = "MD", 
                m1i = gene[, "meanB"], m2i = gene[, "meanA"], 
                n1i = gene[, "sizeB"], n2i = gene[, "sizeA"], 
                sd1i = sqrt(gene[, "varB"]), sd2i = sqrt(gene[, "varA"]))
    t <- summary(t)[, -c(3,4)]
    colnames(t) <- c("estimate", "var", "lower95", "upper95")
    
    # Summary effect size given effect sizes and variance
    summary <- rma(t$estimate, t$var, method = "DL", test = "t")

    # Combine into table
    t <- cbind(t, pvalue = gene[, "pvalue"], weight = weights(summary))
    rbind(t, 'summary' = list(summary$beta, summary$se^2 , summary$ci.lb, summary$ci.ub,
                              summary$pval, sum(weights(summary))))
  }, simplify = FALSE)
}, simplify = FALSE)
save(summaryDiffExpr, file = "resources/summaryDiffExpr_corrected_lm_eg_allsamples.RData")

########## Bar plot of differentially expressed genes between all Braak regions ##########

p <- plot.deg.numbers(summaryDiffExpr)
pdf("diff_expr_barplot_corrected_lm_eg_allsamples.pdf", 6, 4)
print(p)
dev.off()