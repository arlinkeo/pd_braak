# Differential expression based on mean expression
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(metafor)
library(reshape2)
library(ggplot2)
library(plyr)
load("resources/braakInfo.RData")
load("../ABA_Rdata/BrainExpr.RData")
# brainExpr <- readRDS("resources/expr_celltype_corrected.rds")

# Pairwise combinations of Braak regions 1-6
braakPairs <- t(combn(braakNames, 2))
rownames(braakPairs) <- apply(braakPairs, 1, paste, collapse = "-")
colnames(braakPairs) <- c("region_A", "region_B")
# 
# # Function T-test for each gene
# ttestGene <- function(a, b) {
#   test2tail <- t.test(a, b) # two-sided
#   estimate <- test2tail$estimate
#   names(estimate) <- NULL
#   confidence95 <- test2tail$conf.int
#   c('pvalue' = test2tail$p.value,
#     'meanA' = estimate[1], 'varA' = var(unlist(a)),
#     'meanB' = estimate[2], 'varB' = var(unlist(b)),
#     'sizeA' = length(a), 'sizeB' = length(b),
#     'lower95' = confidence95[1], 'upper95' = confidence95[2])
# }

# T-test to get p-values and CI's (only needed for forest plot of meta-analysis)
t.test.table <- dget("PD/t.test.table.R")

ttest <- lapply(donorNames, function(d){
  print(d)
  expr <- brainExpr[[d]]
  exprll <- lapply(braak_idx[[d]], function(b){ # expr. in Braak regions 1-6
    expr[, b]
  })
  tab <- alply(braakPairs, 1, function(r){
    print(r)
    region_a <- exprll[[r[1]]]
    region_b <- exprll[[r[2]]]
    t.test.table(region_a, region_b)
  }, .dims = TRUE) # keep names
  simplify2array(tab) # 3D array: genes x measures x region pairs
})
ttest <- simplify2array(ttest) # 4D array: genes x measures x region pairs x donors
# save(ttest, file = "resources/ttest.RData")

# Number of diff. expr. genes
apply(ttest, c(3,4), function(x){
  meanDiff <- x[, "meanA"] - x[, "meanB"]
  sum(x[, "BH"] < 0.05 & abs(meanDiff) > 1)
})

summaryDiffExpr <- apply(ttest, 3, function(r){
  apply(r, 1, function(gene){
    gene <- t(gene)
    
    # Get effect sizes and its variance (needed for meta-analysis) and confidence intervals (region B vs. A)
    t <- escalc(measure = "MD", 
                m1i = gene[, "meanB"], m2i = gene[, "meanA"], 
                n1i = gene[, "sizeB"], n2i = gene[, "sizeA"], 
                sd1i = sqrt(gene[, "varB"]), sd2i = sqrt(gene[, "varB"]))
    t <- summary(t)[, -c(3,4)]
    colnames(t) <- c("estimate", "var", "lower95", "upper95")
    
    # Summary effect size given effect sizes and variance
    summary <- rma(t$estimate, t$var, method = "DL", test = "t") 
    
    # Combine into table
    t <- cbind(t, pvalue = gene[, "pvalue"], weight = weights(summary))
    rbind(t, 'summary' = list(summary$beta, summary$se^2 , summary$ci.lb, summary$ci.ub,
                              summary$pval, sum(weights(summary))))
  })
})
# 
# # Braak region pairs -> Genes -> Donors (inverted table)
# diffExpr <- sapply(rownames(braakPairs), function(p){
#   sapply(genes, function(g){
#     as.data.frame(t(sapply(donorNames, function(d){
#       res <- unlist(ttest[[d]][[p]][g, ])
#     })))
#   }, simplify = FALSE)
# }, simplify = FALSE)
# save(diffExpr, file = "resources/diffExpr.RData")
# 
# # Meta-analysis of mean expression difference across donors
# summaryDiffExpr <- sapply(names(diffExpr), function(rp){ # For each Braak region pair
#   print(rp)
#   rp <- diffExpr[[rp]]
#   lapply(rp, function(gene){
# 
#     # Get variance and confidence intervals
#     m1i <- gene$meanA
#     m2i <- gene$meanB
#     n1i <- gene$sizeA
#     n2i <- gene$sizeB
#     sd1i <- sqrt(gene$varA)
#     sd2i <- sqrt(gene$varB)
#     # Get mean difference and its variance
#     # Flipped values so that +ve correlated genes have +ve fold-change, and vice versa
#     t <- escalc(measure = "MD", m1i = m2i, m2i = m1i, n1i = n2i, n2i = n1i, sd1i = sd2i, sd2i = sd1i)
#     t <- summary(t)
# 
#     # Get summary estimate
#     summary <- rma(t$yi, t$vi, method = "DL", test = "t") # Summary effect size
# 
#     donors <- sapply(rownames(gene), function(n){ gsub("donor", "Donor ", n)})
#     meanDiff <- as.numeric(t$yi)
#     varDiff <- t$vi
#     lower95 <- t$ci.lb
#     upper95 <- t$ci.ub
#     weight <- round(weights(summary), digits = 2)
#     pvalue <- gene$pvalue
#     t <- data.frame(donors, meanDiff, varDiff, lower95, upper95, pvalue, weight)
# 
#     # Combine into table
#     rbind(t, 'summary' = list("Summary", summary$beta, summary$se^2 , summary$ci.lb, summary$ci.ub,
#                                    summary$pval, sum(weight)))
#   })
# }, simplify = FALSE)
# save(summaryDiffExpr, file = "resources/summaryDiffExpr.RData")

########## Bar plot of differentially expressed genes between all Braak regions ##########

plot.deg.numbers <- dget("PD/plot.deg.numbers.R")
p <- plot.deg.numbers(summaryDiffExpr)
pdf("diff_expr_barplot_uncorrected.pdf", 6, 4)
print(p)
dev.off()
