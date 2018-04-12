# Summary Effect size (mean difference) of differential expressed (all) genes

setwd("C:/Users/dkeo/surfdrive/pd_braak")
library(metafor)
library("metap")

source("PD/base_script.R")
load("resources/diffExpr.RData")
load("resources/braakStages.RData")

sizesB <- t(sapply(braakStages, function(m)apply(m, 2, sum)))
nGenes <- length((diffExpr$`braak1-braak2`))

summaryDiffExpr <- sapply(names(diffExpr), function(rp){ # For each Braak region pair
  print(rp)
  r <- unlist(strsplit(rp, split ="-"))
  rA <- r[1]
  rB <- r[2]
  rp <- diffExpr[[rp]]
  lapply(rp, function(gene){
    
    # Get variance and confidence intervals
    m1i <- gene$meanA
    m2i <- gene$meanB
    n1i <- sizesB[, rA]
    n2i <- sizesB[, rB]
    sd1i <- sqrt(gene$varA)
    sd2i <- sqrt(gene$varB)
    # Get mean difference and its variance
    t <- escalc(measure = "MD", m1i = m1i, m2i = m2i, n1i = n1i, n2i = n2i, sd1i = sd1i, sd2i = sd2i)
    t <- summary(t)
    
    # Get summary estimate
    summary <- rma(t$yi, t$vi, method = "DL", test = "t") # Summary effect size
    
    # tscore <- summary$b/summary$se
    # pval <- 2*pt(-abs(tscore), df = 2)
    # correctedp <- p.adjust(pval, n = 19992)
    
    donors <- sapply(rownames(gene), function(n){ gsub("donor", "Donor ", n)})
    meanDiff <- as.numeric(t$yi)
    varDiff <- t$vi
    lower95 <- t$ci.lb
    upper95 <- t$ci.ub
    weight <- round(weights(summary), digits = 2)
    pvalue <- gene$pvalue
    t <- data.frame(donors, meanDiff, varDiff, lower95, upper95, pvalue, weight)
    
    # Combine into table
    t <- rbind(t, 'summary' = list("Summary", summary$beta, summary$se^2 , summary$ci.lb, summary$ci.ub, 
                                   summary$pval, sum(weight)))
    t
    
  })
}, simplify = FALSE)

save(summaryDiffExpr, file = "resources/summaryDiffExpr.RData")