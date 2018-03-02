# Summary of eigen gene differential expression in BRaak 1 vs. 6

setwd("C:/Users/dkeo/surfdrive/pd_braak")
library(metafor)
library("metap")

source("PD/base_script.R")
load("resources/diffExpr_eigengene.RData")
load("resources/braakStages.RData")

sizesB <- t(sapply(braakStages, function(m)apply(m, 2, sum)))
# nGenes <- length((diffExpr_eigengene$`braak1-braak2`))

braakNames <- braakNames[-c(2:5)]

b= braakNames[1]

summaryDiffExpr <- sapply(braakNames, function(b){ # For each st of modules for Braak region 1 and 6
  print(b)
  rA <- "braak1"
  rB <- "braak6"
  egList <- diffExpr_eigengene[[b]]
  x=lapply(egList, function(eg){
    
    # Get variance and confidence intervals
    m1i <- eg$meanA
    m2i <- eg$meanB
    n1i <- sizesB[, rA]
    n2i <- sizesB[, rB]
    sd1i <- sqrt(eg$varA)
    sd2i <- sqrt(eg$varB)
    # Get mean difference and its variance
    t <- escalc(measure = "MD", m1i = m1i, m2i = m2i, n1i = n1i, n2i = n2i, sd1i = sd1i, sd2i = sd2i)
    t <- summary(t)
    
    # Get summary estimate
    summary <- rma(t$yi, t$vi, method = "DL", test = "t") # Summary effect size

    donors <- sapply(rownames(eg), function(n){ gsub("donor", "Donor ", n)})
    meanDiff <- as.numeric(t$yi)
    varDiff <- t$vi
    lower95 <- t$ci.lb
    upper95 <- t$ci.ub
    weight <- round(weights(summary), digits = 2)
    pvalue <- eg$pvalue
    t <- data.frame(donors, meanDiff, varDiff, lower95, upper95, pvalue, weight)
    
    # Combine into table
    t <- rbind(t, 'summary' = list("Summary", summary$beta, summary$se^2 , summary$ci.lb, summary$ci.ub, 
                                   summary$pval, sum(weight)))
    t
    
  })
  
}, simplify = FALSE)

lapply(summaryDiffExpr, function(tab){
  diffExpr <- do.call(rbind.data.frame, lapply(tab, function(g) g["summary",]))
  
  diffExpr$benjamini_hochberg <- p.adjust(diffExpr$pvalue, method = "BH")
  diffExpr[order(diffExpr$benjamini_hochberg), ]
})
