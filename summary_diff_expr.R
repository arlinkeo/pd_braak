# Summary Effect size (mean difference) of differential expressed (all) genes

setwd("C:/Users/dkeo/surfdrive/Parkinson")
library(metafor)
library("metap")

source("PD/base_script.R")
load("resources/diffExpr.RData")
load("resources/braakStages.RData")

sizesB <- t(sapply(braakStages, function(m)apply(m, 2, sum)))

summaryDiffExpr <- sapply(names(diffExpr), function(rp){ # For each Braak region pair
  print(rp)
  r <- unlist(strsplit(rp, split ="-"))
  rA <- r[1]
  rB <- r[2]
  rp <- diffExpr[[rp]]
  nGenes <- length(rp)
  lapply(rp, function(gene){ # For each gene
    
    meanDiff <- gene$meanA - gene$meanB # mean difference for each donor
    # size <- sizesB[, bs]
    varDiff <- (gene$varA / sizesB[, rA]) + (gene$varB / sizesB[, rB]) # nonpooled variance of mean difference 
    summary <- rma(meanDiff, varDiff, method = "DL", test = "t") # Summary effect size
    weight <- round(weights(summary), digits = 2)
    donors <- sapply(rownames(gene), function(n){ gsub("donor", "Donor ", n)})
    pvalue <- summary$pval
    benjamini_hochberg <- p.adjust(pvalue , method = "BH", n = nGenes)
    tab <- cbind(donors, meanDiff, varDiff, gene[, c("lower95", "upper95", "benjamini_hochberg")], weight)
    tab <- rbind(tab, 'summary' = list("Summary", summary$beta, summary$se^2 , summary$ci.lb, summary$ci.ub, 
                                       benjamini_hochberg, sum(weight)))
    tab
  })
}, simplify = FALSE)

save(summaryDiffExpr, file = "resources/summaryDiffExpr.RData")