# Summary Effect size (mean difference) of differential expressed (all) genes

setwd("C:/Users/dkeo/surfdrive/Parkinson")
library(metafor)
library("metap")

source("PD/base_script.R")
# load("../ABA_Rdata/BrainExprNorm.RData")
load("resources/diffExprBraak.RData")
load("resources/braakStages.RData")

#Combine table of braak and merged braak 
braakStages2 <- lapply(donorNames, function(d){
  braak <- braakStages[[d]]
  braakM1 <- braakMerged1[[d]]
  braakM2 <- braakMerged2[[d]]
  cbind(braak, braakM1, braakM2)
})
braakNames2 <- c(braakNames, braakNamesMerged1, braakNamesMerged2)

sizesB <- t(sapply(braakStages2, function(m)apply(m, 2, sum)))
sizesNB <- sizesNB <- sapply(nonBraak, function(l){sapply(l, sum)})

genes <- rownames(diffExprRef$nonBraakA1$donor9861$braak1)
names(genes) <- genes

summaryMeanDiff <- lapply(names(diffExprRef), function(ref){ # For each non-Braak reference
  dll <- diffExprRef[[ref]]
  lapply(genes, function(g){# for each gene
    print(paste0("Ref: ", ref, ", #", which(genes == g), ", gene id: ", g, ", gene name: ", entrezId2Name(g)))
    tabll <- lapply(braakNames2, function(bs){# result table for each braak stage
      tab <- as.data.frame(t(sapply(donorNames, function(d){ # row for each donor
        unlist(dll[[d]][[bs]][g,]) 
      })))
      meanDiff <- tab$meanB - tab$meanNB # mean difference for each donor
      size <- sizesB[, bs]
      varDiff <- (tab$varB / size) + (tab$varNB / sizesNB[ , ref]) # variance of mean difference for each donor
      summary <- rma(meanDiff, varDiff, method = "DL") # Summary effect size
      weight <- round(weights(summary), digits = 2)
      donors <- sapply(rownames(tab), function(n){ gsub("donor", "Donor ", n)})
      pvalue <- summary$pval
      benjamini_hochberg <- p.adjust(pvalue , method = "BH", n = length(genes))
      bonferroni <- p.adjust(pvalue , method = "bonferroni", n = length(genes))
      tab <- cbind(donors, meanDiff, varDiff, tab[, c("lower95", "upper95", "pvalue", "benjamini_hochberg", "bonferroni")], size, weight)
      tab <- rbind(tab, 'summary' = list("Summary", summary$beta, summary$se^2 , summary$ci.lb, summary$ci.ub, 
                                         pvalue, benjamini_hochberg, bonferroni, sum(size), sum(weight)))
      tab
    })
  })
})
names(summaryMeanDiff) <- names(nonBraak)

save(summaryMeanDiff, file = "resources/summaryMeanDiff.RData")
