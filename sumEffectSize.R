# Summary Effect size (mean difference) of differential expressed (all) genes

setwd("C:/Users/dkeo/surfdrive/Parkinson")
library(metafor)

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
sizesNB <- sapply(nonBraak, function(l){sapply(l, sum)})

genes <- rownames(diffExprRef$nonBraakA1$donor9861$braak1)
names(genes) <- genes

summaryEffect <- lapply(names(diffExprRef), function(ref){ # For each non-Braak reference
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
      donors <- sapply(rownames(tab), function(n){ paste0("Donor ", unlist(strsplit(n, split = "donor"))[2])})
      tab <- cbind(donors, meanDiff, varDiff, tab[, c("lower95", "upper95", "pvalue")], size, weight)
      tab <- rbind(tab, 'summary' = list("Summary", summary$beta, summary$se^2 , summary$ci.lb, summary$ci.ub, 
                                         summary$pval, sum(size), sum(weight)))
      tab
    })
  })
})
names(summaryEffect) <- names(nonBraak)

save(summaryEffect, file = "resources/summaryEffect.RData")

##################################################################################
# genes <- rownames(diffGenesList$braak1$donor9861)
# names(genes) <- genes
# 
# Function to get info table and calculate summary effect for a gene
# sumTable <- function(g){
#   print(paste(g, entrezId2Name(g)))
#   lapply(braakNames, function(bs){
#     donorList <- diffGenesList[[bs]]
#     tab <- as.data.frame(t(sapply(donorList, function(t) unlist(t[g, ]))))# get gene row for each donor
#     tab$meanDiff <- unlist(tab$meanB) - unlist(tab$meanNB)# mean difference for each donor
#     donors <- sapply(rownames(tab), function(n){ paste0("Donor ", unlist(strsplit(n, split = "donor"))[2])})
#     tab$donor <- donors
#     tab$size <- sizes[, bs]
#     tab$varDiff <- (tab$varB / tab$size) + (tab$varNB / sizesNB) # variance of mean difference for each donor
#     tab$pvalue <- tab$benjamini_hochberg # Corrected p-value
#     tab <- tab[, c("lower95", "upper95", "meanDiff", "varDiff", "donor", "size", "pvalue")]
#     sumEffect <- rma(tab$meanDiff, tab$varDiff, method = "DL") # Summary effect size
#     tab$weight <- round(weights(sumEffect), digits = 2)
#     # Add row with summary effect statistics
#     tab <- rbind(tab, 'SummaryEff' = list(sumEffect$ci.lb, sumEffect$ci.ub, sumEffect$beta, sumEffect$se^2, 
#                                           "Summary", sum(tab$size), sumEffect$pval, sum(tab$weight)))
#     tab$isSum <- tab$donor == "Summary"
#     tab$braak <- bs
#     tab
#   })
# }
# 
# #Get values for each gene
# sumEffectSize <- lapply(genes, sumTable)
# save(sumEffectSize, file = "resources/sumEffectSize.RData")
##################################################################################