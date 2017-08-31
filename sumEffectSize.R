# Summary Effect size (mean difference) of differential expressed (all) genes

setwd("C:/Users/dkeo/surfdrive/Parkinson")
library(metafor)

source("PD/base_script.R")
load("resources/diffGenesBraak.RData")

load("resources/braakStages.RData")
braakNames <- names(braakStages)
names(braakNames) <- braakNames
sizes <- sapply(braakStages, function(bs){
  sapply(bs, sum)
})
sizesNB <- sapply(nonBraak, sum)

genes <- rownames(diffGenesList$braak1$donor9861)
names(genes) <- genes

# Function to get info table and calculate summary effect for a gene
sumTable <- function(g){
  print(paste(g, entrezId2Name(g)))
  lapply(braakNames, function(bs){
    donorList <- diffGenesList[[bs]]
    tab <- as.data.frame(t(sapply(donorList, function(t) unlist(t[g, ]))))# get gene row for each donor
    tab$meanDiff <- unlist(tab$meanB) - unlist(tab$meanNB)# mean difference
    donors <- sapply(rownames(tab), function(n){ paste0("Donor ", unlist(strsplit(n, split = "donor"))[2])})
    tab$donor <- donors
    tab$size <- sizes[, bs]
    tab$varDiff <- (tab$varB / tab$size) + (tab$varNB / sizesNB) # variance of mean difference
    tab$pvalue <- tab$benjamini_hochberg
    tab <- tab[, c("lower95", "upper95", "meanDiff", "varDiff", "donor", "size", "pvalue")]
    sumEffect <- rma(tab$meanDiff, tab$varDiff, method = "DL") # Summary effect size
    tab$weight <- round(weights(sumEffect), digits = 2)
    tab <- rbind(tab, 'SummaryEff' = list(sumEffect$ci.lb, sumEffect$ci.ub, sumEffect$beta, sumEffect$se^2, 
                                          "Summary", sum(tab$size), sumEffect$pval, sum(tab$weight)))
    tab$isSum <- tab$donor == "Summary"
    tab$braak <- bs
    tab
  })
}

#Get values for each gene
sumEffectSize <- lapply(genes, sumTable)
save(sumEffectSize, file = "resources/sumEffectSize.RData")