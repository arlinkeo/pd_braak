# select significant genes based on significant summary estimate

setwd("C:/Users/dkeo/surfdrive/Parkinson")

library("metap")

source("PD/base_script.R")
load("resources/summaryMeanDiff.RData")
load("resources/summaryCorr.RData")

#Filter for summary effect
summaryMeanDiff <- lapply(summaryMeanDiff, function(ref){
  lapply(ref, function(g){
    as.data.frame(t(sapply(g, function(bs){
      bs["summary", ]
    })))
  })
})
summaryCorr2 <- as.data.frame(t(sapply(summaryCorr, function(g){
  unlist(g["summary", ])
})))

# Get correlated genes |r>0.5|
corrGenes <- rownames(summaryCorr2)[abs(as.numeric(summaryCorr2$z)) > 0.5]

# save(corrGenes, file = "resources/correlated_genes1.RData")
# load("resources/correlated_genes1.RData")

# Diff. expressed genes in one of the merged braak regions for each non-Braak region
diffGenes <- lapply(summaryMeanDiff, function(ref){
  genes <- sapply(ref, function(g){
    g <- g[braakNamesMerged1, ] # check only in merged braak
    any(abs(unlist(g$meanDiff)) > 1 & g$benjamini_hochberg <0.05)
  })
  names(genes)[genes]
})

similarity <- sapply(diffGenes, function(a){ #jaccard
  sapply(diffGenes, function(b){
    shared <- length(intersect(a, b))
    unshared <- length(union(a, b))
    shared/unshared*100
  })
})#[c(1,3,5,2,4,6), c(1,3,5,2,4,6)]

overlapRef <- Reduce(intersect, diffGenes)
overlapRefA <- Reduce(intersect, diffGenes[c(1,3,5)])
overlapRefB <- Reduce(intersect, diffGenes[c(2,4,6)])

# Correlated genes with diff. expression in one of the braak regions
braakGenesPerRef <- lapply(diffGenes, function(ref){
  intersect(ref, corrGenes)
})
braakGenes <- braakGenesPerRef$nonBraakA2

save(braakGenes, file = "resources/braakGenes.RData")
any(braakGenes== "6622")

braakGenesList <- summaryMeanDiff$nonBraakA2[braakGenes]
#Check which Braak regions show diff. expression in 1-2, 3-4, and/or 5-6
braakGenesList2 <- lapply(braakGenesList, function(x)x[braakNamesMerged1, ])
braakProfile <- sapply(braakGenesList2, function(t) {
  diff <- unlist(t$meanDiff)
  pval <- unlist(t$benjamini_hochberg)
  profile <- as.numeric(abs(diff) > 1 & pval < 0.05)
  neg <- ifelse(diff <0, -1, 1)
  profile <- neg*profile
  paste(as.character(profile), collapse = "")
})
occurences <- as.data.frame(table(braakProfile))
profiles <- as.character(occurences$braakProfile)
names(profiles) <- profiles

profileBraakGenes <- lapply(profiles, function(p){
  idx <- braakProfile %in% p
  names(braakProfile)[idx]
})
save(profileBraakGenes, file = "resources/profileBraakGenes.RData")

#Save as txt-file
profiles <- names(profileBraakGenes)
names(profiles) <- profiles
lapply(profiles, function(p){
  id <- profileBraakGenes[[p]]
  print(length(id))
  name <- entrezId2Name(id)
  print(length(name))
  tab <-  data.frame(id, name)
  write.table(tab, file = paste0("Braak_related_genes/BraakGenes_", p, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
})
#find gene
find.gene <- function(g){
  lapply(profiles, function(p){
    entrezId2Name(profileBraakGenes[p]) %in% g
  })
}


# Correlation of Braak profile genes
profileCorr <- lapply(profileBraakGenes, function(p){
  summaryCorr2[p, ]
})

#Presence of PD-implicated genes
lapply(profileBraakGenes, function(p){
  overlap <- lapply(pdGenesID, function(pd){intersect(pd, p)})
  sapply(overlap, entrezId2Name)
})

