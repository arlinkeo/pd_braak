# select significant genes based on significant summary estimate

setwd("C:/Users/dkeo/surfdrive/Parkinson")

library("metap")

source("PD/base_script.R")
load("resources/summaryMeanDiff.RData")
load("resources/summaryLabelCorr.RData")

# Select non-Braak reference
summaryMeanDiff <- summaryMeanDiff$nonBraakA2

#Filter for summary effect
summaryMeanDiff <- lapply(summaryMeanDiff, function(g){
  as.data.frame(t(sapply(g, function(bs){
    bs["summary", ]
  })))
})

summaryCorr2 <- as.data.frame(t(sapply(summaryCorr, function(g){
  unlist(g["summary", ])
})))

# Get correlated genes |r>0.5|
corrGenes <- rownames(summaryCorr2)[abs(as.numeric(summaryCorr2$z)) > 0.5]

# Reduce table to only merged Braak regions 1-2, 3-4, 5-6
summaryMeanDiff <-lapply(summaryMeanDiff, function(gene){
  gene[braakNamesMerged1, ]
})

# Braak expression profiles  for all genes
braakProfile <- sapply(summaryMeanDiff, function(t) {
  diff <- unlist(t$meanDiff)
  pval <- unlist(t$benjamini_hochberg)
  profile <- as.numeric(abs(diff) > 1 & pval < 0.05)
  neg <- ifelse(diff <0, -1, 1)
  profile <- neg*profile
  paste(as.character(profile), collapse = "")
})
save(braakProfile, file = "resources/braakProfile.RData")

# Binary profiles in data
profiles <- sort(unique(braakProfile))
names(profiles) <- profiles

# Genes per profile
profileList <- lapply(profiles, function(p){
  idx <- braakProfile %in% p
  names(braakProfile)[idx]
})
sapply(profileList, length)
profileList <- profileList[names(profileList) != "000"]

# Diff. expressed genes in one of the merged braak regions for each non-Braak region
diffGenes <- unlist(profileList)

# Correlated genes with diff. expression in one of the braak regions
braakGenes <- intersect(corrGenes, diffGenes)
save(braakGenes, file = "resources/braakGenes.RData")
# any(braakGenes== "6622")

#Occurences of profiles for Braak-related genes (after selection)
braakProfile <- braakProfile[braakGenes]
occurences <- as.data.frame(table(braakProfile))

# Braak genes per profile
profiles <- profiles[profiles != c("000")]
profileList <- lapply(profiles, function(p){
  idx <- braakProfile %in% p
  names(braakProfile)[idx]
})
sapply(profileList, length)
save(profileList, file = "resources/profileList.RData")

#Save as txt-files
lapply(profiles, function(p){
  id <- profileList[[p]]
  name <- entrezId2Name(id)
  print(length(name))
  tab <-  data.frame(id, name)
  write.table(tab, file = paste0("Braak_related_genes/BraakGenes_", p, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
})

#Save as one txt-file
profileListNames <- lapply(profileList, entrezId2Name)
profileVector <- sapply(profileList, function(l) paste(l, collapse = ","))
profileVectorNames <- sapply(profileListNames, function(l) paste(l, collapse = ","))
tab <- data.frame(Profile = names(profileVector), Entrez_id = profileVector, Gene_name = profileVectorNames)
tab <- tab[lapply(tab$Gene_name, nchar) != 0, ]
write.table(tab, file = "Braak_related_genes/Braak_related_genes.txt", quote = FALSE, sep = "\t", row.names = FALSE)

#Presence of PD-implicated genes
profilePDgenes <- unlist(lapply(pdGenesID, function(pd){intersect(pd, braakGenes)}))
statsPD <- cbind(summaryCorr2[profilePDgenes,]$z, braakProfile[profilePDgenes])
rownames(statsPD) <- entrezId2Name(rownames(statsPD))
colnames(statsPD) <- c("r", "profile")
