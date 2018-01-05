# select significant genes based on significant summary estimate

setwd("C:/Users/dkeo/surfdrive/pd_braak")

library("metap")

source("PD/base_script.R")
load("resources/summaryDiffExpr.RData")
load("resources/summaryLabelCorr.RData")

# Select non-Braak reference
diffExpr <- summaryDiffExpr$`braak1-braak6`

#Filter for summary effect
diffExpr <- do.call(rbind.data.frame, lapply(diffExpr, function(g) g["summary",]))
labelCor <- do.call(rbind.data.frame, lapply(summaryLabelCorr, function(g) g["summary",]))

# Get correlated genes |r>0.5|
corrGenes <- rownames(labelCor)[abs(labelCor$r) > 0.5]

# Diff. expressed genes braak 1 vs. braak 6
diffGenes <- rownames(labelCor)[diffExpr$benjamini_hochberg < 0.05]

# Correlated and diff. expressed
braakGenes <- intersect(corrGenes, diffGenes)

# Sort by correlation
r <- labelCor[braakGenes, "r", drop = FALSE]
r <- r[order(r),]

save(braakGenes, file = "resources/braakGenes.RData")
# any(braakGenes== "6622")

#Presence of PD-implicated genes
