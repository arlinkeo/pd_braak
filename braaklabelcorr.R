#Correlate expression  with braakstage labels

options(stringsAsFactors = FALSE)
setwd("C:/Users/dkeo/surfdrive/Parkinson")

load("../polyQ_coexpression/resources/BrainExpr.RData")
donorNames <- names(brainExpr)
names(donorNames) <- donorNames
load("resources/braakLabels.RData")

# Correlate gene expression to Braak labels for each brain (without non-Braak region)
geneLabelCor <- sapply(donorNames, function(d){
  labels <- braakLabels[[d]]
  nbCols <- labels == 0
  labels <- labels[!nbCols]
  expr <- brainExpr[[d]][, !nbCols]
  apply(expr, 1, function(g){# corr. each gene
    g <- unlist(g)
    cor(g, as.numeric(labels)) # Pearson's r by default
  })
})
save(geneLabelCor, file = "resources/geneLabelCor.RData")