#Correlate expression  with braakstage labels

options(stringsAsFactors = FALSE)
setwd("C:/Users/dkeo/surfdrive/Parkinson")

load("../ABA_Rdata/BrainExpr.RData")
donorNames <- names(brainExpr)
names(donorNames) <- donorNames
load("../ABA_Rdata/BrainExprNorm.RData")
load("resources/braakLabels.RData")

braak.cor <- function(dataList) {
  sapply(donorNames, function(d){
    labels <- braakLabels[[d]]
    nbCols <- labels == 0
    labels <- labels[!nbCols]
    expr <- dataList[[d]][, !nbCols]
    apply(expr, 1, function(g){# corr. each gene
      g <- unlist(g)
      cor(g, as.numeric(labels)) # Pearson's r by default
    })
  })
}

# Correlate gene expression to Braak labels for each brain (without non-Braak region)
geneLabelCor <- braak.cor(brainExpr)
geneLabelCorNorm <- braak.cor(brainExprNorm)
save(geneLabelCor, file = "resources/geneLabelCor.RData")
save(geneLabelCorNorm, file = "resources/geneLabelCorNorm.RData")