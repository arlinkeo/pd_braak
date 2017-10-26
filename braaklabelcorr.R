#Correlate expression  with braakstage labels

setwd("C:/Users/dkeo/surfdrive/Parkinson")
source("PD/base_script.R")

load("../ABA_Rdata/BrainExprNorm.RData")
load("resources/braakLabels.RData")

braak.cor <- function(dataList) {
  sapply(donorNames, function(d){
    labels <- braakLabels[[d]]
    cols <- labels != 0
    labels <- as.numeric(labels[cols])
    expr <- dataList[[d]][, cols]
    apply(expr, 1, function(g){# corr. each gene
      g <- unlist(g)
      cor(g, labels) # Pearson's r by default
    })
  })
}

# Correlate gene expression to Braak labels for each brain (without non-Braak region)
geneLabelCor <- braak.cor(brainExprNorm)
save(geneLabelCor, file = "resources/geneLabelCor.RData")