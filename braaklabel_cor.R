#Correlate expression  with braakstage labels
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")

library(metafor)
load("resources/braakInfo.RData")

brainExpr <- readRDS("../AHBA_Arlin/gene_expr.RDS")

##### Correlation on expression data #####
summary.braak.cor <- dget("PD/summary.braak.cor.R")

labels <- lapply(donorNames, function(d){
  idx <- unlist(braak_idx[[d]])
  braakLabels[[d]][idx]
})

expr_braak <- lapply(donorNames, function(d){
  idx <- unlist(braak_idx[[d]])
  brainExpr[[d]][, idx] # Expression in Braak regions
})
summaryLabelCor <- summary.braak.cor(expr_braak, labels)
save(summaryLabelCor, file = "resources/summaryLabelCor.RData")
