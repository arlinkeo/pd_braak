#Correlate expression  with braakstage labels
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(metafor)

load("resources/braakInfo.RData")
labels <- lapply(braakLabels, function(x) x[x != 0] )# All Braak labels

summary.braak.cor <- dget("PD/summary.braak.cor.R")

# For all genes
samples <- lapply(braakLabels, function(x) x != 0)
expr <- readRDS("resources/expr_neuroncorrected.rds")
expr_braak <- lapply(donorNames, function(d){
  e <- expr[[d]]
  labels <- samples[[d]]
  e[, labels]
})
# expr <- select.expr(samples = samples) #AHBA expression
summaryLabelCor <- summary.braak.cor(expr_braak, labels)
save(summaryLabelCor, file = "resources/summaryLabelCor.RData")