#Correlate expression  with braakstage labels
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(metafor)
load("resources/braakInfo.RData")
load("resources/summaryCoef.RData")

##### Correlation on fold-changes #####

# Convert lists: gene -> braak coef

fc_cor <- as.data.frame(t(sapply(ahba.genes(), function(g){
  estimates <- sapply(braakNames[-1], function(b){
    unlist(summaryCoef[[b]][[g]]["summary", "Estimate"])
  })
  fc <- c(braak1 = 0, estimates)
  r <- cor.test(fc, c(1:6))
  c(r = unname(r$estimate), pvalue = r$p.value)
  # plot(c(1:6), fc)
})))
fc_cor$BH <- p.adjust(fc_cor$r, method = "BH")
save(fc_cor, file = "resources/fc_cor.RData")

##### Correlation on expression data #####

labels <- lapply(braakLabels, function(x) x[x != 0] )# All Braak labels

summary.braak.cor <- dget("PD/summary.braak.cor.R")

# For all genes
samples <- lapply(braakLabels, function(x) x != 0)
expr <- readRDS("../ABA_Rdata/BrainExprNorm.rds")#"resources/expr_neuroncorrected.rds")

expr_braak <- lapply(donorNames, function(d){
  e <- expr[[d]]
  labels <- samples[[d]]
  e[, labels]
})
# expr <- select.expr(samples = samples) #AHBA expression
summaryLabelCor <- summary.braak.cor(expr_braak, labels)
save(summaryLabelCor, file = "resources/summaryLabelCor.RData")
