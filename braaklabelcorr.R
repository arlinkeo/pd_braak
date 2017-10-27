#Correlate expression  with braakstage labels

setwd("C:/Users/dkeo/surfdrive/Parkinson")
source("PD/base_script.R")
library(metafor)

load("../ABA_Rdata/BrainExprNorm.RData")
load("resources/braakLabels.RData")

braak.cor <- function(dataList, stat) {
  sapply(donorNames, function(d){
    labels <- braakLabels[[d]]
    cols <- labels != 0
    labels <- as.numeric(labels[cols])
    expr <- dataList[[d]][, cols]
    apply(expr, 1, function(g){# corr. each gene
      g <- unlist(g)
      r <- cor.test(g, labels) # Pearson's r by default
      r[[stat]]
    })
  })
}

# Correlate gene expression to Braak labels for each brain (without non-Braak region)
geneLabelCor <- braak.cor(brainExprNorm, "estimate")
geneLabelPVal <- braak.cor(brainExprNorm, "p.value")


##### Summary effect of correlations
genes <- rownames(brainExprNorm$donor9861)
braakSize <- sapply(braakLabels, function(x){sum(sapply(c(1:6), function(b){sum(x==b)}))})

# Transform correlations to Fisher's z-scale and get corresponding sampling variances and confidence intervals
geneLabelZscore <- apply(geneLabelCor, 1, function(r){
  t <- escalc(measure = "ZCOR", ri = r, ni = braakSize)
  t <- summary(t)
  t <- t[, c(1,2,5,6)]
  colnames(t) <- c("z", "variance", "lower95", "upper95")
  t
})

# Get summary correlation
summaryCorr <- lapply(genes, function(g){
  tab <- geneLabelZscore[[g]]
  pvalue <- geneLabelPVal[g,]
  donors <- sapply(rownames(tab), function(n){ gsub("donor", "Donor ", n)})
  summary <- rma(tab$z, tab$variance, method = "DL")
  weight <- round(weights(summary), digits = 2)
  tab <- cbind(donors, tab, pvalue, braakSize, weight)
  rbind(tab, 'summary' = list("Summary", summary$beta, summary$se^2 , summary$ci.lb, summary$ci.ub, 
                              summary$pval, sum(braakSize), sum(weight)))
})
names(summaryCorr) <- genes
save(summaryCorr, file = "resources/summaryCorr.RData")