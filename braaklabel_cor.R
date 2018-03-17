#Correlate expression  with braakstage labels
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(metafor)

load("resources/braakLabels.RData")

### FUNCTIONS ###

# Correlate labels with genes
braak.cor <- function(x, l, stat){ #x: data, l: labels
  apply(x, 1, function(g){# corr. each gene/row
    g <- unlist(g)
    r <- cor.test(g, as.numeric(l)) # Pearson's r by default
    r[[stat]]
  }) # vector with stats
}

# correlation with Braak labels for all genes and donors
braak.cor.data <- function(matList, stat, labelList) {
  as.data.frame(sapply(names(matList), function(d){
    labels <- labelList[[d]] # of all samples
    expr <- matList[[d]]
    braak.cor(expr, labels, stat)
  }))
}

# Transform correlations to Fisher's z-scale and get corresponding sampling variances and confidence intervals for each row (gene)
z.transform <- function(mat, size){
  apply(mat, 1, function(r){
    t <- escalc(measure = "ZCOR", ri = r, ni = size)
    t <- summary(t)
    t <- t[, c(1,2,5,6)]
    colnames(t) <- c("r", "variance", "lower95", "upper95")
    t
  })
}

# Correlated expression with labels and get summary effect
summary.braak.cor <- function(matList, labelList){
  
  genes <- rownames(matList[[1]])
  size <- sapply(matList, ncol)
  
  # Get correlations, transformed correlation and pvalues
  geneLabelCor <- braak.cor.data(matList, "estimate", labelList)
  geneLabelZscore <- z.transform(geneLabelCor, size)
  geneLabelPVal <- braak.cor.data(matList, "p.value", labelList)
  
  # Summary effect of correlations
  sapply(genes, function(g){
    tab <- geneLabelZscore[[g]]
    pvalue <- unlist(geneLabelPVal[g,])
    donors <- sapply(rownames(tab), function(n){ gsub("donor", "Donor ", n)})
    summary <- rma(tab$r, tab$variance, method = "DL", test = "t")
    weight <- round(weights(summary), digits = 2)
    tab <- cbind(donors, tab, pvalue, size, weight)
    t <- rbind(tab, 'summary' = list("Summary", summary$beta, summary$se^2 , summary$ci.lb, summary$ci.ub, 
                                     summary$pval, sum(size), sum(weight)))
    
    t[, c("r", "variance", "lower95", "upper95")] <- back.transform(t[, c("r", "variance", "lower95", "upper95")])
    t
  }, simplify = FALSE)
  
}

###########

labels <- lapply(braakLabels, function(x) x[x != 0] )# All Braak labels

# For all genes
samples <- lapply(braakLabels, function(x) x != "0")
expr <- select.expr(samples = samples) #AHBA expression
summaryLabelCorr <- summary.braak.cor(expr, labels)
save(summaryLabelCorr, file = "resources/summaryLabelCorr.RData")

# same for eigen genes
load("resources/eigenExpr.RData") # expression of module eigen genes
eigenExpr <- eigenExpr[["braak1-6"]]
summaryLabelCorrEG <- summary.braak.cor(eigenExpr, labels)
save(summaryLabelCorrEG, file = "resources/summaryLabelCorrEG.RData")