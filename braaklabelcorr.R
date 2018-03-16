#Correlate expression  with braakstage labels

setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(metafor)

load("resources/braakLabels.RData")

braakSize <- sapply(braakLabels, function(x){sum(sapply(c(1:6), function(b){sum(x==b)}))})

### FUNCTIONS ###

# correlation with Braak labels for all genes and donors
braak.cor <- function(dataList, stat) {
  as.data.frame(sapply(donorNames, function(d){
    labels <- braakLabels[[d]] # of all samples
    labels <- as.numeric(labels)
    labels <- labels[labels != 0]
    expr <- dataList[[d]]
    apply(expr, 1, function(g){# corr. each gene
      g <- unlist(g)
      r <- cor.test(g, labels) # Pearson's r by default
      r[[stat]]
    })
  }))
}

# Transform correlations to Fisher's z-scale and get corresponding sampling variances and confidence intervals for each row (gene)
z.transform <- function(mat){
  apply(mat, 1, function(r){
    t <- escalc(measure = "ZCOR", ri = r, ni = braakSize)
    t <- summary(t)
    t <- t[, c(1,2,5,6)]
    colnames(t) <- c("r", "variance", "lower95", "upper95")
    t
  })
}

# Correlated expression with labels and get summary effect
summary.braak.cor <- function(matList){
  
  genes <- rownames(matList[[1]])
  
  # Get correlations, transformed correlation and pvalues
  geneLabelCor <- braak.cor(matList, "estimate")
  geneLabelZscore <- z.transform(geneLabelCor)
  geneLabelPVal <- braak.cor(matList, "p.value")
  
  # Summary effect of correlations
  sapply(genes, function(g){
    tab <- geneLabelZscore[[g]]
    pvalue <- unlist(geneLabelPVal[g,])
    donors <- sapply(rownames(tab), function(n){ gsub("donor", "Donor ", n)})
    summary <- rma(tab$r, tab$variance, method = "DL", test = "t")
    weight <- round(weights(summary), digits = 2)
    tab <- cbind(donors, tab, pvalue, braakSize, weight)
    rbind(tab, 'summary' = list("Summary", summary$beta, summary$se^2 , summary$ci.lb, summary$ci.ub, 
                                summary$pval, sum(braakSize), sum(weight)))
  }, simplify = FALSE)
  
}

###########

# For all genes

# Get correlations, transformed correlation and pvalues
samples <- lapply(braakLabels, function(x) x != "0")
expr <- select.expr(samples = samples) #AHBA expression
summaryCor <- summary.braak.cor(expr)

# Back transform correlations, variances and confidence intervals
summaryLabelCorr <- lapply(summaryCor, function(t){
  t[, c("r", "variance", "lower95", "upper95")] <- back.transform(t[, c("r", "variance", "lower95", "upper95")])
  t
})

save(summaryLabelCorr, file = "resources/summaryLabelCorr.RData")


# same for eigen genes

load("resources/eigenExpr.RData") # expression of module eigen genes
eigenExpr <- lapply(eigenExpr, function(x)x[["braak1-6"]])
summaryCorEG <- summary.braak.cor(eigenExpr)

# Back transform correlations, variances and confidence intervals
summaryLabelCorrEG <- lapply(summaryCorEG, function(t){
  t[, c("r", "variance", "lower95", "upper95")] <- back.transform(t[, c("r", "variance", "lower95", "upper95")])
  t
})

save(summaryLabelCorrEG, file = "resources/summaryLabelCorrEG.RData")