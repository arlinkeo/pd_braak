# Function for summary Braak label correlation

# correlation with Braak labels for all genes and donors
braak.cor.data <- function(matList, labelList) {
  correlations <- sapply(names(matList), function(d){
    labels <- labelList[[d]] # of all samples
    expr <- matList[[d]]
    t(apply(expr, 1, function(gene){# corr. each gene/row
      gene <- unlist(gene)
      labels <- as.numeric(labels)
      r <- cor.test(gene, labels) # Pearson's r by default
      c(r = unname(r$estimate), pvalue = r$p.value, size = ncol(expr))
    })) # genes x measures
  }, simplify = FALSE)
  simplify2array(correlations) # genes x measures x donors
  # genell <- apply(correlations, 1, data.frame) # dataframe for each gene
  # lapply(genell, function(x) as.data.frame(t(x)))
}

# Transform correlations to Fisher's z-scale and get corresponding sampling variances and confidence intervals for each row (gene)
z.transform <- function(g){
  t <- escalc(measure = "ZCOR", ri = g$r, ni = g$size)
  t <- summary(t)
  t <- t[, c(1,2,5,6)]
  colnames(t) <- c("r", "variance", "lower95", "upper95")
  rownames(t) <- rownames(g)
  cbind(t, pvalue = g$pvalue, size = g$size)
}

# Back transform correlations
back.transform <- function(x){
  a <- exp(1)^(2*x)
  (a-1)/(a+1)
}

# Correlated expression with labels and get summary effect
summary.braak.cor <- function(matList, labelList){
  geneLabelCor <- braak.cor.data(matList, labelList) # get correlations and p-values
  geneLabelZscore <- lapply(geneLabelCor, z.transform) 
  
  # Summary effect of correlations
  tabList <- lapply(geneLabelZscore, function(t){
    summary <- rma(t$r, t$variance, method = "DL", test = "t")
    t$weight <- round(weights(summary), digits = 2)
    t$donors <- sapply(rownames(t), function(n){ gsub("donor", "Donor ", n)})
    t <- rbind(t, 'summary' = list(summary$beta, summary$se^2 , summary$ci.lb, summary$ci.ub, 
                                     summary$pval, sum(t$size), sum(t$weight), "Summary"))
    t[, c("r", "variance", "lower95", "upper95")] <- back.transform(t[, c("r", "variance", "lower95", "upper95")])
    t
  })
  simplify2array(tabList)
}