#Summary co-expression for all gene pairs in whole Braak region
setwd("/tudelft.net/staff-bulk/ewi/insy/DBL/Arlin/pd_braak")

library(metafor)

source("PD/base_script.R")
load("resources/gene_coexpr.Rdata")

# Genen pairs
genes <- rownames(braakExpr$donor9861)
genepairs <- as.data.frame(t(combn(genes, 2)))
rownames(genepairs) <- apply(genepairs, 1, function(x)paste(x, collapse = "-"))
colnames(genepairs) <- c("gene_A", "gene_B")

#Correlations across donors per gene pair
genepairCor <- t(apply(genepairs, 1, function(g){
  sapply(gene_coexpr, function(m){
    m$r[g[1],g[2]]
  })
}))

# Sampling size
braakSize <- sapply(wholeBraak, sum)

# Transform correlations to Fisher's z-scale and get corresponding sampling variances and confidence intervals
genepairZscore <- apply(genepairCor, 1, function(r){
  t <- escalc(measure = "ZCOR", ri = r, ni = braakSize)
  t <- summary(t)
  t <- t[, c(1,2,5,6)]
  colnames(t) <- c("r", "variance", "lower95", "upper95")
  t
})

# Summary correlation
summaryCor <- lapply(genepairZscore, function(t){
  donors <- sapply(rownames(t), function(n){ gsub("donor", "Donor ", n)})
  summary <- rma(t$r, t$variance, method = "DL")
  weight <- round(weights(summary), digits = 2)
  t <- cbind(donors, t, braakSize, weight)
  rbind(t, 'summary' = list("Summary", summary$beta, summary$se^2 , summary$ci.lb, summary$ci.ub, 
                            sum(braakSize), sum(weight)))
})

# Back transform correlations, variances and confidence intervals
back.transform <- function(x){
  a <- exp(1)^(2*x)
  as.data.frame((a-1)/(a+1))
}

summaryGeneCoexpr <- lapply(summaryCor, function(t){
  t[, c("r", "variance", "lower95", "upper95")] <- back.transform(t[, c("r", "variance", "lower95", "upper95")])
  t
})

save(summaryGeneCoexpr, file = "resources/summaryGeneCoexpr.RData")

#Print memory usage
print(object.size(x=lapply(ls(), get)), units="Mb")
