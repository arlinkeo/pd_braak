# Summary mean expression in each Braak region across six brain donors

setwd("C:/Users/dkeo/surfdrive/pd_braak")
library(metafor)
source("PD/base_script.R")

load("../ABA_Rdata/BrainExprNorm.RData")
load("resources/braakStages.RData")

# Select data in samples
exprBraak <- lapply(donorNames, function(d){
  labels <- braakStages[[d]] == 1
  expr <- brainExprNorm[[d]]
  exprll <- apply(labels, 2, function(v){ # expr. in Braak regions 1-6
    expr[, v]
  })
})


# # Mean expression and variance in each Braak region and brain donor
# stat.expr <- function(fun) {
#   lapply(donorNames, function(d){
#     labels <- braakStages[[d]] == 1
#     expr <- brainExprNorm[[d]]
#     apply(labels, 2, function(v){ # expr. in Braak regions 1-6
#       m <- expr[, v]
#       apply(m, 1, fun)
#     })
#   })
# }
# 
# meanExpr <- stat.expr(mean)
# varExpr <- stat.expr(var)

# Info
genes <- rownames(brainExprNorm$donor9861)
nGenes <- length(genes)
sizesB <- t(sapply(braakStages, function(m)apply(m, 2, sum)))
# 
# # List of genes with mean and variance across donors for Braak regions 1-6
# geneTables <- lapply(braakNames, function(b){
#   sapply(genes, function(g){
#     mean <- sapply(donorNames, function(d){ meanExpr[[d]][g, b] })
#     var <- sapply(donorNames, function(d){ varExpr[[d]][g, b] })
#     data.frame(mean, var)
#   }, simplify = FALSE)
# })

# Summary mean expression
summaryMean <- lapply(braakNames[c(1, 6)], function(b){
  sapply(genes, function(g){
    print(paste0(b, ", ", g))
    
    # Mean, variance, CIs of expr 
    tab <- as.data.frame(t(sapply(donorNames, function(d){
      geneExpr <- unlist(exprBraak[[d]][[b]][g, ])
      t <- t.test(geneExpr)
      c(t$estimate, var(geneExpr), t$conf.int[1], t$conf.int[2], t$p.value)
    })))
    colnames(tab) <- c("mean", "var", "lower95", "upper95", "pvalue")
    tab$size <- sizesB[, b]
    
    # Summary statistics
    summary <- rma(yi = tab$mean, vi = tab$var, method = "DL", test = "t") # Summary effect size
    tab$weight <- weights(summary)
    summaryEffect <- c(summary$b, summary$se^2, summary$ci.lb, summary$ci.ub, summary$pval, sum(tab$size),sum(tab$weight))
    tab <- rbind(tab, 'summary' = summaryEffect)
  }, simplify = FALSE)
})
save(summaryMean, file = "resources/summaryMean.RData")

# Extract summary
summaryMean2 <- lapply(summaryMean, function(b){
  t(sapply(b, function(g){
    unlist(g["summary", ])
  }))
})

# Differential expression between summary means in Braak 1 and 6
lapply(genes, function(g){
  b1 <- summaryMean2$braak1[g, ]
  b6 <- summaryMean2$braak6[g, ]
  m1 <- b1$mean
  m6 <- b6$mean
  v1 <- b1$var
  v6 <- b6$var
  t <- (m1-m6)/sqrt(v1+v2)
})
