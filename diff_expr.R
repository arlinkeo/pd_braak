# Differential expression in braak regions of PD

setwd("C:/Users/dkeo/surfdrive/Parkinson")
library("metap")
library(plyr)

source("PD/base_script.R")
load("../ABA_Rdata/BrainExprNorm.RData")
load("resources/braakStages.RData")

# Gene info
genes <- rownames(brainExprNorm$donor9861)
nGenes <- length(genes)

# Pairwise combinations of Braak regions 1-6
braakPairs <- t(combn(braakNames, 2))
rownames(braakPairs) <- apply(braakPairs, 1, paste, collapse = "-")
colnames(braakPairs) <- c("region_A", "region_B")

# Function T-test for each gene
ttestGene <- function(a, b) {# Braak vs. non-braak
  test2tail <- t.test(a, b) # two-sided
  estimate <- test2tail$estimate
  names(estimate) <- NULL
  confidence95 <- test2tail$conf.int
  c('pvalue' = test2tail$p.value, 
    'meanA' = estimate[1], 'varA' = var(unlist(a)), 
    'meanB' = estimate[2], 'varB' = var(unlist(b)), 
    'lower95' = confidence95[1], 'upper95' = confidence95[2])
}

# Donors -> Braak region pairs -> genes (table)
ttest <- lapply(donorNames, function(d){
  labels <- braakStages[[d]] == 1
  expr <- brainExprNorm[[d]]
  exprll <- apply(labels, 2, function(v){ # expr. in Braak regions 1-6
    expr[, v]
  })
  tabll <- alply(braakPairs, 1, function(p){
    region_a <- exprll[[p[1]]]
    region_b <- exprll[[p[2]]]
    genesTab <- as.data.frame(t(sapply(genes, function(g){
      ttestGene(region_a[g, ], region_b[g, ])
    })))
    genesTab$'benjamini_hochberg' <- p.adjust(genesTab$'pvalue' , method = "BH", n = nGenes) #corrected p
    genesTab
  }, .dims = TRUE) # keep names
  tabll
})



# Braak region pairs -> Genes -> Donors (table)
diffExpr <- sapply(rownames(braakPairs), function(p){
  lapply(genes, function(g){
    as.data.frame(t(sapply(donorNames, function(d){
      res <- ttest[[d]][[p]][g, ]
    })))
  })
}, simplify = FALSE)

save(diffExpr, file = "resources/diffExpr.RData")