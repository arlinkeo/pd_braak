# Eigen gene differential expression

setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
load("../ABA_Rdata/BrainExprNorm.RData")
load("resources/modules.RData")
load("resources/braakStages.RData")

braakNames <- braakNames[-c(2:5)]

# PCA first component of subselection expr. matrices
eigenExpr <- lapply(donorNames, function(d){
  labels <- braakStages[[d]][, braakNames] == 1
  expr <- brainExprNorm[[d]]
  
  lapply(modules, function(m){ # For list of modules found in Braak 1 and 6 resp.
    as.data.frame(t(sapply(m, function(genes){ # For each module with genes (grouped gene rows)
      modExpr <- t(expr[genes, ]) # Expr of genes in modules across samples in one braak region
      pc <- prcomp(modExpr)#$x[, 1:2]# 1st PC (eigen gene expr)
      # plot(pc[,1], pc[,2])
      pc[,1]
    })))
  })
})
save(eigenExpr, file = "resources/eigenExpr.RData")

# Function T-test for each gene
ttestGene <- function(a, b) {
  test2tail <- t.test(a, b) # two-sided
  estimate <- test2tail$estimate
  names(estimate) <- NULL
  confidence95 <- test2tail$conf.int
  c('pvalue' = test2tail$p.value, 
    'meanA' = estimate[1], 'varA' = var(unlist(a)), 
    'meanB' = estimate[2], 'varB' = var(unlist(b)), 
    'lower95' = confidence95[1], 'upper95' = confidence95[2])
}

# T-test eigen expression in Braak 1 vs. 6
ttest <- lapply(donorNames, function(d){
  labels <- braakStages[[d]]
  modulesets <- eigenExpr[[d]]
  lapply(modulesets, function(expr){ # For list of modules found in Braak 1 and 6 resp.
    genesTab <- as.data.frame(t(apply(expr, 1, function(eg){ # For each eigen gene row
      region_a <- eg[as.logical(labels[, "braak1"])]
      region_b <- eg[as.logical(labels[, "braak6"])]
      ttestGene(region_a, region_b)
    })))
    genesTab$'benjamini_hochberg' <- p.adjust(genesTab$'pvalue' , method = "BH", n = nrow(genesTab)) #corrected p
    genesTab[order(genesTab$benjamini_hochberg),]
  })
})

# Saved datastructure for meta-analysis: module set -> Eigen genes -> table of Donors
diffExpr_eigengene <- lapply(braakNames, function(b){
  mods <- names(modules[[b]])
  list <- sapply(mods, function(eg){
    as.data.frame(t(sapply(donorNames, function(d){
      unlist(ttest[[d]][[b]][eg, ])
    })))
  }, simplify = FALSE)
})
save(diffExpr_eigengene, file = "resources/diffExpr_eigengene.RData")