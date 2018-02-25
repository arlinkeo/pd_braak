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
      modExpr <- expr[genes, ] # Expr of genes in modules across all samples
      prcomp(modExpr)$rotation[, 1] # 1st PC (eigen gene expr)
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
diffExpr_eigengene <- lapply(donorNames, function(d){
  labels <- braakStages[[d]]
  modulesets <- eigenExpr[[d]]
  lapply(modulesets, function(expr){ # For list of modules found in Braak 1 and 6 resp.
    genesTab <- as.data.frame(t(apply(expr, 1, function(eg){ # For each eigen gene row
      region_a <- eg[as.logical(labels[, "braak1"])]
      region_b <- eg[as.logical(labels[, "braak6"])]
      ttestGene(region_a, region_b)
    })))
    genesTab$'benjamini_hochberg' <- p.adjust(genesTab$'pvalue' , method = "BH", n = nrow(genesTab)) #corrected p
    genesTab
  })
})
save(diffExpr_eigengene, file = "resources/diffExpr_eigengene.RData")