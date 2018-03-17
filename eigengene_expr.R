# Eigen gene differential expression
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")

load("../ABA_Rdata/BrainExprNorm.RData")
load("resources/modules.RData")
load("resources/braakLabels.RData")

regions <- c("braak1", "braak6", "braak1-6")

samples <- lapply(braakLabels, function(labels){
  braak1 <- labels == "1"
  braak6 <- labels == "6"
  braak1to6 <- labels != "0"
  list(braak1 = braak1, braak6 = braak6, 'braak1-6' = braak1to6)
})

# eigen gene function for data matrix (samples x genes)
eigen.gene <- function(x){
  eg <- prcomp(x)$x[, 1]# 1st PC (eigen gene expr)
  mean <- apply(x, 1, mean)
  if (cor(eg, mean) > 0) eg else -eg # flip sign of eigen gene based on the data
}

# Function for eigen gene expression for each module in the same region
# x: data, l: list of modules with genes, s: logical vector of samples (columns)
eigen.data <- function(x, l, s){ 
  df <- as.data.frame(t(sapply(l, function(genes){ # For each module with genes (grouped gene rows)
    eigen.gene(t(x[genes, s]))
  })))
  colnames(df) <- names(s)[s]
  df
}

# PCA first component of subselection expr. matrices
eigenExpr <- sapply(regions, function(r){ # For list of modules found in different brain regions  
  m <-modules[[r]] # modules with lists of genes
  lapply(donorNames, function(d){
    s <- samples[[d]][[r]] # logical vector
    expr <- brainExprNorm[[d]]
    eigen.data(expr, m, s)
  })
}, simplify = FALSE)
save(eigenExpr, file = "resources/eigenExpr.RData")