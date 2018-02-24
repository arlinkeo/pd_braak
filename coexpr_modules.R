# Co-expression modules in Braak 1 and 6

setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(WGCNA)

# Load co-expression matrices
braakCoexpr <- list(
  'braak1' = readRDS("resources/avgCoexpr_braak1.rds"), 
  'braak6' = readRDS("resources/avgCoexpr_braak6.rds")
)

### Hierarchical clustering and module detection ###
modules <- sapply(braakCoexpr, function(coexpr){
  # Distance matrix
  dist <- as.dist(1 - abs(coexpr))
  sapply(c("single", "complete", "average"), function(method){
    # Hierarchical clustering Braak 1
    tree <- hclust(dist, method = method)
    # Co-expression modules per Braak region
    tree$module <- cutreeDynamicTree(tree)
    table(tree$module) # count
    tree$color <- labels2colors(tree$module)
    tree;
  }, simplify = FALSE)
}, simplify = FALSE)
save(modules, file = "resources/modules.RData")

# Plot dendrogram and modules 
pdf("coexpr_modules.pdf", 12, 3)
lapply(names(modules), function(b) {
  treeList <- modules[[b]]
  lapply(names(treeList), function(method){
    title <- paste0("Co-expression modules in Braak ", b, ", ", method, " linkage")
    tree <- treeList[[method]]
    plotDendroAndColors(tree, colors = tree$color, dendroLabels = FALSE, main = title)
  })
})
dev.off()

# count modules
lapply(modules, function(b){
  lapply(b, function(t){
    table(t$module)
  })
})

# Modules with gene names
modules_genes <- lapply(modules, function(b){
  tree <- b$average
  modNames <- sort(unique(tree$module))[-1] # Unique module names, remove module 0 (gray)
  names(modNames) <- modNames
  groups <-  lapply(modNames, function(m){
    tree$labels[tree$module == m]
  })
})
save(modules_genes, file = "resources/module_genes.RData")