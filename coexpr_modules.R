# Co-expression modules in Braak 1 and 6

setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(WGCNA)

### Hierarchical clustering and module detection ###
regions <- c("braak1", "braak6", "braak1-6")
hierclust_tree <- sapply(regions, function(r){
  file <- paste0("resources/avgCoexpr_", r, ".rds")
  coexpr <- readRDS(file)
  dist <- as.dist(1 - coexpr)
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
save(hierclust_tree, file = "resources/hierclust_tree.RData")

# Plot dendrogram and modules 
pdf("coexpr_modules.pdf", 12, 2)
lapply(names(hierclust_tree), function(b) {
  treeList <- hierclust_tree[[b]]
  lapply(names(treeList), function(method){
    title <- paste0("Co-expression modules in Braak ", b, ", ", method, " linkage")
    tree <- treeList[[method]]
    plotDendroAndColors(tree, colors = tree$color, dendroLabels = FALSE, marAll = c(0, 4, 0.2, 0), main = title)
  })
})
dev.off()

# count modules and missed genes
lapply(hierclust_tree, function(b){
  t(sapply(b, function(t){
    n <- max(t$module)
    missed <- table(t$module)["0"]
    c(modules= n, missed_genes = missed)
  }))
})

# Modules with gene names
modules <- lapply(hierclust_tree, function(b){
  tree <- b$complete#b$average
  modNames <- sort(unique(tree$module)) # Unique module names, remove module 0 (gray)
  modNames <- modNames[modNames!="0"]
  names(modNames) <- modNames
  groups <-  lapply(modNames, function(m){
    tree$labels[tree$module == m]
  })
})
save(modules, file = "resources/modules.RData")

# overlap of genes between modules in B1 and B6
overlap <- sapply(modules$braak1, function(m1){
  sapply(modules$braak6, function(m6){
    overlap <- length(intersect(m1, m6))
    jaccard <- overlap/(length(m1) + length(m6) - overlap)
    round(jaccard*100, digits = 0)
  })
})
library(gplots)
rampcols <- colorRampPalette(colors = c("blue", "white", "red"), space="Lab")(200)
heatmap.2(overlap, col = rampcols, 
          Rowv=FALSE, Colv=FALSE, 
          cexCol = 1, cexRow = 1,
          scale = "none", trace = "none", dendrogram = "none", 
          main = "Jaccard index module overlap")
#, lwid = c(1, 10), lhei = c(1,10))
