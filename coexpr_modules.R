# Co-expression modules in Braak 1 and 6

# avgCoexpr <- readRDS("output/avgCoexpr_wholeBraak.rds")

### Hierarchical clustering and module detection ###
hier.clust <- function(coexpr){
  dist <- as.dist(1 - coexpr)
  sapply(c("single", "complete", "average"), function(method){
    print(method)
    # Hierarchical clustering Braak 1
    tree <- hclust(dist, method = method)
    # Co-expression modules per Braak region
    tree$module <- cutreeDynamicTree(tree)
    tree$color <- labels2colors(tree$module)
    table(tree$module) # count
    tree;
  }, simplify = FALSE)
}
hierclust_tree <- hier.clust(avgCoexpr)
saveRDS(hierclust_tree, file = "output/hierclust_tree.rds")

# Plot dendrogram and modules 
pdf("output/coexpr_modules.pdf", 12, 2)
lapply(names(hierclust_tree), function(method){
  title <- ""#paste0("Co-expression modules in Braak 1-6, ", method, " linkage")
  tree <- hierclust_tree[[method]]
  plotDendroAndColors(tree, colors = tree$color, dendroLabels = FALSE, marAll = c(0, 4, 1, 0), main = title)
})
dev.off()

# count modules and missed genes
t(sapply(hierclust_tree, function(t){
  n <- max(t$module)
  missed <- table(t$module)["0"]
  c(modules= n, missed_genes = missed)
}))

# Modules with gene names
t <- hierclust_tree$average
modNames <- sort(unique(t$module)) # Unique module names, remove module 0 (gray)
modNames <- modNames[modNames!="0"]
names(modNames) <- paste0("M", modNames)
modules <- lapply(modNames, function(m){
  t$labels[t$module == m]
})
 