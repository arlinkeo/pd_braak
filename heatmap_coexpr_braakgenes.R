# Heatmap co-expression of genes
setwd("C:/Users/dkeo/surfdrive/pd_braak")
library(gplots)
library(RColorBrewer)

source("PD/base_script.R")
load("resources/avgCor.RData")
load("resources/braakGenes.RData")

# Get data frame with id, label, membership, and color
cluster.coexpr <- function(x){
  distance <- dist(x)
  t <- hclust(distance)
  membership <- cutree(t, h = 4)
  colPal <- brewer.pal(length(unique(membership)), "Set1")
  color <- sapply(membership, function(x) colPal[x] )
  id <- names(membership)
  label <- entrezId2Name(id)
  data.frame(id, label, membership, color)
}

heatmap.coexpr <- function(x, membership, title){
  heatmap.2(x, col = rev(heat.colors(100)), 
            trace = "none", dendrogram = "row", labRow = NA, labCol = NA, key = TRUE, 
            breaks = seq(0,1,by=0.01), RowSideColors = membership$color)
  title(title)
}

# Get coexpression and modules
coexpr <- lapply(braakGenes, function(g){  x <- avgCor[g, g] })
modules <- lapply(coexpr, cluster.coexpr)

# Plot heatmaps of pos. and neg. correlated Braak genes
pdf("heatmap_coexpr_braakgenes.pdf", 8, 8)
lapply(names(braakGenes), function(r) {
  heatmap.coexpr(coexpr[[r]], modules[[r]], r)
})
dev.off()

# Save data frame with module info
save(modules, file = "resources/modules.RData")