# Heatmap co-expression of genes
setwd("C:/Users/dkeo/surfdrive/pd_braak")
library(gplots)
library(RColorBrewer)

source("PD/base_script.R")
load("resources/avgCor.RData")
load("resources/braakGenes.RData")

#Functions
cluster.genes <- function(x){
  distance <- dist(x)
  t <- hclust(distance)
  cutree(t, h = 4)
}

heatmap.coexpr <- function(n){
  g <- braakGenes[[n]]
  x <- avgCor[g, g]
  rownames(x) <- entrezId2Name(rownames(x))
  colnames(x) <- rownames(x)
  membership <- cluster.genes(x)
  colPal <- brewer.pal(length(unique(membership)), "Set1")
  # names(colPal) <- c(1:4)
  clustCol <- sapply(membership, function(x) colPal[x] )
  
  heatmap.2(x, col = rev(heat.colors(100)), 
            # col = colorRampPalette(c("blue", "white", "red"))(n = 100), 
            trace = "none", dendrogram = "row", labRow = NA, labCol = NA, key = TRUE, 
            breaks = seq(0,1,by=0.01), RowSideColors = clustCol)
  title(n)
}

# Plot heatmaps of pos. and neg. correlated Braak genes
pdf("heatmap_coexpr_braakgenes.pdf", 8, 8)
lapply(names(braakGenes), heatmap.coexpr)
dev.off()

# Modules with gene entrez IDs
modules <- lapply(names(braakGenes), function(r){
  g <- braakGenes[[r]]
  x <- avgCor[g, g]
  membership <- cluster.genes(x)
  membership <- sapply(membership, function(x) paste0(substr(r, 1, 3), x))
  data.frame(id = names(membership), membership = membership)
  # cNames <- unique(membership)
  # clusters <- lapply(cNames, function(c) {
  #   names(membership)[membership == c]
  # })
  # names(clusters) <- cNames
  # clusters
})
modules <- Reduce(rbind, modules)
save(modules, file = "resources/modules.RData")

#Print table

lapply(modules, function(r){lapply(r, entrezId2Name )})
