# Gephi plot of Braak-related genes

setwd("/tudelft.net/staff-bulk/ewi/insy/DBL/Arlin/pd_braak")
# setwd("C:/Users/dkeo/surfdrive/pd_braak")
library("rgexf")
library(RColorBrewer)

source("PD/base_script.R")
load("resources/avgCor.RData")
load("resources/modules.RData")
load("resources/summaryLabelCorr.RData")

# modules <- lapply(modules, function(m){ m[sample(m$id, 200), ] })

# Simplify label correlations to summary estimate only
labelCor <- sapply(summaryLabelCorr, function(t) unlist(t["summary", "r"]) )

# Write gephi files for pos. and neg. Braak-related genes
lapply(names(modules), function(r){
  m <- modules[[r]]
  
  # Nodes
  genes <- m$id
  nodes <- m[,c("id", "label")]
  cor <- rep(substr(r, 1, 3), nrow(m))
  nodeSize <- abs(labelCor[genes]) * 100
  nodesAtt <- data.frame(membership = m$membership, cor = cor)
  nodeCol <- as.data.frame(t(col2rgb(m$color)))
  nodeCol$alpha <- rep(1, nrow(nodeCol))
  nodesVizAtt <- list(color = nodeCol, size = nodeSize)
  
  # Edges
  coexpr <- avgCor[genes, genes]
  # lower.tri(coexpr) <- 0 # rm duplicate links
  diag(coexpr) <- 0
  edges <- as.data.frame(as.table(coexpr))
  colnames(edges) <- c("source", "target", "coexpr")
  edges <- edges[edges$coexpr > 0.5, ] # filter edges 
  # print(min(edges$coexpr))
  colPal <- rev(heat.colors(100))
  edgeCol <- colPal[as.numeric(cut(edges$coexpr, breaks = seq(0,1, by =0.01)))]
  edgeCol <- as.data.frame(t(col2rgb(edgeCol)))
  edgeCol$alpha <- rep(1, nrow(edgeCol))
  edgesVizAtt <- list(color = edgeCol)
  
  
  write.gexf(nodes = nodes, edges = edges[, c("source", "target")], edgesWeight = edges$coexpr, 
             edgesVizAtt = edgesVizAtt,
             # edgesAtt = edgesAtt,
             nodesVizAtt = nodesVizAtt,
             nodesAtt = nodesAtt,
             defaultedgetype = "undirected", output = paste0("coexpr_", r, ".gexf"))
})