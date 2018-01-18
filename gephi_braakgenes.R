# Gephi plot of Braak-related genes

setwd("C:/Users/dkeo/surfdrive/pd_braak")
library("rgexf")
library(RColorBrewer)

source("PD/base_script.R")
load("resources/avgCor.RData")
load("resources/modules.RData")

modules <- modules[sample(modules$id, 100), ]

# Nodes
genes <- modules$id
nodes <- data.frame(id = genes, label = entrezId2Name(genes))
colors <- col2hex(colors()[c(652, 630, 635, 640, 644, 610, 614, 27, 657, 204, 156, 525, 424, 386, 1)])
names(colors) <- unique(modules$membership)
memberColor <- sapply(modules$membership, function(x) colors[x])
nodesAtt <- data.frame(membership = modules$membership, color = memberColor, row.names = rownames(modules))

# Edges
coexpr <- avgCor[genes, genes]
edges <- as.data.frame(as.table(coexpr))
colnames(edges) <- c("source", "target", "coexpr")
edges <- edges[which(edges$coexpr >0.65), ]
# edges$weight.sign <- sapply(edges$coexpr, function(x) if (x>0) "+" else if (x<0) "-" else "0")
# edges$coexpr <- abs(edges$coexpr)

write.gexf(nodes = nodes, edges = edges[, c("source", "target")], edgesWeight = edges[, c("coexpr")], 
           # edgesVizAtt = edgesVizAtt, 
           # edgesAtt = as.data.frame(edges[, c("weight.sign")]), 
           # nodesVizAtt = nodeViz, 
           nodesAtt = nodesAtt,
           defaultedgetype = "undirected", output = "Images/coexpr_gephi/coexpr.gexf")
