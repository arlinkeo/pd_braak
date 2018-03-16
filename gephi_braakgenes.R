# Gephi plot of modules with Braak-related genes

# setwd("/tudelft.net/staff-bulk/ewi/insy/DBL/Arlin/pd_braak")
source("PD/base_script.R")
setwd("C:/Users/dkeo/surfdrive/pd_braak")
library("rgexf")

# load("resources/braakGenes.RData")
# braakGenes <- unlist(braakGenes)

load("resources/modules.RData")
# modules_braak <- lapply(modules_braak, unlist)
# modules_braak <- lapply(modules_braak, function(m){ sample(m, 200) })

load("resources/hierclust_tree.RData")

# Summary label correlations
load("resources/summaryLabelCorr.RData")
labelCor <- sapply(summaryLabelCorr, function(t) unlist(t["summary", "r"]) )

# Write gephi files for modules in Braak 1 and 6
lapply(names(modules), function(r){  
  m <- hierclust_tree[[r]][["average"]]
  genes <- modules[[r]]
  gene_idx <- match(genes, m$labels)
  
  # Nodes
  nodes <- data.frame(id = genes, label = entrezId2Name(genes))
  cor <- labelCor[gene_idx]
  corSign <- ifelse(cor > 0, "+", "-")
  # nodeSize <- ifelse(genes %in% braakGenes, 2, 1)#abs(cor) * 100
  braak_gene <- ifelse(genes %in% braakGenes, 1, 0)
  nodesAtt <- data.frame(membership = m$module[gene_idx], cor, braak_gene, corSign)
  nodeCol <- data.frame(t(col2rgb(m$color[gene_idx])))
  nodeCol$alpha <- rep(1, nrow(nodeCol))
  nodesVizAtt <- list(color = nodeCol)#, size = nodeSize)
  
  # Edges
  file <- paste0("resources/avgCoexpr_", r, ".rds")
  coexpr <- readRDS(file)
  coexpr <- coexpr[genes, genes]
  coexpr[upper.tri(coexpr)] <- 0 # rm duplicate links
  diag(coexpr) <- 0
  edges <- as.data.frame(as.table(coexpr))
  colnames(edges) <- c("source", "target", "coexpr")
  edges <- edges[edges$coexpr != 0, ] # filter edges 
  edges <- edges[abs(edges$coexpr) > 0.9, ] #
  # filter edges 
  # print(min(edges$coexpr))
  edgeCol <- c("blue", "red")[as.numeric(edges$coexpr>0)+1]
  edgeCol <- as.data.frame(t(col2rgb(edgeCol)))
  edgeCol$alpha <- rep(1, nrow(edgeCol))
  edgesVizAtt <- list(color = edgeCol)
  edges$coexpr <- abs(edges$coexpr)
  
  write.gexf(nodes = nodes, edges = edges[, c("source", "target")], edgesWeight = edges$coexpr, 
             edgesVizAtt = edgesVizAtt,
             # edgesAtt = edgesAtt,
             nodesVizAtt = nodesVizAtt,
             nodesAtt = nodesAtt,
             defaultedgetype = "undirected", output = paste0("coexpr_", r, "_absCor.gexf"))

})
