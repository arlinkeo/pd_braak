# Heatmap Expression of modules
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")

load("../ABA_Rdata/BrainExprNorm.RData")
load("resources/modules_braak_absCor.RData")
modules_braak <- lapply(modules_braak, unlist)
load("resources/braakStages.RData")
braakNames <- braakNames[-c(2:5)]
load("resources/hierclust_tree_absCor.RData")

# Plotting method 1 using ggplot2 and facets to indicate modules and braak regions
# library(ggplot2)
# library(reshape2)
# 
# lapply(braakNames, function(b){
#   genes <- modules_braak[[b]]
#   tree <- hierclust_tree[[b]][["average"]]
#   
#   gene_idx <- match(genes, tree$labels)
#   modules <- tree$module[gene_idx]
#   names(modules) <- genes
#   
#   lapply(donorNames, function(d){
#     
#     # Subselect expression matrices
#     labels <- braakLabels[[d]]
#     samples <- names(labels)[labels == 1 | labels == 6]
#     exprMat <- brainExprNorm[[d]][genes, samples]
#     
#     # Plot heatmap
#     mat <- as.data.frame(t(exprMat))
#     labels <- labels[samples]
#     tab <- cbind(sample = rownames(mat), braakstage = labels, mat)
#     tab.m <- melt(tab, id.vars = c("sample", "braakstage"), variable.name = "gene", value.name = "expr")
#     tab.m$braakstage <- factor(tab.m$braakstage, levels = c(1:6)) # keep order
#     tab.m$modules <- factor(modules[tab.m$gene], levels = unique(modules))
#     
#     p <- ggplot(tab.m, aes(sample, gene, group = braakstage)) +
#       geom_tile(aes(fill = expr)) +
#       facet_grid(modules~braakstage, scales = "free", space = "free", switch = "y") +
#       scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", name = "Expression") +
#       theme(axis.text = element_blank(),
#             axis.ticks = element_blank(),
#             panel.border = element_rect(fill = NA, colour = "black", size = 0.2)
#       ) +
#       labs(x = "Braak stage", y = "Module") +
#       ggtitle(paste("Expression of ", b," modules in Donor", tail(unlist(strsplit(d, split = "donor")), 1)))
#     print(p)
#     
#   })
# })

# Plotting method 2 using gplot2 and colors to indicate modules and braak regions
# d= donorNames[1]
# b= braakNames[1]

library(gplots)
ontology <- read.csv("../ABA_human_processed/Ontology_edited.csv")

# IDs and AHBA colors for each sample per donor
sampleInfo <- lapply(donorNames, function(d){
  sampleIds <- read.csv(paste("../ABA_human_processed/sample_info_normalized_microarray_", d, "_2014-11-11.csv", sep = ""))[ , 1]
  info <- ontology[match(sampleIds, ontology$id), ]
  info$color_hex_triplet <- sapply(info$color_hex_triplet, function(c){
    if (nchar(c) == 5) {paste("#0", c, sep = "")} else {paste("#", c, sep = "")}
  })
  info
})

pdf("heatmap_expr_modules_absCor.pdf", 8, 4)
# layout(matrix(c(1:12), 2, 6), weights = rep(2, 6), heights = rep(1, 2))

lapply(braakNames, function(b){
  genes <- modules_braak[[b]]
  tree <- hierclust_tree[[b]][["average"]]
  
  gene_idx <- match(genes, tree$labels)
  modules <- tree$module[gene_idx]
  names(modules) <- genes
  colors <- tree$color[gene_idx]
  
  lapply(donorNames, function(d){
    
    # Subselect expression matrices
    samples <- braakStages[[d]]
    samples <- as.logical(bitwOr(samples[, "braak1"], samples[, "braak6"]))
    df <- sampleInfo[[d]][samples,]
    exprMat <- as.matrix(brainExprNorm[[d]][genes, samples])
    
    colOrder <- rev(order(df$graph_order))
    df <- df[colOrder, ]
    exprMat <- exprMat[, colOrder]
    
    colPal <- c("darkblue", "white", "darkred")
    rampcols <- colorRampPalette(colors = colPal, space="Lab")(100)
    heatmap.2(exprMat, col = rampcols, 
              labRow = NA, labCol = df$acronym, 
              Rowv=FALSE, Colv=FALSE, 
              cexCol = 1, cexRow = 1,
              scale = "none", trace = "none", dendrogram = "none", key = FALSE, 
              RowSideColors = colors, ColSideColors = df$color_hex_triplet,
              main = paste0("Expression of modules in ", b, " and ", d),
              margins = c(8, 5))
              #, lwid = c(1, 10), lhei = c(1,10))
    
  })
})

dev.off()
