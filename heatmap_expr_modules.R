# Heatmap Expression of modules
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")

load("../ABA_Rdata/BrainExprNorm.RData")
load("resources/modules_braak_absCor.RData")
modules_braak <- lapply(modules_braak, unlist)
load("resources/braakStages.RData")
braakNames <- braakNames[-c(2:5)]
load("resources/hierclust_tree_absCor.RData")

# Plotting method 2 using gplot2 and colors to indicate modules and braak regions
d= donorNames[1]
b= braakNames[1]

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

load("resources/braakLabels.RData")

pdf("heatmap_expr_modules_absCor.pdf", 8, 4)
# layout(matrix(c(1:12), 2, 6), weights = rep(2, 6), heights = rep(1, 2))

# lapply(braakNames, function(b){
lapply(c("braak1-6"), function(b){
  genes <- modules_braak[[b]]
  tree <- hierclust_tree[[b]][["average"]]
  
  gene_idx <- match(genes, tree$labels)
  modules <- tree$module[gene_idx]
  names(modules) <- genes
  colors <- tree$color[gene_idx]
  
  lapply(donorNames, function(d){
    
    # Subselect expression matrices
    samples <- braakStages[[d]]
    # samples <- as.logical(bitwOr(samples[, "braak1"], samples[, "braak6"]))
    # samples <- as.logical(samples[, b])
    samples <- as.logical(apply(samples, 1, sum))
    df <- sampleInfo[[d]][samples,]
    labels <- braakLabels[[d]][samples]
    exprMat <- as.matrix(brainExprNorm[[d]][genes, samples])
    
    colOrder <- order(labels, -df$graph_order)
    df <- df[colOrder, ]
    exprMat <- exprMat[, colOrder]
    
    colPal <- c("darkblue", "white", "darkred")
    rampcols <- colorRampPalette(colors = colPal, space="Lab")(200)
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
