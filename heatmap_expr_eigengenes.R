# Heatmap Expression of eigen gene
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")

load("resources/eigenExpr.RData")
eigenExpr <- eigenExpr[["braak1-6"]]
load("resources/braakStages.RData")
load("resources/braakLabels.RData")
# load("resources/hierclust_tree.RData")
load("resources/modules.RData")
library(gplots)
ontology <- read.csv("../ABA_human_processed/Ontology_edited.csv")
load("resources/summaryLabelCorrEG.RData")
labelCor <- do.call(rbind.data.frame, lapply(summaryLabelCorrEG, function(g) g["summary",]))

# heatmap
eigenExpr <- lapply(eigenExpr, function(m){
  rownames(m) <- paste0("M", rownames(m))
  m
})

# IDs and AHBA colors for each sample per donor
sampleInfo <- lapply(donorNames, function(d){
  sampleIds <- read.csv(paste("../ABA_human_processed/sample_info_normalized_microarray_", d, "_2014-11-11.csv", sep = ""))[ , 1]
  info <- ontology[match(sampleIds, ontology$id), ]
  info$color_hex_triplet <- sapply(info$color_hex_triplet, function(c){
    if (nchar(c) == 5) {paste("#0", c, sep = "")} else {paste("#", c, sep = "")}
  })
  info
})

modules <- names(modules[["braak1-6"]])
rowOrder <- order(labelCor$r)
colPal <- c("darkblue", "white", "darkred")
rampcols <- colorRampPalette(colors = colPal, space="Lab")(201)
rowColor <- rampcols[as.numeric(cut(labelCor$r, breaks = 201))][rowOrder]
rowsep <- which(labelCor$r[rowOrder] > 0)[1]

pdf("heatmap_expr_eigengenes.pdf", 8, 6)
lapply(donorNames, function(d){
  # Subselect expression matrices
  samples <- braakStages[[d]]
  samples <- as.logical(apply(samples, 1, sum))
  df <- sampleInfo[[d]][samples,]
  labels <- braakLabels[[d]][samples]
  exprMat <- as.matrix(eigenExpr[[d]])
  colColor <- df$color_hex_triplet
  
  colOrder <- order(labels, -df$graph_order)
  df <- df[colOrder, ]
  exprMat <- exprMat[rowOrder, colOrder]
  # rownames(exprMat)[seq(2, nrow(exprMat), by= 3)] <- paste0("            ", rownames(exprMat)[seq(2, nrow(exprMat), by= 3)])
  # rownames(exprMat)[seq(3, nrow(exprMat), by= 3)] <- paste0("                        ", rownames(exprMat)[seq(3, nrow(exprMat), by= 3)])
  rownames(exprMat)[-which(rownames(exprMat) %in% c("M85", "M10", "M2", "M48"))] <- "" 
  colColor <- colColor[colOrder]
  
  heatmap.2(exprMat, col = rampcols, 
            labCol= df$name, 
            rowsep = rowsep, sepcolor = "black",
            Rowv=FALSE, Colv=FALSE, 
            cexCol = 0.1, cexRow = 1,
            scale = "none", trace = "none", dendrogram = "none", key = FALSE, 
            RowSideColors = rowColor, ColSideColors = colColor,
            main = d,
            margins = c(15, 10))
  
})
dev.off()
