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
rowOrder <- order(-labelCor$r)
colPal <- c("darkblue", "white", "darkred")
rampcols <- colorRampPalette(colors = colPal, space="Lab")(200)
  
pdf("heatmap_expr_eigengenes.pdf", 8, 6)

lapply(donorNames, function(d){
  # Subselect expression matrices
  samples <- braakStages[[d]]
  # samples <- as.logical(bitwOr(samples[, "braak1"], samples[, "braak6"]))
  # samples <- as.logical(samples[, b])
  samples <- as.logical(apply(samples, 1, sum))
  df <- sampleInfo[[d]][samples,]
  labels <- braakLabels[[d]][samples]
  exprMat <- as.matrix(eigenExpr[[d]])
  
  colOrder <- order(labels, -df$graph_order)
  df <- df[colOrder, ]
  exprMat <- exprMat[rowOrder, colOrder]
  

  heatmap.2(exprMat, col = rampcols, 
            labCol = df$name, 
            Rowv=FALSE, Colv=FALSE, 
            cexCol = 0.1, cexRow = 0.5,
            scale = "none", trace = "none", dendrogram = "none", key = FALSE, 
            ColSideColors = df$color_hex_triplet,
            main = paste0("Expression of modules in Braak 1-6"),
            margins = c(15, 5))
})

dev.off()
