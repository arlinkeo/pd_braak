# Correction for cell-type expression
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(ggplot2)
library(reshape2)
load("resources/braakInfo.RData")
load("../ABA_Rdata/BrainExprNorm.RData")

# Cell-type genes
celltypes <- sapply(c("Astrocytes", "Endothelial_cells", "Microglia", "Neurons", "Oligodendrocytes"), function(type){
  file = paste0("brainscope_celltypes/", type, ".txt")
  as.character(read.csv(file, header = TRUE)$entrez_id)
}, simplify = FALSE)

# eigen gene function for data matrix (samples x genes)
eigen.gene <- function(x){
  eg <- prcomp(x)$x[, 1]# 1st PC (eigen gene expr)
  mean <- apply(x, 1, mean)
  if (cor(eg, mean) > 0) eg else -eg # flip sign of eigen gene based on the data
}

# Cell-type eigen gene expression
eg_celltype <- lapply(donorNames, function(d){
  l <- t(sapply(celltypes, function(ct){
    x <- brainExprNorm[[d]][ct, ]
    eigen.gene(t(x))
  }))
})

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

# Plot celltype eigengene expression
plots <- lapply(donorNames, function(d){
  samples <- braakLabels[[d]] != 0
  expr <- eg_celltype[[d]][, samples]
  # expr <- rbind(expr, sum = apply(expr, 2, sum))
  labels <- braakLabels[[d]][samples]
  
  graph_order <- sampleInfo[[d]][samples, "graph_order"]
  colOrder <- order(labels, -graph_order)
  expr <- expr[, colOrder]
  colnames(expr) <- make.names(colnames(expr), unique = TRUE)
  
  df <- melt(expr)
  colnames(df) <- c("celltype", "sample", "expr")
  df$sample <- factor(df$sample, levels = unique(df$sample))
  
  intercepts <- match(c(1:6), labels[colOrder])[-1]
  
  ggplot(df, aes(x=sample, y=expr, color=celltype)) + 
    geom_point() +
    geom_smooth(aes(group=celltype)) +
    geom_vline(xintercept = intercepts) +
    ggtitle(d) +
    theme(panel.background = element_blank(), 
          panel.grid = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
})

pdf("celltype_eigengene.pdf", 8, 6)
plots
dev.off()