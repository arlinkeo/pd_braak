# Correction for cell-type expression
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(ggplot2)
library(gplots)
library(reshape2)
load("resources/braakInfo.RData")
load("../ABA_Rdata/BrainExpr.RData")

# Cell-type genes
celltypes <- sapply(c("Astrocytes", "Endothelial_cells", "Microglia", "Neurons", "Oligodendrocytes"), function(type){
  file = paste0("brainscope_celltypes/", type, ".txt")
  as.character(read.csv(file, header = TRUE)$entrez_id)
}, simplify = FALSE)

# Neurons and non-neurons
# celltypes <- list(nonneurons = Reduce(c, celltypes[-4]), neurons = celltypes[[4]])

# eigen gene function for data matrix (samples x genes)
eigen.gene <- function(x){
  eg <- prcomp(x, center = FALSE)$x[, 1]# 1st PC (eigen gene expr)
  mean <- apply(x, 1, mean)
  if (cor(eg, mean) > 0) eg else -eg # flip sign of eigen gene based on the data
}

# Cell-type eigen gene expression
eg_celltype <- lapply(donorNames, function(d){
  t(sapply(celltypes, function(ct){
    x <- brainExpr[[d]][ct, ]
    eigen.gene(t(x))
  }))
})

# Plot celltype eigengene expression
plots <- lapply(donorNames, function(d){
  samples <- braakLabels[[d]] != 0
  expr <- eg_celltype[[d]][, samples]
  labels <- braakLabels[[d]][samples]
  
  graph_order <- sampleInfo[[d]][samples, "graph_order"]
  colOrder <- order(labels, -graph_order)
  expr <- expr[, colOrder]
  colnames(expr) <- make.names(colnames(expr), unique = TRUE)
  
  df <- melt(as.matrix(expr))
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

# Fit linear model and use residuals as neuron-corrected expression
expr_neuroncorrected <- lapply(donorNames, function(d){
  expr_ct <- as.data.frame(t(eg_celltype[[d]]))
  expr <- brainExpr[[d]]
  expr <-t(expr)
  fit <- lm(expr ~ 
              expr_ct$Astrocytes + 
              expr_ct$Endothelial_cells + 
              expr_ct$Microglia + 
              expr_ct$Neurons + 
              expr_ct$Oligodendrocytes)
  t(fit$residuals)
})
saveRDS(expr_neuroncorrected, file = "resources/expr_neuroncorrected.rds")

########## Heatmap Before and after correction ##########
genes <- unlist(celltypes)
celltypeColors <- sapply(genes, function(g){
  if (g %in% celltypes$Astrocytes) "red"
  else if (g %in% celltypes$Endothelial_cells) "yellow"
  else if (g %in% celltypes$Microglia) "pink"
  else if (g %in% celltypes$Neurons) "blue"
  else if (g %in% celltypes$Oligodendrocytes) "green"
  else "gray"
})

lapply(donorNames, function(d){
  info <- sampleInfo[[d]]
  expr <- scale(t(brainExpr[[d]]))
  expr2 <- scale(t(expr_neuroncorrected[[d]]))
  
  rowOrder <- order(-info$graph_order)
  info <- info[rowOrder, ]
  expr <- expr[rowOrder, genes]
  expr2 <- expr2[rowOrder, genes]
  
  colPal <- c("blue", "white", "red")
  rampcols <- colorRampPalette(colors = colPal, space="Lab")(100)
  rampcols <- c(rep(col2hex(colPal[1]), 50), rampcols, rep(col2hex(colPal[3]), 50))
  
  file = paste0("celltype_correction_heatmap_", d, ".pdf")
  pdf(file, 8, 9)
  heatmap.2(expr, col = rampcols, 
            labRow = info$acronym, labCol = entrezId2Name(colnames(expr)), 
            Rowv=FALSE, Colv=FALSE, 
            cexCol = .5, cexRow = .5,
            scale = "none", trace = "none", dendrogram = "none", #key = FALSE, 
            RowSideColors = info$color_hex_triplet, ColSideColors = celltypeColors,
            main = paste0("Expression of cell types in ", d, " (uncorrected)"),
            margins = c(5, 5))
  heatmap.2(expr2, col = rampcols, 
            labRow = info$acronym, labCol = entrezId2Name(colnames(expr)), 
            Rowv=FALSE, Colv=FALSE, 
            cexCol = .5, cexRow = .5,
            scale = "none", trace = "none", dendrogram = "none", #key = FALSE, 
            RowSideColors = info$color_hex_triplet, ColSideColors = celltypeColors,
            main = paste0("Expression of cell types in ", d, " (corrected)"),
            margins = c(5, 5))
  dev.off()
})