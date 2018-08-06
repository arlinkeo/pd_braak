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

# Cell-type mean gene expression for all samples (whole brain)
mean_celltype <- lapply(donorNames, function(d){
  t(sapply(celltypes, function(ct){
    x <- brainExpr[[d]][ct, ]
    apply(x, 2, mean)
  }))
})

# Plot celltype eigengene expression
plots <- lapply(donorNames, function(d){
  idx <- Reduce(c, braak_idx[[d]]) # for samples
  expr <- mean_celltype[[d]][, idx]
  colnames(expr) <- make.names(colnames(expr), unique = TRUE)
  
  df <- melt(as.matrix(expr))
  colnames(df) <- c("celltype", "sample", "expr")
  df$sample <- factor(df$sample, levels = unique(df$sample))
  
  intercepts <- match(c(1:6), braakLabels[[d]][idx])[-1]
  
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
pdf("celltype_mean.pdf", 8, 6)
plots
dev.off()

# Fit linear model
lm_fit <- lapply(donorNames, function(d){
  idx <- Reduce(c, braak_idx[[d]]) # idx for samples
  braak <- braakLabels[[d]][idx] # braak labels
  ct <- t(mean_celltype[[d]][, idx]) # samples x cell-types
  expr <- t(brainExpr[[d]][, idx]) # samples x genes
  lm(expr ~ ct : braak)
})

# Use residuals as cell-corrected expression
expr_celltype_corrected <- lapply(donorNames, function(d){
  fit <- lm_fit[[d]]
  t(fit$residuals)
})
saveRDS(expr_celltype_corrected, file = "resources/expr_celltype_corrected.rds")

########## Heatmap Before and after correction ##########
ct_genes <- unlist(celltypes)
celltypeColors <- sapply(ct_genes, function(g){
  if (g %in% celltypes$Astrocytes) "red"
  else if (g %in% celltypes$Endothelial_cells) "yellow"
  else if (g %in% celltypes$Microglia) "pink"
  else if (g %in% celltypes$Neurons) "blue"
  else if (g %in% celltypes$Oligodendrocytes) "green"
  else "gray"
})

lapply(donorNames, function(d){
  idx <- Reduce(c, braak_idx[[d]]) # for samples
  info <- sampleInfo[[d]][idx, ]
  
  # Scale expression across all Braak samples
  expr <- scale(t(brainExpr[[d]][ct_genes, idx])) # samples x genes
  expr2 <- scale(t(expr_celltype_corrected[[d]][ct_genes,]))
  
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

# Scatter plots
lapply(donorNames, function(d){
  samples <- Reduce(c, braak_idx[[d]])
  label <- braakLabels[[d]][samples]
  ct <- mean_celltype[[d]][, samples]
  expr <- brainExpr[[d]][, samples]
  
  g <- "114"#name2EntrezId("SNCA")
  
  lapply(names(celltypes), function(c){
    ref <- ct[c, ] # Cell-type expression
    gene <- unlist(expr[g, ])
    df <- data.frame(ref, gene, label)
    ggplot(df, aes(x=ref, y=gene, color = label, shape = label)) +
      geom_point() +
      geom_smooth(method = lm, se = FALSE, fullrange = TRUE) +
      labs(x = c, y = entrezId2Name(g)) +
      ggtitle(d) +
      theme(panel.background = element_blank(), 
            panel.grid = element_blank(), 
            axis.line = element_line(colour = "black"))
  })
})
