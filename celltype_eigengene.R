# Correction for cell-type expression
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(ggplot2)
library(gplots)
library(reshape2)
load("resources/braakInfo.RData")
load("../ABA_Rdata/BrainExpr.RData")

# Cell-type genes
celltypes <- sapply(c("Neurons", "Astrocytes", "Oligodendrocytes", "Microglia", "Endothelial_cells"), function(type){
  file = paste0("brainscope_celltypes/", type, ".txt")
  as.character(read.csv(file, header = TRUE)$entrez_id)
}, simplify = FALSE)

# # PSEA cell-type markers
# celltypes <- list(Neuron = c("NEFL", "ENO2", "SLC12A5", "KCNQ2", "SCN3A"),
#                   Astrocyte = c("GFAP", "AQP4", "GJA1"),
#                   Oligodendrocyte = c("MOG", "MAG", "MOBP", "MBP"),
#                   Microglia = c("CD37", "CD53"),
#                   Endothelial = c("PECAM1"))
# celltypes <- lapply(celltypes, name2EntrezId)

# Cell-type mean gene expression for all samples (whole brain)
mean_celltype <- lapply(donorNames, function(d){
  t(sapply(celltypes, function(ct){
    x <- brainExpr[[d]][ct, ]
    apply(x, 2, mean)
  }))
})

# Cell-type eigengene
eg_celltype <- lapply(donorNames, function(d){
  t(sapply(celltypes, function(ct){
    x <- t(brainExpr[[d]][ct, ])
    eg <- prcomp(x, center = FALSE)$x[, 1]# 1st PC (eigen gene expr)
    mean <- apply(x, 1, mean)
    if (cor(eg, mean) > 0) eg else -eg # flip sign of eigen gene based on the data
  }))
})

plot.ref.signal <- dget("PD/plot.ref.signal.R")

# Plot celltype eigengene expression
plots <- lapply(donorNames, function(d){
  idx <- Reduce(c, braak_idx[[d]]) # for samples
  expr <- eg_celltype[[d]][, idx]
  plot.ref.signal(expr)
})
pdf("celltype_eigengene_centered_scaled.pdf", 8, 6)
plots
dev.off()

# Fit linear model all samples
lm_fit <- lapply(donorNames, function(d){
  braak <- braakLabels[[d]] # braak labels
  ct <- t(eg_celltype[[d]]) # samples x cell-types
  expr <- t(brainExpr[[d]]) # samples x genes
  fit <- lm(expr ~ ct)
})

# Fit linear model with Braak samples only
lm_fit <- lapply(donorNames, function(d){
  idx <- unlist(braak_idx[[d]]) # idx for samples
  braak <- braakLabels[[d]][idx] # braak labels
  ct <- t(eg_celltype[[d]][, idx]) # samples x cell-types
  expr <- t(brainExpr[[d]][, idx]) # samples x genes
  fit <- lm(expr ~ ct)
})

# Use residuals as cell-corrected expression
expr_celltype_corrected <- lapply(donorNames, function(d){
  fit <- lm_fit[[d]]
  t(fit$residuals)
})
saveRDS(expr_celltype_corrected, file = "resources/expr_corrected_lm_eg_braaksamples.rds")

########## Heatmap Before and after correction ##########
ct_genes <- unlist(celltypes)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
ggcolor <- gg_color_hue(5)
celltypeColors <- sapply(ct_genes, function(g){
  if (g %in% celltypes$Neurons) ggcolor[1]
  else if (g %in% celltypes$Astrocytes) ggcolor[2]
  else if (g %in% celltypes$Oligodendrocytes) ggcolor[3]
  else if (g %in% celltypes$Microglia) ggcolor[4]
  else if (g %in% celltypes$Endothelial_cells) ggcolor[5]
  else "gray"
})

lapply(donorNames, function(d){
  idx <- unlist(braak_idx[[d]]) # for samples
  info <- sampleInfo[[d]][idx, ]
  
  # Scale expression across all Braak samples
  expr <- scale(t(brainExpr[[d]][ct_genes, idx])) # samples x genes
  expr2 <- scale(t(expr_celltype_corrected[[d]][ct_genes,]))
  
  colPal <- c("blue", "white", "red")
  rampcols <- colorRampPalette(colors = colPal, space="Lab")(100)
  rampcols <- c(rep(col2hex(colPal[1]), 50), rampcols, rep(col2hex(colPal[3]), 50))
  
  file = paste0("celltype_correction_heatmap/neuroncorrected_", d, ".pdf")
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
  samples <- unlist(braak_idx[[d]])
  label <- braakLabels[[d]][samples]
  ct <- mean_celltype[[d]][, samples]
  expr <- brainExpr[[d]][, samples]
  
  g <- "114"#name2EntrezId("SNCA")
  
  fit <- lm_fit[[d]]
  coef <- fit$coefficients[, "6622"]
  
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
