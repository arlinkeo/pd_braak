# Volcano plot for differential expression analysis

setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")

library(ggplot2)
library(gridExtra)
library(ggrepel)
load("resources/braakGenes.RData")
braakGenes <- unlist(braakGenes)
load("resources/summaryDiffExpr.RData")

# Extract summary statistics
diffExpr <- lapply(summaryDiffExpr, function(rp){
  do.call(rbind.data.frame, lapply(rp, function(g) g["summary",]))
})

# # Number of diff. expr. genes for each region pair after correction
# diffGenes <- sapply(diffExpr, function(rp){
#   order <- order(rp$benjamini_hochberg) # order absolute corr.
#   diffGenes1 <- rownames(rp)[order[1:1992]] # top 10% genes
#   order <- rev(order(abs(rp$meanDiff))) # order absolute corr.
#   diffGenes2 <- rownames(rp)[order[1:1992]] # top 10% genes
#   intersect(diffGenes1, diffGenes2)
# })

##### Volcano plot #####

theme <- theme(legend.position = "none",
               panel.background = element_blank(),
               axis.line = element_line(colour = "black"),
               axis.title =  element_text(size = 12),
               plot.title = element_text(size = 12, face = "bold"))

# corrected and log-transformed p-values
diffExpr <- lapply(diffExpr, function(tab){
  tab$benjamini_hochberg <- p.adjust(tab$pvalue, method = "BH")
  tab$'logp' <- -log10(tab$benjamini_hochberg)
  tab
})

# axis limits for all plots
ctab <- Reduce(rbind, diffExpr)
xmax <- max(ctab$meanDiff)
xmin <- min(ctab$meanDiff)
ymax <- ceiling(max(ctab$'logp'[is.finite(tab$'logp')]))

plotll <- lapply(names(diffExpr), function(rp){
  tab <- diffExpr[[rp]]
  tab$info <- as.numeric(rownames(tab) %in% braakGenes)
  
  # tab$labels <- entrezId2Name(rownames(tab))
  # pd <- c("DNAJC13","SNCA","GCH1","INPP5F", "ASH1L", "ITPKB", "ELOVL7", "ZNF184", "SCARB2", "DDRGK1")
  # tab$labels[!tab$labels %in% pd] <- ""
  # tab[name2EntrezId(pd), "info"] <- 2
  
  # Plotting order of data points 
  tab$info <- as.factor(tab$info)
  order <- order(tab$info)
  tab <- tab[order, ]

  p <- ggplot(tab, aes(meanDiff, logp, colour = info)) +
    geom_point(size = 0.25, alpha = 0.3) +
    scale_colour_manual(values = c("0"="#F8766D", "1"="#00BFC4", "2"="black")) +
    # geom_text_repel(label = tab$labels, fontface = "italic", colour = "black", size = 3, nudge_x = 0.2) +
    labs(x = "fold-change", y = "-log10 p-value") +
    scale_x_continuous(limits = c(xmin, xmax), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, ymax), expand = c(0,0)) +
    ggtitle(paste("Braak stage ", gsub("braak", "", gsub("-", " vs ", rp)))) +
    theme
  
  name <- paste0("DiffExpr_braak/volcano_ttest/volcanoplot_", rp, ".pdf")
  pdf(file = name, 4, 3)
  print(p)
  dev.off()
})
# 
# pdf(file = "volcanoplot_summary_diff.pdf", 12, 12)
# grid.arrange(grobs = plotll, top = "Summary differential expression", nrow = 5)
# dev.off()