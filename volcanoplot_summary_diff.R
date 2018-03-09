# Volcano plot for differential expression analysis

setwd("C:/Users/dkeo/surfdrive/pd_braak")

library(ggplot2)
library(gridExtra)
library(ggrepel)

source("PD/base_script.R")
load("resources/summaryDiffExpr.RData")

# Extract summary statistics, correct p-values
diffExpr <- lapply(summaryDiffExpr, function(rp){
  # as.data.frame(t(sapply(rp, function(g) {
  #   unlist(g["summary", !colnames(g) %in% c("donors", "weight")])
  # })))
  rp <- do.call(rbind.data.frame, lapply(rp, function(g) g["summary",]))
  rp$benjamini_hochberg <- p.adjust(rp$pvalue, method = "BH")
  rp$'logp' <- -log10(rp$benjamini_hochberg)
  rp
})

# Number of diff. expr. genes for each region pair after correction
diffGenes <- sapply(diffExpr, function(rp){
  order <- order(rp$benjamini_hochberg) # order absolute corr.
  diffGenes1 <- rownames(rp)[order[1:1992]] # top 10% genes
  order <- rev(order(abs(rp$meanDiff))) # order absolute corr.
  diffGenes2 <- rownames(rp)[order[1:1992]] # top 10% genes
  intersect(diffGenes1, diffGenes2)
})

##### Volcano plot #####

theme <- theme(legend.position = "none",
               panel.background = element_blank(),
               axis.line = element_line(colour = "black"),
               axis.title =  element_text(size = 12),
               plot.title = element_text(size = 12, face = "bold"))

xmax <- max(sapply(diffExpr, function(rp){max(rp$mean)}))
xmin <- min(sapply(diffExpr, function(rp){min(rp$mean)}))
ymax <- max(sapply(diffExpr, function(rp){  max(rp$'logp'[is.finite(rp$'logp')])}))


volcano.plot <- function(tab, rp){
  ggplot(tab, aes(meanDiff, logp, colour = info)) +
    geom_point(size = 0.25) +
    scale_colour_manual(values = c("0"="#F8766D", "1"="#00BFC4", "2"="black")) +
    geom_text_repel(label = tab$labels, fontface = "italic", colour = "black", size = 3, nudge_x = 0.2) +
    labs(x = "fold-change", y = "-log10 p-value") +
    scale_x_continuous(limits = c(xmin, xmax), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, ymax), expand = c(0,0)) +
    ggtitle(paste("Braak stage ", gsub("braak", "", gsub("-", " vs ", rp)))) +
    theme
}

plotll <- lapply(names(diffExpr), function(rp){
  tab <- diffExpr[[rp]]
  deg <- diffGenes[[rp]]
  tab$info <- as.numeric(rownames(tab) %in% deg)
  # tab$info <- as.numeric(tab$benjamini_hochberg < 0.05 & abs(tab$meanDiff) > 1.54)
  tab$labels <- entrezId2Name(rownames(tab))
  tab$labels[!tab$labels %in% c("SNCA")] <- ""
  tab[name2EntrezId(c("SNCA")), "info"] <- 2
  tab$info <- as.factor(tab$info)
  order <- order(tab$info)
  tab <- tab[order, ]
  
  p <- volcano.plot(tab, rp)
  name <- paste0("DiffExpr_braak/volcano_ttest/volcanoplot_", rp, ".pdf")
  pdf(file = name, 4, 3)
  print(p)
  dev.off()
})
# 
# pdf(file = "volcanoplot_summary_diff.pdf", 12, 12)
# grid.arrange(grobs = plotll, top = "Summary differential expression", nrow = 5)
# dev.off()