# Volcano plot for differential expression analysis

setwd("C:/Users/dkeo/surfdrive/pd_braak")

library(ggplot2)
library(gridExtra)
library(ggrepel)

source("PD/base_script.R")
load("resources/summaryDiffExpr_ztest.RData")

# Extract summary statistics
summaryStats <- lapply(summaryDiffExpr, function(rp){
  as.data.frame(t(sapply(rp, function(g) {
    unlist(g["summary", !colnames(g) %in% c("donors", "weight")])
  })))
})

# Number of diff. expr. genes for each region pair
sapply(summaryStats, function(rp){
  sum(rp$benjamini_hochberg < 0.05 & abs(rp$meanDiff) > 1) # significance and 2-fold change
})

##### Volcano plot #####

theme <- theme(legend.position = "none",
               panel.background = element_blank(),
               axis.line = element_line(colour = "black"),
               axis.title =  element_text(size = 12),
               plot.title = element_text(size = 12, face = "bold"))

volcano.plot <- function(tab, rp){
  tab$'logp' <- -log10(tab$benjamini_hochberg)
  # tab$info <- tab$benjamini_hochberg < 0.05 & abs(tab$meanDiff) > 1
  ggplot(tab, aes(meanDiff, logp, colour = info)) +
    geom_point(size = 0.5) +
    scale_colour_manual(values = c("0"="#F8766D", "1"="#00BFC4", "2"="black")) +
    geom_text_repel(label = tab$labels, fontface = "italic", colour = "black", size = 3, nudge_x = 0.2) +
    labs(x = "fold-change", y = "-log10 p-value") +
    ggtitle(paste("Braak stage ", gsub("braak", "", gsub("-", " vs ", rp)))) +
    theme
}

plotll <- lapply(names(summaryStats), function(rp){
  tab <- summaryStats[[rp]]
  tab$info <- as.numeric(tab$benjamini_hochberg < 0.05 & abs(tab$meanDiff) > 1)
  tab$labels <- entrezId2Name(rownames(tab))
  tab$labels[!tab$labels %in% c("SNCA", "GCH1")] <- ""
  tab[name2EntrezId(c("SNCA", "GCH1")), "info"] <- 2
  tab$info <- as.factor(tab$info)
  
  p <- volcano.plot(tab, rp)
  name <- paste0("DiffExpr_braak/volcano_ztest/volcanoplot_", rp, ".pdf")
  pdf(file = name, 4, 3)
  print(p)
  dev.off()
})
# 
# pdf(file = "volcanoplot_summary_diff.pdf", 12, 12)
# grid.arrange(grobs = plotll, top = "Summary differential expression", nrow = 5)
# dev.off()