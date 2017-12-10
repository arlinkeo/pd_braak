# Volcano plot for differential expression analysis

setwd("C:/Users/dkeo/surfdrive/Parkinson")

library(ggplot2)
library(gridExtra)
library(ggrepel)

source("PD/base_script.R")
load("resources/summaryDiffExprNonBraak.RData")

genesPerBraak <- lapply(c(braakNames, braakNamesMerged1), function(bs){
  as.data.frame(t(sapply(summaryMeanDiff$nonBraakA2, function(g) {
    unlist(g[[bs]]["summary", c("meanDiff", "varDiff", "lower95", "upper95", "benjamini_hochberg")])
  })))
})

theme <- theme(legend.position = "none",
               axis.title =  element_text(size = 12),
               plot.title = element_text(size = 12, face = "bold"))

volcano.plot <- function(tab, bs, labels){
  tab <- tab[order(tab$info), ]
  ggplot(tab, aes(meanDiff, logp, colour = info)) +
    geom_point(size = 1) +
    scale_colour_manual(values = c("0"="#F8766D", "1"="#00BFC4", "2"="black")) +
    geom_text_repel(label = tab$labels, fontface = "italic", colour = "black", size = 3, nudge_x = 0.2) +
    labs(x = "fold-change", y = "-log10 p-value") +
    ggtitle(paste("Braak stage", tail(unlist(strsplit(bs, split = "braak")), 1))) +
    theme
}

plotll <- lapply(braakNamesMerged1, function(bs){
  tab <- genesPerBraak[[bs]]
  tab$'logp' <- -log10(tab$benjamini_hochberg)
  tab$info <- as.numeric(tab$benjamini_hochberg < 0.05 & abs(tab$meanDiff) > 1)
  tab$labels <- entrezId2Name(rownames(tab))
  tab$labels[!tab$labels %in% c("SNCA", "GCH1")] <- ""
  tab[name2EntrezId(c("SNCA", "GCH1")), "info"] <- 2
  tab$info <- as.factor(tab$info)
  volcano.plot(tab, bs)
})

pdf(file = "volcanoplot_summary_diff_nonbraak.pdf", 12, 4)
grid.arrange(grobs = plotll, top = "Summary differential expression", nrow = 1)
dev.off()
