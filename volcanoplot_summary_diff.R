# Volcano plot for differential expression analysis

setwd("C:/Users/dkeo/surfdrive/Parkinson")

library(ggplot2)
library(gridExtra)

source("PD/base_script.R")
load("resources/summaryMeanDiff.RData")

genesPerBraak <- lapply(c(braakNames, braakNamesMerged1), function(bs){
  as.data.frame(t(sapply(summaryMeanDiff$nonBraakA2, function(g) {
    unlist(g[[bs]]["summary", c("meanDiff", "varDiff", "lower95", "upper95", "benjamini_hochberg")])
  })))
})

theme <- theme(legend.position = "none",
               panel.background = element_blank(),
               axis.line = element_line(colour = "black"),
               axis.title =  element_text(size = 12),
               plot.title = element_text(size = 12, face = "bold"))

volcano.plot <- function(tab, bs){
  ggplot(tab, aes(meanDiff, logp, colour = info)) +
  geom_point(size = 0.5) +
  labs(x = "fold-change", y = "-log10 p-value") +
  ggtitle(paste("Braak stage", tail(unlist(strsplit(bs, split = "braak")), 1))) +
  theme
}

plotll <- lapply(braakNamesMerged1, function(bs){
  tab <- genesPerBraak[[bs]]
  tab$'logp' <- -log10(tab$benjamini_hochberg)
  tab$info <- tab$benjamini_hochberg < 0.05 & abs(tab$meanDiff) > 1
  volcano.plot(tab, bs)
})

pdf(file = "volcanoplot_summary_diff.pdf", 12, 4)
grid.arrange(grobs = plotll, top = "Summary differential expression", nrow = 1)
dev.off()