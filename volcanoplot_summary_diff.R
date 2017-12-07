# Volcano plot for differential expression analysis

setwd("C:/Users/dkeo/surfdrive/Parkinson")

library(ggplot2)
library(gridExtra)

source("PD/base_script.R")
load("resources/summaryDiffExpr.RData")

# Extract summary statistics
summaryStats <- lapply(summaryDiffExpr, function(rp){
  as.data.frame(t(sapply(rp, function(g) {
    unlist(g["summary", !colnames(g) %in% c("donors", "weight")])
  })))
})

theme <- theme(legend.position = "none",
               panel.background = element_blank(),
               axis.line = element_line(colour = "black"),
               axis.title =  element_text(size = 12),
               plot.title = element_text(size = 12, face = "bold"))

volcano.plot <- function(tab, rp){
  tab$'logp' <- -log10(tab$benjamini_hochberg)
  tab$info <- tab$benjamini_hochberg < 0.05 & abs(tab$meanDiff) > 1
  ggplot(tab, aes(meanDiff, logp, colour = info)) +
  geom_point(size = 0.5) +
  labs(x = "fold-change", y = "-log10 p-value") +
  ggtitle(paste("Braak stage ", gsub("braak", "", gsub("-", " vs ", rp)))) +
  theme
}

plotll <- lapply(names(summaryStats), function(rp){
  tab <- summaryStats[[rp]]
  p <- volcano.plot(tab, rp)
  name <- paste0("DiffExpr_braak/volcanoplot_", rp, ".pdf")
   pdf(file = name, 4, 3)
   print(p)
   dev.off()
})
# 
# pdf(file = "volcanoplot_summary_diff.pdf", 12, 12)
# grid.arrange(grobs = plotll, top = "Summary differential expression", nrow = 5)
# dev.off()