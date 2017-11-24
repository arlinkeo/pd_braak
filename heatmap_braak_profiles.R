# Heatmap of expression of genes grouped by Braak profile
setwd("C:/Users/dkeo/surfdrive/Parkinson")
source("PD/base_script.R")

library(ggplot2)
library(reshape2)
# library(gridExtra)
# library(grid)

load("../ABA_Rdata/BrainExprNorm.RData")
load("resources/braakProfile.RData")
# load("resources/profileList.RData")
load("resources/braakLabels.RData") # Braak stage label vectors
load("resources/braakGenes.RData")

pdf("heatmap_braak_genes.pdf", 8, 9)

lapply(donorNames, function(d){
  
  # Subselect expression matrices
  labels <- braakLabels[[d]]
  samples <- names(labels)[labels != 0]
  exprMat <- brainExprNorm[[d]]
  exprMat <- exprMat[braakGenes, samples]
  
  # Plot heatmap
  mat <- as.data.frame(t(exprMat))
  labels <- labels[samples]
  tab <- cbind(sample = rownames(mat), braakstage = labels, mat)
  tab.m <- melt(tab, id.vars = c("sample", "braakstage"), variable.name = "gene", value.name = "expr")
  tab.m$braakstage <- factor(tab.m$braakstage, levels = c(1:6)) # keep order
  profiles <- braakProfile[braakGenes]
  tab.m$profile <- factor(rep(profiles, sum(samples)), levels = unique(profiles))
  
  p <- ggplot(tab.m, aes(sample, gene, group = braakstage)) +
    geom_tile(aes(fill = expr)) +
    facet_grid(profile~braakstage, scales = "free", space = "free") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", name = "Expression") +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black", size = 0.2)
    ) +
    labs(x = "Braak stage", y = "Expression profile") +
    ggtitle(paste("Expression of Braak related genes in Donor", tail(unlist(strsplit(d, split = "donor")), 1)))
  print(p)
  
})

dev.off()