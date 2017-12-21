# Heatmap subset 500 genes of average co-expression and summary co-expression
setwd("C:/Users/dkeo/surfdrive/Parkinson")
source("PD/base_script.R")

library(ggplot2)
library(reshape2)
load("resources/summaryGeneCoexpr_subset.RData")
load("resources/avgCor.RData")
load("resources/subsetGenes.RData")

subsetNames <- entrezId2Name(subsetGenes)

avgCor <- avgCor[subsetGenes, subsetGenes]
rownames(avgCor) <- subsetNames
colnames(avgCor) <- subsetNames
diag(avgCor) <- 1

sumCor <-  sapply(subsetGenes, function(g1){
  sapply(subsetGenes, function(g2){
    pair <- paste(g1, g2, sep = "-")
    ifelse(pair %in% names(summaryGeneCoexpr), summaryGeneCoexpr[[pair]]["summary", "r"], NA)
  })
})
rownames(sumCor) <- subsetNames
colnames(sumCor) <- subsetNames
diag(sumCor) <- 1
m <- sumCor
m[upper.tri(m)] <- t(m)[upper.tri(m)]

distance <- dist(avgCor)
t <- hclust(distance)
order <- subsetNames[t$order]

plot.heatmap <- function(mat, order, t){
  tab <- melt(mat[order, order])
  tab$Var1 <- factor(tab$Var1, levels = order)
  tab$Var2 <- factor(tab$Var2, levels = order)
  ggplot(tab, aes(Var1, Var2)) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", name = "r") +
    labs(x = "", y = "") +
    ggtitle(paste0("Gene co-expression (", t, " correlation", ")")) +
    theme(axis.text = element_text(face = "italic", size = 3), axis.text.x = element_text(angle = 90))
}

p1 <- plot.heatmap(avgCor, order, "average")
p2 <- plot.heatmap(m, order, "summary")

pdf("heatmap_coexpr500.pdf", 8, 6.5)
p1
p2
dev.off()
