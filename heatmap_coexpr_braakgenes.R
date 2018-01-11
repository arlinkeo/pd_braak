# Heatmap co-expression of genes
setwd("C:/Users/dkeo/surfdrive/pd_braak")
library(ggplot2)
library(reshape2)

source("PD/base_script.R")
load("resources/avgCor.RData")
load("resources/braakGenes.RData")

# Get correlated genes |r>0.5|
genes <- braakGenes
braakGenesCoexpr <- avgCor[genes, genes]
rownames(braakGenesCoexpr)

#Cluster
distance <- dist(braakGenesCoexpr)
t <- hclust(distance)
order <- genes[t$order]

# Cut tree
membership <- cutree(t, k = 2)
sum(membership == 1)

# Heat stats-function with dendogram
heatmap(braakGenesCoexpr, col = colorRampPalette(c("blue", "white", "red"))(n = 200))

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

p1 <- plot.heatmap(braakGenesCoexpr, order, "average")

pdf("heatmap_coexpr_braakgenes.pdf", 8, 6.5)
p1
dev.off()

