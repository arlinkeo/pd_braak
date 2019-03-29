# Line plot of dopaminergic genes
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(ggplot2)
# library(reshape2)
brainExpr <- readRDS("../AHBA_Arlin/gene_expr.RDS")
load("resources/braakInfo.RData")

# Dopaminergic genes
genes <- c("SNCA", "SNCB", "GCH1", "TH", "SLC6A3", "SLC18A2")
genes_id <- name2EntrezId(genes)

# Plot function
plot.matrix <- function(expr, title){
  df <- melt(as.matrix(expr))
  colnames(df) <- c("gene", "sample", "expr")
  df$gene <- factor(df$gene, levels = unique(df$gene))
  df$sample <- factor(df$sample, levels = unique(df$sample))
  
  ggplot(df, aes(x=sample, y=expr, group=gene, color=gene)) + 
    geom_point() +
    geom_line() +
    # geom_smooth(aes(group=gene), span = 0.1) +
    # geom_vline(xintercept = vline) +
    ggtitle(title) +
    theme_classic()
}

# Data per donor
expr <- lapply(donorNames, function(d){
  idx <- braak_idx[[d]]
  e <- sapply(idx, function(i){
    e <- brainExpr[[d]][genes_id, i]
    apply(e, 1, median)
  })
  rownames(e) <- genes
  e
})
# labels <- lapply(donorNames, function(d){
#   idx <- unlist(braak_idx[[d]])
#   braakLabels[[d]][idx]
# })

# Plot
pdf("lineplot_dopaminergic_genes.pdf",6, 4)
# Per donor
lapply(donorNames, function(d){
  e <- expr[[d]]
  # l <- labels[[d]]
  # intercepts <- match(c(1:6), l)[-1]-1
  plot.matrix(e, d)
})
# Average median across donor
mean_median <- apply(simplify2array(expr), c(1,2), mean)
plot.matrix(mean_median, "mean median across donors")
median_median <- apply(simplify2array(expr), c(1,2), median)
plot.matrix(median_median, "median median across donors")
max_median <- apply(simplify2array(expr), c(1,2), max)
plot.matrix(max_median, "max median across donors")
min_median <- apply(simplify2array(expr), c(1,2), min)
plot.matrix(min_median, "min median across donors")
dev.off()

# Concatenate data
expr.concat <- Reduce(cbind, expr)
labels.concat <- unlist(labels)
colorder <- order(labels.concat)
expr.concat <- expr.concat[, colorder]
labels.concat <- labels.concat[colorder]
intercepts.concat <- match(c(1:6), labels.concat)[-1]-1

plot.concat <- plot.matrix(expr.concat, intercepts.concat, "dopaminergic genes")
plot.concat
