# Line plot of dopaminergic genes
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(ggplot2)
library(reshape2)
brainExpr <- readRDS("../AHBA_Arlin/gene_expr.RDS")
load("resources/braakInfo.RData")
load("resources/eigenExpr.RData")

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

# Function get expression of genes and modules
get.expr <- function(g, m){
  lapply(donorNames, function(d){
    idx <- braak_idx[[d]]
    labels <- braakLabels[[d]][unlist(idx)]
    idx_eg <- sapply(names(idx), function(i){ #index of sorted labels
      which(labels == gsub("R","", i))
    }, simplify = FALSE) 
    e <- sapply(names(idx), function(i){
      e <- brainExpr[[d]][g, idx[[i]]] # expr of genes
      em <- eigenExpr[[d]][m, idx_eg[[i]]]# expr of module
      e <- rbind(em, e)
      apply(e, 1, median)
    })
    rownames(e) <- c(m, g)
    e
  })
}

# Data of dopaminergic genes & module
genes <- c("SNCA", "GCH1", "TH", "SLC6A3", "SLC18A2")
genes_id <- name2EntrezId(genes)
module <- c("M127")
expr_dopa <- get.expr(genes_id, module)

pdf("lineplot_dopaminergic_genes.pdf", 6, 4)
# # Per donor
# lapply(donorNames, function(d){
#   e <- expr[[d]]
#   plot.matrix(e, d)
# })
# Average median across donor
mean_median <- apply(simplify2array(expr_dopa), c(1,2), mean)
plot.matrix(mean_median, "mean median across donors")
# median_median <- apply(simplify2array(expr_dopa), c(1,2), median)
# plot.matrix(median_median, "median median across donors")
# max_median <- apply(simplify2array(expr_dopa), c(1,2), max)
# plot.matrix(max_median, "max median across donors")
# min_median <- apply(simplify2array(expr_dopa), c(1,2), min)
# plot.matrix(min_median, "min median across donors")
dev.off()

# Data of blood oxygen genes
genes <- c("HBD", "HBB", "HBA1", "HBA2", "OASL")
genes_id <- name2EntrezId(genes)
module <- c("M47")
expr_bloodoxy <- get.expr(genes_id, module)

pdf("lineplot_bloodoxygen_genes.pdf", 6, 4)
mean_median <- apply(simplify2array(expr_bloodoxy), c(1,2), mean)
plot.matrix(mean_median, "mean median across donors")
dev.off()