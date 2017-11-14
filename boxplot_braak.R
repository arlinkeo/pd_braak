# Boxplot of expression in Braak regions for a single gene
setwd("C:/Users/dkeo/surfdrive/Parkinson")
library(ggplot2)

source("PD/base_script.R")
load("../ABA_Rdata/BrainExprNorm.RData")
load("resources/braakLabels.RData") # Braak stage label vectors
load("resources/summaryCorr.RData") # Expression - Braak label correlation
# load("resources/profileBraakGenes.RData") # Binary Braak expression profile

# Get gene's expression
gene.expr <- function(g, d){
  unlist(brainExprNorm[[d]][gene, ])
}

# Default theme for boxplot
theme <- theme(panel.background = element_blank(), panel.grid = element_blank(), 
               axis.line = element_line(colour = "black"),
               legend.title = element_blank())

# Boxplot for a gene
boxplot.gene <- function(g){
  exprll <- lapply(donorNames, function(d) {
    expr <- gene.expr(g, d)
    label <- braakLabels[[d]]
    donor <- paste("Donor", tail(unlist(strsplit(d, split = "donor")), 1))
    donor <- rep(donor, length(expr))
    df <- data.frame(expr, label, donor)
  })
  exprll <- Reduce(rbind, exprll)
  
  r <- format(summaryCorr[[g]]["summary", c("z", "pvalue")], digits = 2)
  title <- paste0(entrezId2Name(g), ", r=", r$z, ", p=", r$pvalue)
  
  p <- ggplot(exprll) + geom_boxplot(aes(x = factor(label), y = expr, fill = factor(donor))) +
    labs(x = "Braak stage", y = "Expression") +
    ggtitle(title) +
    theme
  p
}

plot.pdf <- function(name, genes){
  pdf(name, 8, 6)
  lapply(genes, function(gene){
    p <- boxplot.gene(gene)
    print(p)
  })
  dev.off()
}

# Boxplot for each PD-implicated gene
plot.pdf("boxplot_high_impact_genes.pdf", pdGenesID$hiImpact)
plot.pdf("boxplot_susceptible_genes.pdf", pdGenesID$susceptible)
plot.pdf("boxplot_HLA_genes.pdf", pdGenesID$hla)
plot.pdf("boxplot_lysosome_genes.pdf", pdGenesID$lysosome)