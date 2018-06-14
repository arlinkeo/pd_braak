# Boxplot of expression in Braak regions for a single gene
setwd("C:/Users/dkeo/surfdrive/pd_braak")
library(ggplot2)
source("PD/base_script.R")
load("../ABA_Rdata/BrainExpr.RData")
# load("../ABA_Rdata/BrainExprNorm.RData")
load("resources/braakLabels.RData") # Braak stage label vectors
load("resources/summaryLabelCorr.RData")

# Get gene's expression
gene.expr <- function(g, d){
  unlist(brainExpr[[d]][g, ])
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
  exprll <- exprll[exprll$label != "0", ]#Remove Braak 0
  
  r <- format(summaryLabelCorr[[g]]["summary", c("r", "pvalue")], digits = 2)
  title <- paste0(entrezId2Name(g), ", r=", r$r)
  
  p <- ggplot(exprll) + geom_boxplot(aes(x = factor(label), y = expr, fill = factor(donor))) +
    labs(x = "Braak stage", y = "Expression (log2-transformed)") +
    ggtitle(title) +
    theme
  p
}

plot.pdf <- function(name, genes){
  pdf(name, 8, 5)
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
