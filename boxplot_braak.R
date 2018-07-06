# Boxplot of expression in Braak regions for a single gene
setwd("C:/Users/dkeo/surfdrive/pd_braak")
library(ggplot2)
source("PD/base_script.R")
load("../ABA_Rdata/BrainExpr.RData")
# load("../ABA_Rdata/BrainExprNorm.RData")
load("resources/braakInfo.RData") # Braak stage label vectors
load("resources/summaryLabelCorr.RData")

# Get gene's expression
gene.expr <- function(g, d){
  unlist(brainExpr[[d]][g, ])
}

# Default theme for boxplot
theme <- theme(panel.background = element_blank(), panel.grid = element_blank(), 
               axis.line = element_line(colour = "black"),
               legend.title = element_blank(),
               legend.key = element_blank())

# Boxplot for a gene
boxplot.gene <- function(g){
  df <- lapply(donorNames, function(d) {
    expr <- gene.expr(g, d)
    label <- braakLabels[[d]]
    donor <- d
    donor <- rep(donor, length(expr))
    data.frame(expr, label, donor)
  })
  df <- Reduce(rbind, df)
  df <- df[df$label != "0", ]#Remove Braak 0
  df$label <- factor(df$label, levels = sort(unique(df$label)))
  df$donor <- factor(df$donor, levels = unique(df$donor))
  
  r <- format(summaryLabelCorr[[g]]["summary", c("r", "pvalue")], digits = 2)
  title <- paste0(entrezId2Name(g), ", r=", r$r)
  
  p <- ggplot(df) + 
    geom_boxplot(aes(x = label, y = expr, alpha = donor, fill = label)) +
    labs(x = "Braak stage", y = "Expression (log2-transformed)") +
    guides(alpha=guide_legend(override.aes=list(fill=hcl(c(15,195),100,0), colour=NA))) +
    scale_alpha_discrete(labels = gsub("donor", "Donor ", donorNames)) +
    scale_fill_manual(values = unname(braakColors), guide = FALSE) +
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