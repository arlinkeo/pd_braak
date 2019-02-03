# Boxplot of expression in Braak regions for a single gene
setwd("C:/Users/dkeo/surfdrive/pd_braak")
library(ggplot2)
library(reshape2)
source("PD/base_script.R")
brainExpr <- readRDS("../AHBA_Arlin/gene_expr.RDS")
load("resources/braakInfo.RData") # Braak stage label vectors
load("resources/summaryLabelCor.RData")
load("resources/braakGenes.RData")

# Default theme for boxplot
theme <- theme(panel.background = element_blank(), panel.grid = element_blank(), 
               axis.line = element_line(colour = "black"),
               legend.title = element_blank(),
               legend.key = element_blank())

# Boxplot function
box.plot <- function(df, title){
  p <- ggplot(df) + 
    geom_boxplot(aes(x = label, y = expr, alpha = donor, fill = label), outlier.size = 1) +
    labs(x = "Brain region", y = "Expression (log2-transformed)") +
    guides(alpha=guide_legend(override.aes=list(fill=hcl(c(15,195),100,0), colour=NA))) +
    scale_alpha_discrete(labels = gsub("donor", "Donor ", donorNames)) +
    scale_fill_manual(values = unname(braakColors), guide = FALSE) +
    ggtitle(title) +
    theme
  p
}

# Prepare data (mean across genes)
prepare.data <- function(g){
  df <- lapply(donorNames, function(d) {
    expr <- brainExpr[[d]][g, ]
    expr <- apply(expr, 2, mean) # Mean across genes
    label <- braakLabels[[d]]
    # donor <- rep(d, length(expr))
    data.frame(expr, label)#, donor)
  })
  df <- melt(df, value.name = "expr")
  colnames(df)[4] <- "donor"
  df <- df[df$label != "0", ]#Remove Braak 0
  df$label <- paste0("R", df$label)
  df$label <- factor(df$label, levels = sort(unique(df$label)))
  df$donor <- factor(df$donor, levels = unique(df$donor))
  df
}

# Boxplot for a gene
boxplot.gene <- function(g, title){
  df <- prepare.data(g)
  box.plot(df, title)
}

plot.pdf <- function(name, genes){
  pdf(name, 6, 4)
  lapply(genes, function(g){
    r <- format(summaryLabelCor[[g]]["summary", c("r", "pvalue")], digits = 2)
    title <- paste0(entrezId2Name(g), ", r=", r$r)
    df <- prepare.data(g)
    p <- box.plot(df, title)
    print(p)
  })
  dev.off()
}

# Boxplot for each PD-implicated gene
plot.pdf("boxplot_high_impact_genes.pdf", pdGenesID$hiImpact)
plot.pdf("boxplot_susceptible_genes.pdf", pdGenesID$jansen2017)
plot.pdf("boxplot_HLA_genes.pdf", pdGenesID$hla)

# Boxplot for mean expression of -ve and +ve Braak genes
bg <- list(
  down = braakGenes$entrez_id[braakGenes$r < 0],
  up = braakGenes$entrez_id[braakGenes$r > 0]
)
meanExpr <- lapply(bg, prepare.data) 
df <- melt(meanExpr)
colnames(df) <- c("label", "variable", "donor", "expr", "dir")
y_max <- max(sapply(meanExpr, function(x) max(x$expr)))
y_min <- min(sapply(meanExpr, function(x) min(x$expr)))

pdf("boxplot_AHBA.pdf", 6, 4)
box.plot(df, "Mean BRGs") + facet_grid(.~dir, space = "free", scales = "free") +
  scale_y_continuous(limits = c(y_min, y_max))
dev.off()
