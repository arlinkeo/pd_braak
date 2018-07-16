# Boxplot of expression in Braak regions for a single gene
setwd("C:/Users/dkeo/surfdrive/pd_braak")
library(ggplot2)
source("PD/base_script.R")
load("../ABA_Rdata/BrainExpr.RData")
load("resources/braakInfo.RData") # Braak stage label vectors
load("resources/summaryLabelCorr.RData")
load("resources/braakGenes.RData")

# Default theme for boxplot
theme <- theme(panel.background = element_blank(), panel.grid = element_blank(), 
               axis.line = element_line(colour = "black"),
               legend.title = element_blank(),
               legend.key = element_blank())

# Boxplot function
box.plot <- function(df, title){
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

# Prepare data
prepare.data <- function(g){
  df <- lapply(donorNames, function(d) {
    expr <- brainExpr[[d]][g, ]
    expr <- apply(expr, 2, mean) # Mean across genes
    label <- braakLabels[[d]]
    donor <- rep(d, length(expr))
    data.frame(expr, label, donor)
  })
  df <- Reduce(rbind, df)
  df <- df[df$label != "0", ]#Remove Braak 0
  df$label <- factor(df$label, levels = sort(unique(df$label)))
  df$donor <- factor(df$donor, levels = unique(df$donor))
  df
}

plot.pdf <- function(name, genes){
  pdf(name, 8, 5)
  lapply(genes, function(g){
    r <- format(summaryLabelCorr[[g]]["summary", c("r", "pvalue")], digits = 2)
    title <- paste0(entrezId2Name(g), ", r=", r$r)
    df <- prepare.data(g)
    p <- box.plot(df, title)
    print(p)
  })
  dev.off()
}

# Boxplot for each PD-implicated gene
plot.pdf("boxplot_high_impact_genes.pdf", pdGenesID$hiImpact)
plot.pdf("boxplot_susceptible_genes.pdf", pdGenesID$susceptible)
plot.pdf("boxplot_HLA_genes.pdf", pdGenesID$hla)

# Boxplot for mean expression of -ve and +ve Braak genes
braak <- list(
  'r<0' = braak_neg <- braakGenes$entrez_id[braakGenes$braak_r < 0],
  'r>0' = braakGenes$entrez_id[braakGenes$braak_r > 0])

df <- lapply(names(braak), function(r){
  g <- braak[[r]]
  title <- r
  df <- prepare.data(g)
  df$r <- r
  df
})
df <- Reduce(rbind, df)

pdf("boxplot_AHBA.pdf", 10, 4)
box.plot(df, "Braak genes") + facet_grid(.~r, space = "free", scales = "free")
dev.off()