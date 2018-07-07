# Differential expression in GTEx data base

setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(reshape)
library(ggplot2)
load("resources/braakInfo.RData") # Braak colors

df <- readRDS("../GTEX/gtex_expr.rds")
# df <- readRDS("../GTEX/gtex_expr_norm.rds")

# Select samples corresponding to brain tissues of Braak regions
sample.annot <- read.table("../GTEX/GTEx_v7_Annotations_SampleAttributesDS.txt", sep = "\t", header = TRUE, fill = T, quote = "")
roi <- c('3' = "Substantia nigra", '4' = "Amygdala", '5' = "Anterior cingulate cortex", '6' = "Frontal Cortex")
samples <- lapply(roi, function(r){
  rows <- grep(r, sample.annot$SMTSD)
  ids <- gsub("-", ".", sample.annot[rows, "SAMPID"])
  intersect(ids, colnames(df))
})

# remove genes with 0's in all selected regions
zerorows <- lapply(samples, function(s){
  apply(df[, s], 1, sum) == 0
})
zerorows <- Reduce("|", zerorows)

df <- df[!zerorows, ]
genes <- rownames(df)

# T-test function for 2 regions
ttest.regions <- function(a,b){
  s1 <- samples[[a]]
  s2 <- samples[[b]]
  df1 <- df[, s1]
  df2 <- df[, s2]
  as.data.frame(t(sapply(genes, function(g){
    print(paste0(a, "-", b, "; ", which(genes %in% g), "; ", g))
    r1 <- df1[g, ]
    r2 <- df2[g, ]
    t <- t.test(r1, r2)
    meanDiff <- t$estimate[1]-t$estimate[2]
    pvalue <- t$p.value
    c(meanDiff=unname(meanDiff), pvalue=unname(pvalue))
  })))
}

# T-test
regionpairs <- combn(names(roi), 2)
colnames(regionpairs) <- apply(regionpairs, 2, function(x) paste0(x[1], "-",  x[2]))
ttest.all <- apply(regionpairs, 2, function(x){
  r1 <- x[1]
  r2 <- x[2]
  ttest.regions(r1, r2)
})
save(ttest.all, file = "Resources/ttest.all_gtex.RData")

# Number of diff. genes
diffGenes <- lapply(ttest.all, function(rp){
  rp$BH <- p.adjust(rp$pvalue)
  rp[abs(rp$meanDiff) > 2 & rp$BH < 0.01, ]
})
sapply(diffGenes, nrow)

# Overlap with progression genes
# load("resources/braakGenes.RData")

# Correct p-value for number of all genes
# braak_gtex <- lapply(ttest.all, function(rp){
#   rp$BH <- p.adjust(rp$pvalue)
#   sapply(names(braak_ensembl), function(r){ # For + and - correlated progression genes
#     rows <- braak_ensembl[[r]]
#     t <- rp[rows,]
#     if (r== "positive_r") t[t$meanDiff < 2 & t$BH < 0.05, ]
#     else t[t$meanDiff > 2 & t$BH < 0.05, ]
#   }, simplify = FALSE)
# })
# 
# # Correct p-value for number of progression genes
# braak_gtex <- lapply(ttest.all, function(rp){
#   sapply(names(braak_ensembl), function(r){ # For + and - correlated progression genes
#     rows <- braak_ensembl[[r]]
#     t <- rp[rows,]
#     t$BH <- p.adjust(t$pvalue)
#     if (r== "positive_r") t[t$meanDiff < 2 & t$BH < 0.05, ]
#     else t[t$meanDiff > 2 & t$BH < 0.05, ]
#   }, simplify = FALSE)
# })
# 
# sapply(braak_gtex, function(rp)sapply(rp,  nrow))
# gtex_pos <- lapply(braak_gtex, function(x){
#   rownames(x$positive_r)
# })
# gtex_neg <- lapply(braak_gtex, function(x){
#   rownames(x$negative_r)
# })
# length(Reduce(union, gtex_pos[1:3]))
# length(Reduce(union, gtex_neg[c(1,3)]))

# Box plot of all Braak genes

# Braak genes
conversion_tab <- read.csv("ahba_entrez2ensembl_braak.txt", sep = "\t")
braak_pos <- conversion_tab$ensembl_gene_id[conversion_tab$dir == "pos"]
braak_neg <- conversion_tab$ensembl_gene_id[conversion_tab$dir == "neg"]
braak_pos <- intersect(rownames(df), braak_pos)
braak_neg <- intersect(rownames(df), braak_neg)

# Mean expression across Braak genes
meanExpr <- lapply(samples, function(s){
  t <- sapply(list("r<0" = braak_neg, "r>0" = braak_pos), function(g){
    apply(df[g, s], 2, mean)
  })
  t=melt(t)
})
meanExpr <- melt.list(meanExpr)
colnames(meanExpr) <- c("sample", "r", "variable", "expr", "region")
meanExpr$region <- factor(meanExpr$region, levels = unique(meanExpr$region))

# meanExpr <- lapply(braak_ensembl, function(r){
#   t <- sapply(samples, function(s){
#     d <- df[, s]
#     apply(d[r, ], 1, mean) # mean expression of progression genes
#   })
#   melt(t)
# })
# meanExpr <- melt.list(meanExpr)
# colnames(meanExpr) <- c("gene", "region", "value", "mean_expr", "r")
# meanExpr$region <- factor(meanExpr$region, levels = names(roi))

theme <- theme(panel.background = element_blank(), panel.grid = element_blank(), 
               axis.line = element_line(colour = "black"))

names(braakColors) <- gsub("braak", "", names(braakColors))

box.plot <- function(df, title){
  ggplot(df) + 
    geom_boxplot(aes(y = expr, x = region, fill = region)) +
    labs(x = "Braak stage", y = "Expression (TPM)") +
    ggtitle(title) +
    scale_x_discrete(expand=c(0.2,0)) +
    scale_fill_manual(values = braakColors, guide = FALSE) +
    theme
}
p1 <- box.plot(meanExpr, "Expression of Braak genes in GTEx") +
  facet_grid(.~r, scales = 'free', space = 'free', switch = "y")

#boxplot of SNCA
ens_snca <- conversion_tab$ensembl_gene_id[which(conversion_tab$ahba_symbol == "SNCA")]
expr_snca <- lapply(samples, function(s) unlist(df[ens_snca,s]))

sapply(expr_snca, median)

expr_snca <- melt.list(expr_snca)
colnames(expr_snca) <- c("expr", "region")
expr_snca$region <- factor(expr_snca$region, levels = unique(expr_snca$region))

p2 <- box.plot(expr_snca, "Expression of SNCA in GTEx")
p2
pdf("boxplot_GTEX.pdf", 4, 3)
print(p1)
print(p2)
dev.off()