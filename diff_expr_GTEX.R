# Differential expression in GTEx data base

setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(reshape)
library(ggplot2)

df <- read.table("../GTEX/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct", sep = "\t", comment.char = "#", header = TRUE, skip = 2)
rownames(df) <- df$Name

# Filter IDs for protein coding genes only
gene.annot <- read.table("../GTEX/gencode.v19.genes.v7.patched_contigs.gtf", sep = "\t")
gene.annot <- gene.annot[which(gene.annot$V3 == "gene"), ]
gene_info <- strsplit(gene.annot$V9, split = ";")
gene_type <- sapply(gene_info,  function(x){
  gt <- x[[grep("gene_type", x)]]
  gsub(" gene_type ", "", gt)
})
gene_id <- sapply(gene_info,  function(x){
  id <- x[[grep("gene_id", x)]]
  gsub("gene_id ", "", id)
})
genes <- gene_id[gene_type == "protein_coding"]

df <- df[genes, ]
genes <- sapply(genes, function(s) {# remove version from ID
  unlist(strsplit(s, split = "\\."))[1]
})
names(genes) <- genes
rownames(df) <- genes

# Select samples corresponding to brain tissues of Braak regions
sample.annot <- read.table("../GTEX/GTEx_v7_Annotations_SampleAttributesDS.txt", sep = "\t", header = TRUE, fill = T, quote = "")

roi <- c("Substantia nigra", "Amygdala", "Anterior cingulate cortex", "Frontal Cortex")
names(roi) <- c("SN", "AMG", "ACC", "FC")
samples <- lapply(roi, function(r){
  rows <- grep(r, sample.annot$SMTSD)
  ids <- gsub("-", ".", sample.annot[rows, "SAMPID"])
  intersect(ids, colnames(df))
})

# select samples for each region
dfList <- lapply(samples, function(s){ 
  df[, s]
})

# T-test function for 2 regions
ttest.regions <- function(a,b){
  as.data.frame(t(sapply(genes, function(g){
    print(paste0(a, "-", b, "; ", which(genes %in% g), "; ", g))
    r1 <- unlist(dfList[[a]][g, ])
    r2 <- unlist(dfList[[b]][g, ])
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
load("resources/braakGenes.RData")
library(biomaRt)
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=91)
id_conversion <- lapply(braakGenes, function(x){ # Get ensemble ID and HGNC gene symbol
  getBM(c('ensembl_gene_id', 'entrezgene', 'hgnc_symbol'), filters=c('entrezgene'), mart=ensembl, values=x)
})
id_conversion <- lapply(id_conversion, function(x){ # Add AHBA gene symbols
  x$AHBAsymbol <- entrezId2Name(x$entrezgene)
  x
})

braak_ensembl <- lapply(id_conversion, function(x){ # For + and - correlated progression genes
  ens <- x$ensembl_gene_id
  intersect(ens, genes)
})

# Correct p-value for number of progression genes
braak_gtex <- lapply(ttest.all, function(rp){
  sapply(names(braak_ensembl), function(r){ # For + and - correlated progression genes
    rows <- braak_ensembl[[r]]
    t <- rp[rows,]
    t$BH <- p.adjust(t$pvalue)
    if (r== "positive_r") t[t$meanDiff < 2 & t$BH < 0.05, ]
    else t[t$meanDiff > 2 & t$BH < 0.05, ]
  }, simplify = FALSE)
})

sapply(braak_gtex, function(rp)sapply(rp,  nrow))
gtex_pos <- lapply(braak_gtex, function(x){
  rownames(x$positive_r)
})
gtex_neg <- lapply(braak_gtex, function(x){
  rownames(x$negative_r)
})
length(Reduce(union, gtex_pos[1:3]))
length(Reduce(union, gtex_neg[c(1,3)]))

# Box plot of all Braak genes
meanExpr <- lapply(braak_ensembl, function(r){
  t <- sapply(dfList, function(df){
    apply(df[r, ], 1, mean) # mean expression of progression genes
  })
  melt(t)
})
meanExpr <- melt.list(meanExpr)
colnames(meanExpr) <- c("gene", "region", "value", "mean_expr", "r")
meanExpr$region <- factor(meanExpr$region, levels = names(roi))

theme <- theme(panel.background = element_blank(), panel.grid = element_blank(), 
               axis.line = element_line(colour = "black"),
               legend.title = element_blank())

p <- ggplot(meanExpr) + 
  geom_boxplot(aes(y = mean_expr, x = region, fill = region)) +
  labs(x = "Braak stage region", y = "Expression (TPM)") +
  ggtitle("Expression of Braak-related genes in GTEx") +
  scale_x_discrete(expand=c(0.2,0), labels = c(3:6)) +
  facet_grid(.~r, scales = 'free', space = 'free', switch = "y") +
  theme
p

#boxplot of SNCA
ens_snca <- id_conversion$positive_r$ensembl_gene_id[which(id_conversion$positive_r$entrezgene== name2EntrezId("SNCA"))]
expr_snca <- lapply(dfList, function(x) unlist(x[ens_snca,]))

sapply(expr_snca, median)

expr_snca <- melt.list(expr_snca)
colnames(expr_snca) <- c("expr", "region")
expr_snca$region <- factor(expr_snca$region, levels = names(roi))

p <- ggplot(expr_snca) + 
  geom_boxplot(aes(y = expr, x = region, fill = region)) +
  labs(x = "Braak stage region", y = "Expression (TPM)") +
  ggtitle("Expression of SNCA in GTEx") +
  scale_x_discrete(expand=c(0.2,0), labels = c(3:6)) +
  theme
pdf("snca_gtex.pdf", 4, 4)
print(p)
dev.off()
