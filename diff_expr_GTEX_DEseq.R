# Differential expression in GTEx data base
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(reshape2)
library(ggplot2)
library(plyr)
library(DESeq2)
load("resources/braakInfo.RData") # Braak colors
source("PD/t.test.table.R")
load("resources/braakGenes.RData")

# Load and filter data
gtex_expr <- read.table("../GTEX/Gene_read_counts/All_Tissue_Site_Details.combined.reads.gct", sep = "\t", comment.char = "#", header = TRUE, skip = 2)
rownames(gtex_expr) <- gtex_expr$Name
gtex_expr <- gtex_expr[ ,-c(1,2)]

# Subject info
subject_info <- read.table("../GTEX/GTEx_v7_Annotations_SubjectPhenotypesDS.txt", sep = "\t", header = TRUE)
total <- nrow(subject_info)
male <- (sum(subject_info$SEX == "1") / total)*100

# Filter ensembl IDs for protein coding genes only
gene.annot <- read.table("../GTEX/gencode.v19.genes.v7.patched_contigs.gtf", sep = "\t")
gene.annot <- gene.annot[which(gene.annot$V3 == "gene"), ] # Filter for genes
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

# Expression of only protein coding genes
gtex_expr <- gtex_expr[genes, ]
genes2 <- sapply(genes, function(s) {# remove version from ID
  unlist(strsplit(s, split = "\\."))[1]
})
rownames(gtex_expr) <- genes2

saveRDS(gtex_expr, file = "resources/gtex_expr.rds")
gtex_expr <- readRDS("resources/gtex_expr.rds")

################################################################################

# Select samples corresponding to brain tissues of Braak regions
sample.annot <- read.table("../GTEX/GTEx_v7_Annotations_SampleAttributesDS.txt", sep = "\t", header = TRUE, fill = T, quote = "")
roi <- c('3' = "Substantia nigra", '4' = "Amygdala", '5' = "Anterior cingulate cortex", '6' = "Frontal Cortex")
samples <- lapply(roi, function(r){
  rows <- grep(r, sample.annot$SMTSD)
  ids <- gsub("-", ".", sample.annot[rows, "SAMPID"])
  intersect(ids, colnames(gtex_expr))
})

zerorows <- lapply(samples, function(s){
  apply(gtex_expr[, s], 1, sum) == 0
})
zerorows <- Reduce("|", zerorows)
gtex_expr <- gtex_expr[!zerorows, ]

samples <- melt(samples)
colnames(samples) <- c("sample", "region")
samples$region <- factor(samples$region, levels = unique(samples$region))

# DEseq
regionpairs <- combn(names(roi), 2)
colnames(regionpairs) <- apply(regionpairs, 2, function(x) paste0(x[1], "-",  x[2]))
regionpairs <- as.list(as.data.frame(apply(regionpairs, 2, function(x)c("region", x))))

gtex_expr <- gtex_expr[, samples$sample]
dds <- DESeqDataSetFromMatrix(countData = gtex_expr, samples, design = ~region)
deseq_test <- DESeq(dds)
deseq <- lapply(regionpairs, function(x){
  res <- results(deseq_test, contrast = x)
  t <- as.matrix(as.data.frame(res@listData))
  rownames(t) <- res@rownames
  t
})
deseq <- simplify2array(deseq) # 3D array: genes x measures x region pairs
save(deseq, file = "Resources/deseq_gtex.RData")

# Number of diff. genes
apply(deseq, c(3), function(x){
  sum(x[, "padj"] < 0.05 & abs(x[, "log2FoldChange"]) > 1, na.rm = TRUE)
})

# Test Overlap with Braak genes
bg <- list( # Braak genes selected WIHTOUT cell-type correction
  down = braakGenes$entrez_id[braakGenes$r < 0],
  up = braakGenes$entrez_id[braakGenes$r > 0]
)

# Convert entrez IDs to ensembl IDs
conversion_tab <- read.csv("ahba_entrez2ensembl.txt", sep = "\t")
bg <- lapply(bg, function(g){
  rows <- which(conversion_tab$entrezgene %in% g)
  id <- conversion_tab[rows, "ensembl_gene_id"]
  intersect(id, rownames(gtex_expr))
})

# Number of diff. Braak genes
apply(deseq, c(3), function(x){
  sapply(bg, function(g){
    sum(x[g, "padj"] < 0.05 & abs(x[g, "log2FoldChange"]) > 1, na.rm = TRUE)
  })
})

########## Box plot of Braak stage related genes ##########
####
# Plot TPM values

# Load TPM values
gtex_tpm <- read.table("../GTEX/Gene_TPMs/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct", sep = "\t", comment.char = "#", header = TRUE, skip = 2)
rownames(gtex_tpm) <- gtex_tpm$Name
gtex_tpm <- gtex_tpm[ ,-c(1,2)]
gtex_tpm <- gtex_tpm[genes, ]# Expression of only protein coding genes
rownames(gtex_tpm) <- genes2
gtex_tpm <- gtex_tpm[!zerorows, samples$sample]

theme <- theme(panel.background = element_blank(), 
               panel.grid = element_blank(), 
               axis.line = element_line(colour = "black"))

box.plot <- function(df, title){
  ggplot(df) + 
    geom_boxplot(aes(y = expr, x = region, fill = region)) +
    labs(x = "Brain region", y = "Expression (TPM)") +
    ggtitle(title) +
    scale_x_discrete(expand=c(0.2,0)) +
    scale_fill_manual(values = braakColors, guide = FALSE) +
    theme
}

# Mean expression across Braak genes
df <- lapply(bg, function(g){
  t <- apply(gtex_tpm[g, ], 2, mean)
  t <- data.frame(sample = names(t), expr = t, region = samples$region)
  melt(t)
})
df <- melt(df)
colnames(df) <- c("sample", "region", "variable", "expr", "dir")
df$region <- paste0("R", df$region)
df$region <- factor(df$region, levels = sort(unique(df$region)))
df$expr <- log2(df$expr)

pdf("boxplot_GTEX.pdf", 2.5, 4)
box.plot(df, "Mean expression across Braak genes") +
  facet_grid(.~dir, scales = 'free', space = 'free', switch = "y")
dev.off()

#boxplot of PD genes
prepare.data <- function(g){
  ens <- conversion_tab$ensembl_gene_id[which(conversion_tab$entrezgene == g)]
  df <- lapply(samples, function(s) unlist(gtex_expr[ens, s]))
  df <- melt(df)
  colnames(df) <- c("expr", "region")
  df$region <- paste0("R", df$region)
  df$region <- factor(df$region, levels = sort(unique(df$region)))
  # df$region <- factor(df$region, levels = unique(df$region))
  df
}

plot.pdf <- function(name, genes){
  pdf(name, 2, 3)
  lapply(genes, function(g){
    title <- paste0(entrezId2Name(g))
    df <- prepare.data(g)
    p <- box.plot(df, title)
    print(p)
  })
  dev.off()
}

plot.pdf("boxplot_GTEX_PD_variant_genes.pdf", 
         name2EntrezId(c("DNAJC13", "SNCA", "GCH1", "INPP5F", "ASH1L", "ZNF184", "DDRGK1", "ITPKB", "ELOVL7", "SCARB2")))