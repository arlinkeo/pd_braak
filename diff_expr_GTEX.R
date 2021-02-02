# Differential expression in GTEx data base
# First download the RNA-seq data and annotations (V7) from https://gtexportal.org/home/datasets
# Files: 
# GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct.gz
# GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz
# GTEx_v7_Annotations_SubjectPhenotypesDS.txt
# GTEx_v7_Annotations_SampleAttributesDS.txt

library(DESeq2)
library(ggsignif)

# Load and filter data
gtex_expr <- read.table("../../GTEX/Gene_read_counts/All_Tissue_Site_Details.combined.reads.gct", sep = "\t", comment.char = "#", header = TRUE, skip = 2)
rownames(gtex_expr) <- gtex_expr$Name
gtex_expr <- gtex_expr[ ,-c(1,2)]

# Subject info
subject_info <- read.table("../../GTEX/GTEx_v7_Annotations_SubjectPhenotypesDS.txt", sep = "\t", header = TRUE)
total <- nrow(subject_info)
male <- (sum(subject_info$SEX == "1") / total)*100

# Filter ensembl IDs for protein coding genes only
gene.annot <- read.table("../../GTEX/gencode.v19.genes.v7.patched_contigs.gtf", sep = "\t")
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

saveRDS(gtex_expr, file = "output/gtex_expr.rds")

################################################################################

# Select samples corresponding to brain tissues of Braak regions
sample.annot <- read.table("../../GTEX/GTEx_v7_Annotations_SampleAttributesDS.txt", sep = "\t", header = TRUE, fill = T, quote = "")
roi <- c('3' = "Substantia nigra", '4' = "Amygdala", '5' = "Anterior cingulate cortex", '6' = "Frontal Cortex")
sample_list <- lapply(roi, function(r){
  rows <- grep(r, sample.annot$SMTSD)
  ids <- gsub("-", ".", sample.annot[rows, "SAMPID"])
  intersect(ids, colnames(gtex_expr))
})

zerorows <- lapply(sample_list, function(s){
  apply(gtex_expr[, s], 1, sum) == 0
})
zerorows <- Reduce("|", zerorows)
gtex_expr <- gtex_expr[!zerorows, ]

samples <- melt(sample_list)
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
saveRDS(deseq, file = "output/deseq_gtex.rds")

# Number of diff. genes
apply(deseq, c(3), function(x){
  sum(x[, "padj"] < 0.05 & abs(x[, "log2FoldChange"]) > 1, na.rm = TRUE)
})

# Test Overlap with Braak genes

# Convert entrez IDs to ensembl IDs and intersect with GTEx genes
# id_conversion <- read.csv("output/ahba_entrez2ensembl.txt", sep = "\t")
bg_gtex <- lapply(bg, function(g){
  rows <- which(id_conversion$entrezgene %in% g)
  id <- id_conversion[rows, "ensembl_gene_id"]
  intersect(id, rownames(gtex_expr))
})

# Number of diff. Braak genes
degs_num_gtex <- apply(deseq, c(3), function(x){
  sapply(bg_gtex, function(g){
    sum(x[g, "padj"] < 0.05 & abs(x[g, "log2FoldChange"]) > 1, na.rm = TRUE)
  })
})
degs_percentage_gtex <- apply(degs_num_gtex, 2, function(x) round(x/sapply(bg_gtex,length)*100, digits = 2))
degs_num_gtex
degs_percentage_gtex

########## Box plot of Braak stage related genes ##########

# Load TPM values
gtex_tpm <- read.table("../../GTEX/Gene_TPMs/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct", sep = "\t", comment.char = "#", header = TRUE, skip = 2)
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
    geom_boxplot(aes(x = region, y = expr, fill = region)) +
    labs(x = "Brain region", y = "Expression (TPM)") +
    ggtitle(title) +
    scale_x_discrete(expand=c(0.2,0)) +
    scale_fill_manual(values = braakColors, guide = FALSE) +
    theme
}

# Mean expression across Braak genes
df <- lapply(bg_gtex, function(g){
  t <- apply(gtex_tpm[g, ], 2, mean)
  t <- data.frame(sample = names(t), expr = t, region = samples$region)
  # melt(t)
})
df <- melt(df)
colnames(df) <- c("sample", "region", "variable", "expr", "dir")
df$region <- paste0("R", df$region)
df$region <- factor(df$region, levels = sort(unique(df$region)))

pdf("output/boxplot_GTEX.pdf", 2.5, 4)
box.plot(df, "BRGs in GTEx") +
  facet_grid(.~dir, scales = 'free', space = 'free', switch = "y")
dev.off()

########## Heatmap expression of BRGs ##########

expr <- sapply(sample_list, function(s){
  e <- gtex_tpm[unlist(bg_gtex), s]
  apply(e, 1, mean)
})
colnames(expr) <- paste0("R", colnames(expr))
expr <- t(scale(t(expr))) # expr. is scaled across samples

pdf("output/heatmap_expr_BRGs_GTEx_col1.2.pdf", 2.7, 10)
Heatmap(expr, name = 'Z-Score\nexpression',
        col = col_fun,
        row_split = rep(names(lengths(bg_gtex)), lengths(bg_gtex)),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        column_names_gp = gpar(fontsize = 10),
        column_names_rot = 0,
        column_names_centered = TRUE,
        row_title_rot = 0,
        width = unit(ncol(expr), "lines"),
        height = unit(nrow(expr)*.05, "lines")
)
dev.off()

# #boxplot of PD genes
# prepare.data <- function(g){
#   df <- unlist(gtex_tpm[g, ])
#   df <- data.frame(sample = names(df), expr = df, region = samples$region)
#   df <- melt(df)
#   colnames(df) <- c("sample", "region", "variable", "expr")
#   df$region <- paste0("R", df$region)
#   df$region <- factor(df$region, levels = sort(unique(df$region)))
#   df
# }
# 
# plot.pdf <- function(name, genes){
#   pdf(name, 2, 3)
#   lapply(genes, function(g){
#     title <- paste0(entrezId2Name(g))
#     g <- conversion_tab$ensembl_gene_id[which(conversion_tab$entrezgene == g)]
#     df <- prepare.data(g)
#     p <- box.plot(df, title)# +
#     #   geom_signif(comparisons = list(c("R3", "R4")), 
#     #               step_increase = 0.1)
#     # p
#     # pval <- deseq[g, "padj", ]
#     # pg <- ggplot_build(p)
#     # pg$data[[1]]$annotation <- pval
#     # pg$data[[2]]$textsize <- 2.5
#     # pg$data[[2]]$colour <- ifelse(as.numeric(pg$data[[2]]$annotation) < 0.05, "red", "black")
#     # q1 <- ggplot_gtable((pg))
#     # p1 <- plot(q1)
#     print(p)
#   })
#   dev.off()
# }
# 
# plot.pdf("boxplot_GTEX_PD_variant_genes.pdf", 
#          name2EntrezId(c("DNAJC13", "SNCA", "GCH1", "INPP5F", "ASH1L", "ZNF184", "DDRGK1", "ITPKB", "ELOVL7", "SCARB2")))