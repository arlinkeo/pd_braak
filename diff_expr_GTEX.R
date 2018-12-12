# Differential expression in GTEx data base
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(reshape2)
library(ggplot2)
library(plyr)
load("resources/braakInfo.RData") # Braak colors
source("PD/t.test.table.R")
load("resources/braakGenes.RData")
# load("resources/braakGenes2.RData")
# load("resources/braakGenes3.RData")

# Load and filter data

gtex_expr <- read.table("../GTEX/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct", sep = "\t", comment.char = "#", header = TRUE, skip = 2)
rownames(gtex_expr) <- gtex_expr$Name
# genes <- gtex_expr[, c(1,2)]
gtex_expr <- gtex_expr[ ,-c(1,2)]

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
genes <- sapply(genes, function(s) {# remove version from ID
  unlist(strsplit(s, split = "\\."))[1]
})
rownames(gtex_expr) <- genes

# Filter out genes with low TPM values
maxTPM <- apply(gtex_expr, 1, max)
gtex_expr <- gtex_expr[which(maxTPM >= 1), ]
# saveRDS(gtex_expr, file = "resources/gtex_expr.rds")
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

# remove genes with 0's in all selected regions
zerorows <- lapply(samples, function(s){
  apply(gtex_expr[, s], 1, sum) == 0
})
zerorows <- Reduce("|", zerorows)
gtex_expr <- gtex_expr[!zerorows, ]

# T-test
regionpairs <- combn(names(roi), 2)
colnames(regionpairs) <- apply(regionpairs, 2, function(x) paste0(x[1], "-",  x[2]))
ttest <- alply(regionpairs, 2, function(x){
  s1 <- samples[[x[1]]]
  s2 <- samples[[x[2]]]
  df1 <- gtex_expr[, s1]
  df2 <- gtex_expr[, s2]
  t.test.table(df1,df2)
}, .dims = TRUE)
ttest <- simplify2array(ttest) # 3D array: genes x measures x region pairs
save(ttest, file = "Resources/ttest_gtex.RData")

# Number of diff. genes
apply(ttest, c(3), function(x){
  sum(x[, "BH"] < 0.05 & abs(x[, "FC"]) > 1)
})

# T-test Overlap with Braak genes

bg <- list(
  bg1 = list( # Braak genes selected WIHTOUT cell-type correction
    down = braakGenes$entrez_id[braakGenes$r < 0],
    up = braakGenes$entrez_id[braakGenes$r > 0]
  )#,
  # bg2 = list( # Braak genes selected WITH cell-type correction
  #   down = braakGenes2$entrez_id[braakGenes2$braak6 < 0],
  #   up = braakGenes2$entrez_id[braakGenes2$braak6 > 0]
  # ),
  # bg3 = braakGenes3 # Intersection of corrected and uncorrected results
)

# Convert entrez IDs to ensembl IDs
conversion_tab <- read.csv("ahba_entrez2ensembl.txt", sep = "\t")
bg <- lapply(bg, function(s){
  lapply(s, function(g){
    rows <- which(conversion_tab$entrezgene %in% g)
    id <- conversion_tab[rows, "ensembl_gene_id"]
    intersect(id, rownames(gtex_expr))
  })
})

apply(ttest, c(3), function(x){
  sapply(bg$bg1, function(g){
    sum(x[g, "BH"] < 0.05 & abs(x[g, "FC"]) > 1)
  })
})

# Box plot of Braak stage related genes

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
mean_df <- lapply(bg, function(s){
  t <- lapply(samples, function(n){
    sapply(s, function(g){
      apply(gtex_expr[g, n], 2, mean)
    })
  })
  df <- melt(t)
  colnames(df) <- c("sample", "dir", "expr", "region")
  df$region <- paste0("R", df$region)
  df$region <- factor(df$region, levels = sort(unique(df$region)))
  df
})

y_max <- max(sapply(mean_df, function(x) max(x$expr)))

pdf("boxplot_GTEX.pdf", 4, 3)
lapply(names(mean_df), function(n){
  df <- mean_df[[n]]
  box.plot(df, n) +
    facet_grid(.~dir, scales = 'free', space = 'free', switch = "y") +
    scale_y_continuous(limits = c(0, y_max))
})
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
