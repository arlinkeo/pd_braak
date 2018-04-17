# Differential expression in GTEx data base

setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")

df <- read.table("../GTEX/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct", sep = "\t", comment.char = "#", header = TRUE, skip = 2)
rownames(df) <- df$Name
# load("../GTEX/gtex_expr_norm.RData")
# rownames(gtex_expr_norm) <- gtex_expr_norm$Name

annot <- read.table("../GTEX/GTEx_v7_Annotations_SampleAttributesDS.txt", sep = "\t", header = TRUE, fill = T, quote = "")

roi <- c("Substantia nigra", "Amygdala", "Anterior cingulate cortex", "Frontal Cortex")
names(roi) <- c("SN", "AMG", "ACC", "FC")
samples <- lapply(roi, function(r){
  rows <- grep(r, annot$SMTSD)
  ids <- gsub("-", ".", annot[rows, "SAMPID"])
  intersect(ids, colnames(gtex_expr_norm))
})

tab <- lapply(samples, function(s){
  gtex_expr_norm[, s]
})

# T-test
genes <-gtex_expr_norm$Name

# T-test function for 2 regions
ttest.regions <- function(a,b){
  as.data.frame(t(sapply(genes, function(g){
    print(paste0(a, "-", b, "; ", which(genes %in% g), "; ", g))
    r1 <- unlist(tab[[a]][g, ])
    r2 <- unlist(tab[[b]][g, ])
    t <- t.test(r1, r2)
    meanDiff <- t$estimate[1]-t$estimate[2]
    pvalue <- t$p.value
    c(meanDiff=unname(meanDiff), pvalue=unname(pvalue))
  })))
}

regionpairs <- combn(names(roi), 2)
colnames(regionpairs) <- apply(regionpairs, 2, function(x) paste0(x[1], "-",  x[2]))
ttest.all <- apply(regionpairs, 2, function(x){
  r1 <- x[1]
  r2 <- x[2]
  ttest <- ttest.regions(r1, r2)
  file <- paste0("Resources/ttest_gtex_", r1, "_", r2, ".RData")
  print(file)
  save(ttest, file = file)
  ttest
})
# names(ttest.all) <- colnames(regionpairs)
save(ttest.all, file = "Resources/ttest.all_gtex_norm.RData")

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
id_conversion <- lapply(braakGenes, function(x){
  getBM(c('ensembl_gene_id', 'entrezgene', 'hgnc_symbol'), filters=c('entrezgene'), mart=ensembl, values=x)
})
braak_gtex <- sapply(diffGenes, function(rp){
  sapply(id_conversion, function(x){
    ens <- x$ensembl_gene_id
    overlap <- sapply(ens, function(id){
      sum(grepl(id, row.names(rp)))
    })
    sum(overlap)
  })
})
  
  


