# Convert Entrez IDs to Ensembl IDs of all AHBA genes
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(biomaRt)

ahba <- ahba.genes()
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version = 92)

# Get ensemble ID and HGNC gene symbol
id_conversion <- getBM(c('ensembl_gene_id', 'entrezgene', 'hgnc_symbol'), filters=c('entrezgene'), mart=ensembl, values=ahba)
id_conversion$ahba_symbol <- entrezId2Name(id_conversion$entrezgene)
id_conversion <- id_conversion[order(id_conversion$entrezgene), ]
write.table(id_conversion, file = "ahba_entrez2ensembl.txt", quote = FALSE, row.names = FALSE, sep = "\t")

# load("resources/braakGenes.RData")
# rows <- which(id_conversion$entrezgene %in% braakGenes$entrez_id)
# id_conversion_pd <- id_conversion[rows, ]
# id_conversion_pd$r <- braakGenes$r[match(id_conversion_pd$entrezgene, braakGenes$entrez_id)]
# write.table(id_conversion_pd, file = "ahba_entrez2ensembl_braak.txt", quote = FALSE, row.names = FALSE, sep = "\t")