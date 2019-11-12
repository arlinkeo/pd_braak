# Convert Entrez IDs to Ensembl IDs of all AHBA genes
# Conversion of IDs is used for analysis of PD RNA-seq dataset and GTEx RNAseq dataset
library(biomaRt)

ahba <- ahba.genes()
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version = 92)

# Get ensemble ID and HGNC gene symbol
id_conversion <- getBM(c('ensembl_gene_id', 'entrezgene', 'hgnc_symbol'), filters=c('entrezgene'), mart=ensembl, values=ahba)
id_conversion$ahba_symbol <- entrezId2Name(id_conversion$entrezgene)
id_conversion <- id_conversion[order(id_conversion$entrezgene), ]
write.table(id_conversion, file = "output/ahba_entrez2ensembl.txt", quote = FALSE, row.names = FALSE, sep = "\t")