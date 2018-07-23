# Map probe to genes UKBEC dataset
setwd("/tudelft.net/staff-bulk/ewi/insy/DBL/Arlin/pd_braak")
library(biomaRt)
load("../UKBEC/expr.maps.rda",verbose=T)

# Map affy ID to entrez IDs
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version = 92)
affyID <- expr.map$exprID
system.time({
  ukbecGeneID <- getBM(c('ensembl_gene_id', 'entrezgene', 'hgnc_symbol', "affy_huex_1_0_st_v2"), 
                       filters=c("affy_huex_1_0_st_v2"), mart=ensembl, values=affyID)
})
save(ukbecGeneID, file = "resources/ukbecGeneID.RData")
write.table(ukbecGeneID, file = "ukbecGeneID.txt", row.names = FALSE, quote = FALSE, sep = "\t")