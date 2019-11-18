setwd("M:/doorgeefluik/Arlin Keo/Arrays data_coded for Arlin")
options(stringsAsFactors = FALSE)
library("WGCNA")
library(hgu133plus2.db)

expr <- read.csv("rma_CODED_no commas_correct order.csv")
rownames(expr) <- unlist(expr[,1])
expr <- expr[, -1]

# Link probes to entrez IDs
probes <- rownames(expr)
entrezIDs <- mapIds(hgu133plus2.db, probes, c("ENTREZID"), keytype = "PROBEID")
na_rows <- which(is.na(entrezIDs))

# Remove NA rows
entrezIDs <- entrezIDs[-na_rows]
probes <- probes[-na_rows]
expr <- expr[-na_rows, ]

# Map probes to genes
expr2 <- collapseRows(expr, entrezIDs, probes, method = "maxRowVariance", connectivityBasedCollapsing = TRUE)
save(expr2, file = "expr2.RData")

expr_new <- expr2$datETcollapsed
expr_new <- expr_new[order(as.numeric(rownames(expr_new))), ] # sort rows by numeric entrez ID
write.csv(expr_new, file = "rma_genes.csv", quote = FALSE)