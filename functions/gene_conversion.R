# Functions for gene conversino based on gene annotations from AHBA

entrezId2Name <- function (x) {probeInfo$gene_symbol[match(x, probeInfo$entrez_id)]} #Input is vector
name2EntrezId <- function (x) {as.character(probeInfo$entrez_id[match(x, probeInfo$gene_symbol)])} #Input is vector
