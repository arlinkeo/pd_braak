#Base script

options(stringsAsFactors = FALSE)

probeInfo <- read.csv("../ABA_human_processed/probe_info_2014-11-11.csv")

entrezId2Name <- function (x) {probeInfo$gene_symbol[match(x, probeInfo$entrez_id)]} #Input is vector

name2EntrezId <- function (x) {as.character(probeInfo$entrez_id[match(x, probeInfo$gene_symbol)])} #Input is vector


donorNames <- c("donor9861", "donor10021", "donor12876", "donor14380", "donor15496", "donor15697")
names(donorNames) <- donorNames