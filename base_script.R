#Base script

options(stringsAsFactors = FALSE)

probeInfo <- read.csv("../ABA_human_processed/probe_info_2014-11-11.csv")

entrezId2Name <- function (x) {probeInfo$gene_symbol[match(x, probeInfo$entrez_id)]} #Input is vector

name2EntrezId <- function (x) {as.character(probeInfo$entrez_id[match(x, probeInfo$gene_symbol)])} #Input is vector


donorNames <- c("donor9861", "donor10021", "donor12876", "donor14380", "donor15496", "donor15697")
names(donorNames) <- donorNames

braakNames <- c("braak1","braak2","braak3","braak4","braak5","braak6", "braak1-2", "braak3-4", "braak5-6", "braak1-3", "braak4-6")
names(braakNames) <- braakNames

# PD-implicated genes
pdGenes <- list(hiImpact = c("GBA", "LRRK2", "PINK1", "PARK7", "SNCA", "VPS35", "DNAJC13", "CHCHD2"),
                susceptible = c("INPP5F", "TMEM175", "ASH1L", "MAPT", "RIT1", "C14orf83", "STK39", "GPNMB", "BST1", 
                                "SIPA1L2", "DLG2", "NUCKS1", "GCH1", "MCCC1", "FAM47E", "BCKDK", "TMPRSS9", "UBOX5", 
                                "CCDC62", "SYNJ1", "EIF4G1", "FBXO7", "C20orf30", "POLG", "VSP13C", "PLA2G6"),
                hla = c("HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DQB1"),
                lysosome = unlist(read.table("lysosome_geneset.txt", header = FALSE, comment.char = "#", sep = "\n", row.names = NULL))
)

pdGenesID <- lapply(pdGenes, function(l){
  geneInfo <- probeInfo[probeInfo$gene_symbol %in% l, ]
  as.character(geneInfo[ , 6])
})
