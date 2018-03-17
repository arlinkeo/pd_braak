#Base script

options(stringsAsFactors = FALSE)

##### AHBA data and functions #####

# Function for gene conversion
probeInfo <- read.csv("../ABA_human_processed/probe_info_2014-11-11.csv")
entrezId2Name <- function (x) {probeInfo$gene_symbol[match(x, probeInfo$entrez_id)]} #Input is vector
name2EntrezId <- function (x) {as.character(probeInfo$entrez_id[match(x, probeInfo$gene_symbol)])} #Input is vector

# Brain donor names
donorNames <- c("donor9861", "donor10021", "donor12876", "donor14380", "donor15496", "donor15697")
names(donorNames) <- donorNames

# Function to get  AHBA genes
ahba.genes <- function(random = NULL){
  genes <- probeInfo$entrez_id
  genes <- if (!is.null(random)) sample(genes, random) else genes
  as.character(genes)
}

# Function to get list of subselected expression matrices for each donor
select.expr <- function(genes = NULL, samples = NULL){
  expr <- readRDS("../ABA_Rdata/BrainExprNorm.rds")
  lapply(donorNames, function(d){
    e <- expr[[d]]
    if (!is.null(samples)){
      labels <- as.logical(samples[[d]])
      if (is.null(genes)) e[, labels]
      else e[genes, labels]
    } else {
      if (is.null(genes)) e
      else e[genes, ]
    }
  })
}

##### PD info about brain regions and genes #####

braakNames <- c("braak1","braak2","braak3","braak4","braak5","braak6")
names(braakNames) <- braakNames
braakNamesMerged1 <- c("braak1-2", "braak3-4", "braak5-6")
names(braakNamesMerged1) <- braakNamesMerged1
braakNamesMerged2 <- c("braak1-3", "braak4-6")
names(braakNamesMerged2) <- braakNamesMerged2

pdGenes <- list(hiImpact = c("GBA", "LRRK2", "PINK1", "PARK7", "SNCA", "VPS35", "ATP13A2", "PLA2G6", "FBXO7", "DNAJC6", "SYNJ1", "EIF4G1", "DNAJC13", "CHCHD2", "TMEM230", "RIC3"),
                susceptible = c("INPP5F", "TMEM175", "ASH1L", "MAPT", "RIT1", "C14orf83", "STK39", "GPNMB", "BST1", 
                                "SIPA1L2", "DLG2", "NUCKS1", "GCH1", "MCCC1", "FAM47E", "BCKDK", "TMPRSS9", "UBOX5", 
                                "CCDC62", "SYNJ1", "EIF4G1", "FBXO7", "C20orf30", "POLG", "VPS13C", "PLA2G6"),
                hla = c("HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DQB1"),
                lysosome = read.table("lysosome_geneset.txt", header = FALSE, comment.char = "#", sep = "\n", row.names = NULL)[, 1],
                chang2017 = read.table("chang2017_riskgenes.txt", comment.char = "#", sep = "\n", row.names = NULL)[, 1] 
)

pdGenesID <- lapply(pdGenes, function(l){
  geneInfo <- probeInfo[probeInfo$gene_symbol %in% l, ]
  as.character(geneInfo[ , 6])
})


##### Basic functions #####

# Back transform fischer z-score correlations
back.transform <- dget("PD/back.transform.R")
gene.coexpr <- dget("PD/gene.coexpr.R")