#Base script
options(stringsAsFactors = FALSE)

##### AHBA data and functions #####

# Function for gene conversion
probeInfo <- read.csv("../AHBA_Arlin/probe_info_2018-11-18.csv")
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

# IDs and AHBA colors for each sample per donor
ontology <- read.csv("../AHBA_Arlin/Ontology.csv")
sampleInfo <- lapply(donorNames, function(d){
  sampleIds <- read.csv(paste0("../AHBA_Arlin/sample_info_", d, "_2018-11-18.csv"))[ , 1]
  info <- ontology[match(sampleIds, ontology$id), ]
  # info$color_hex_triplet <- sapply(info$color_hex_triplet, function(c){
  #   if (nchar(c) == 5) {paste("#0", c, sep = "")} else {paste("#", c, sep = "")}
  # })
  info$color_hex_triplet <- paste0("#", info$color_hex_triplet)
  info
})

# Function to get list of subselected expression matrices for each donor
# select.expr <- function(genes = NULL, samples = NULL){
#   # expr <- readRDS("../ABA_Rdata/BrainExprNorm.rds")
#   expr <- readRDS("resources/expr_neuroncorrected.rds")
#   lapply(donorNames, function(d){
#     e <- expr[[d]]
#     if (!is.null(samples)){
#       labels <- samples[[d]]
#       if (is.null(genes)) e[, labels]
#       else e[genes, labels]
#     } else {
#       if (is.null(genes)) e
#       else e[genes, ]
#     }
#   })
# }

##### PD info about brain regions and genes #####

braakNames <- c("braak1","braak2","braak3","braak4","braak5","braak6")
names(braakNames) <- braakNames

pdGenes <- list(hiImpact = c("SNCA", "LRRK2", "GBA", "VPS35", "PARK2", "UCHL1", "PINK1", "PARK7", "ATP13A2", "PLA2G6", "FBXO7", "DNAJC6", "SYNJ1", 
                             "EIF4G1", "DNAJC13", "CHCHD2", "C20orf30", "RIC3", "LRP10"), #TMEM230 is C20orf30
                jansen2017 = c("INPP5F", "TMEM175", "ASH1L", "MAPT", "RIT1", "C14orf83", "STK39", "GPNMB", "BST1", 
                               "SIPA1L2", "DLG2", "NUCKS1", "GCH1", "MCCC1", "FAM47E", "BCKDK", "TMPRSS9", "UBOX5", 
                               "CCDC62", "SYNJ1", "EIF4G1", "FBXO7", "C20orf30", "POLG", "VPS13C", "PLA2G6"),
                hla = c("HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DQB1"),
                'Chang et al. 2017' = read.table("chang2017_riskgenes.txt", comment.char = "#", sep = "\n", row.names = NULL, stringsAsFactors = FALSE)[, 1], 
                'Nalls et al. 2014' = read.table("nalls2014_riskgenes.txt", comment.char = "#", sep = "\n", row.names = NULL, stringsAsFactors = FALSE)[, 1],
                parkingenes = read.table("parkin_genes_hgnc.txt", sep = "\t", header = TRUE)$Approved.Symbol
                # lysosome = read.table("lysosome_geneset.txt", header = FALSE, comment.char = "#", sep = "\n", row.names = NULL)[, 1],
                # liscovitch2014 = read.table("ifn_signaling_genes.txt", comment.char = "#", sep = "\n", row.names = NULL)[, 1]
)
pdGenesID <- lapply(pdGenes, name2EntrezId)