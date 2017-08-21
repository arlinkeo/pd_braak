# Correlated genes in DAVID

setwd("C:/Users/dkeo/surfdrive/Parkinson")
options(stringsAsFactors = FALSE)
library("RDAVIDWebService")

load("resources/geneLabelCor.RData")

# Mapping entrez IDs to gene symbols and vice versa
probeInfo <- read.csv("../ABA_human_processed/probe_info_2014-11-11.csv")
entrezId2Name <- function (x) {probeInfo$gene_symbol[match(x, probeInfo$entrez_id)]} #Input is vector
name2entrezId <- function (x) {probeInfo$entrez_id[match(x, probeInfo$gene_symbol)]} #Input is vector

# Average correlation across brains
avgCor <- apply(geneLabelCor, 1, mean)
avgCor <- sort(avgCor)
topMin <- head(avgCor, 10)
topMax <- tail(avgCor, 10)
topMinNames <- names(topMin)
topMaxNames <- names(topMax)
names(topMin) <- entrezId2Name(topMinNames)
names(topMax) <- entrezId2Name(topMaxNames)

# Functional enrichment of up- and down-regulated genes
david<-DAVIDWebService$new(email="D.L.Keo@tudelft.nl",
                           url="https://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
bg_list <- probeInfo$entrez_id
bg <- addList(david, bg_list, idType = "ENTREZ_GENE_ID", listName = "AHBA background", listType = "Background")
bg
t <- 0.05 # EASE p-value threshold
setTimeOut(david, 200000)

result <- addList(david, topMinNames, idType = "ENTREZ_GENE_ID", listName = "top10Min", listType = "Gene")
print(result)
setCurrentBackgroundPosition(david, 1)
getFunctionalAnnotationChartFile(david, "top10MinCor_goterms.txt", threshold=t, count=2L)

result <- addList(david, topMaxNames, idType = "ENTREZ_GENE_ID", listName = "top10Max", listType = "Gene")
print(result)
setCurrentBackgroundPosition(david, 1)
getFunctionalAnnotationChartFile(david, "top10MaxCor_goterms.txt", threshold=t, count=2L)