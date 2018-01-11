# Functional enrichment of genes which expression correlate with Braak regions

setwd("C:/Users/dkeo/surfdrive/pd_braak")

source("PD/base_script.R")

# R version should be 3.25 instead of 3.3 for RDavid
library("ggplot2")
library("RDAVIDWebService")
load("resources/braakGenes.RData")

#Functional enrichment of  genes correlated greater or smaller than 0
david<-DAVIDWebService$new(email="D.L.Keo@tudelft.nl",
                           url="https://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
bg_list <- probeInfo$entrez_id
bg <- addList(david, bg_list, idType = "ENTREZ_GENE_ID", listName = "AHBA background", listType = "Background")
bg
t <- 0.05 # EASE p-value threshold
setTimeOut(david, 200000)

# Enrichment of genes correlated across braak stages
lapply(names(braakGenes), function(l){
  genes <- braakGenes[[l]]
  result <- addList(david, genes, idType = "ENTREZ_GENE_ID", listName = l, listType = "Gene")
  print(result)
  setCurrentBackgroundPosition(david, 1)
  getFunctionalAnnotationChartFile(david, paste0("Functional_analyses/", l, "_goterms.txt"), threshold=t, count=2L)
  getClusterReportFile(david, paste0("Functional_analyses/", l, "_termclusters.txt"), type = c("Term"))
})

#Read DAVID output


#Function to read Rdavid output
read.RdavidOutput <- function(fileName){
  if (file.exists(fileName)){
    terms <- read.csv(fileName, header = TRUE, sep = "\t", colClasses = "character")
    if (nrow(terms) == 0){
      print("...Removed")
      file.remove(fileName)
      NULL
    } else {
      terms
    }
  } else {
    NULL
  }
}

#Benjamini-corrected GO terms
correctedTerms <- lapply(names(braakGenes), function(l){
  file <- paste0("Functional_analyses/", l, "_goterms.txt")
  terms <- read.csv(file, header = TRUE, sep = "\t", colClasses = "character")
  rows <- which(as.numeric(terms$Benjamini) < 0.05)
  terms <- terms[rows, c("Category", "Term", "Bonferroni", "Benjamini")]
  terms
})
sapply(correctedTerms, nrow)

#print top 3 significant terms in merged braak stages
top3terms <- lapply(correctedTerms, function(p){
  p[c(1:3),]
})
top3merged <- do.call(rbind, top3terms)
top3merged <- top3merged[apply(!is.na(top3merged), 1, any),]
top3merged$Benjamini <- format(as.numeric(top3merged$Benjamini), digits = 3, scientific = TRUE)
