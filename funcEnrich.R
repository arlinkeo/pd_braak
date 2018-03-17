# Functional enrichment of genes which expression correlate with Braak regions

setwd("C:/Users/dkeo/surfdrive/pd_braak")

source("PD/base_script.R")
load("resources/modules.RData")
load( "resources/braakGenes.RData")
# R version should be 3.25 instead of 3.3 for RDavid
library("RDAVIDWebService")
# library("plyr")

############################################################################

#Functional enrichment of  genes correlated greater or smaller than 0
david<-DAVIDWebService$new(email="D.L.Keo@tudelft.nl",
                           url="https://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
bg_list <- probeInfo$entrez_id
bg <- addList(david, bg_list, idType = "ENTREZ_GENE_ID", listName = "AHBA background", listType = "Background")
bg
t <- 0.05 # EASE p-value threshold
setTimeOut(david, 200000)

# Enrichment of positively and negatively correlated progression genes
lapply(names(braakGenes), function(r){
  genes <- braakGenes[[r]]
  result <- addList(david, genes, idType = "ENTREZ_GENE_ID", listName = r, listType = "Gene")
  print(result)
  setCurrentBackgroundPosition(david, 1)
  getFunctionalAnnotationChartFile(david, paste0("Functional_analyses/", r, "_goterms.txt"), threshold=t, count=2L)
  getClusterReportFile(david, paste0("Functional_analyses/", r, "_termclusters.txt"), type = c("Term"))
})

# Enrichment of genes in modules enriched for progression genes
lapply(names(modules[3]), function(r){
  m <- modules[[r]]
  lapply(names(m), function(l){
    genes <- m[[l]]
    result <- addList(david, genes, idType = "ENTREZ_GENE_ID", listName = paste(r, ",", l), listType = "Gene")
    print(result)
    setCurrentBackgroundPosition(david, 1)
    getFunctionalAnnotationChartFile(david, paste0("Functional_analyses/", r, "_modules/", l, "_goterms.txt"), threshold=t, count=2L)
    getClusterReportFile(david, paste0("Functional_analyses/", r, "_modules/", l, "_termclusters.txt"), type = c("Term"))
  })
})

##############################################################################

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
correctedTerms <- sapply(names(modules[[3]]), function(r){
  m <- modules[[r]]
  sapply(names(m), function(l){
    file <- paste0("Functional_analyses/", r, "_modules/", l, "_goterms.txt")
    terms <- read.csv(file, header = TRUE, sep = "\t", colClasses = "character")
    rows <- which(as.numeric(terms$Benjamini) < 0.05)
    terms <- terms[rows, c("Category", "Term", "Bonferroni", "Benjamini")]
    terms
  }, simplify = FALSE)
}, simplify = FALSE)
sapply(correctedTerms, function(x)sapply(x, nrow))

# Write text corrected go-terms
lapply(names(correctedTerms), function(r){
  lapply(names(correctedTerms[[r]]), function(l){
    df <- correctedTerms[[r]][[l]]
    write.table(df[, c("Term", "Benjamini")], 
                file = paste0("Functional_analyses/", r, "_modules/", l, "_correctedgo", ".txt"), 
                sep ="\t", row.names = FALSE)
  })
})

#print top 3 significant terms in merged braak stages
top3terms <- lapply(correctedTerms, function(r){
  dfll <- lapply(r, function(l){
    l[c(1:3),]
  })
  data <- do.call(rbind, dfll)
  data[rowSums(is.na(data)) == 0,]
})
# top3merged <- 
top3merged <- top3merged[apply(!is.na(top3merged), 1, any),]
top3merged$Benjamini <- format(as.numeric(top3merged$Benjamini), digits = 3, scientific = TRUE)
