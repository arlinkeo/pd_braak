# Differentially expressed genes in DAVID

setwd("C:/Users/dkeo/surfdrive/Parkinson")
options(stringsAsFactors = FALSE)
library("RDAVIDWebService")

load("resources/diffGenesBraak.RData") # for looking up p-value (2 diff. data structures)
