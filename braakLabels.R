# Braak labels
setwd("C:/Users/dkeo/surfdrive/pd_braak")

source("PD/base_script.R")
load("resources/braakStages.RData")

label.vector <- function(m){
  apply(m, 1, function(v){
    s <- which(v == 1)
    ifelse(length(s) == 0, 0, tail(unlist(strsplit(names(s), split = "braak")), 1))
  })
}

braakLabels <- lapply(braakStages, function(m) label.vector(m))
braakLabelsMerged1 <- lapply(braakMerged1, function(m) label.vector(m))
braakLabelsMerged2 <- lapply(braakMerged2, function(m) label.vector(m))

save(braakLabels, braakLabelsMerged1, braakLabelsMerged2, file = "resources/braakLabels.RData")
