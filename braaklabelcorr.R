#Correlate expression  with braakstage labels

options(stringsAsFactors = FALSE)
setwd("C:/Users/dkeo/surfdrive/Parkinson")
load("resources/braakStages.RData")
braakNames <- names(braakStages)
names(braakNames) <- braakNames
donorNames <- names(braakStages$braak1)
names(donorNames) <- donorNames

bsPerDonor <- lapply(donorNames, function(d){
  t(sapply(braakStages, function(bs){
    bs[[d]]
  }))
})

#Count non-zero labels, if >1 there are multiple labels assigned to a sample
countLabels <- sapply(bsPerDonor, function(d){apply(d, 2, function(v){sum(v!=0)})})
lapply(countLabels, function(d){which(d > 1)})

braakLabels <- sapply(bsPerDonor, function(d){
  (apply(d, 2, function(v){
    s <- which(v == 1)
    ifelse(length(s) == 0, 0, s)
  }))
})
t(sapply(braakLabels, function(bl) {
  sapply(c(1:6), function(s){
    sum(bl == s)
  })
}))