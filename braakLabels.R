# Braak labels
setwd("C:/Users/dkeo/surfdrive/Parkinson")

source("PD/base_script.R")
load("resources/braakStages.RData")

########################################################################
###### Functions to convert binary vectors to single vectors with Braak labels

# #Revert list in list
# revert.list <- function(ll){
#   lapply(donorNames, function(d){
#     t(sapply(ll, function(bs){
#       bs[[d]]
#     }))
#   })
# }
# #generate label vector for all samples
# label.vector <- function(ll){
#   lapply(donorNames, function(dn){
#     d <- ll[[dn]]
#     labels <- apply(d, 2, function(v){
#       s <- which(v == 1)
#       ifelse(length(s) == 0, 0, tail(unlist(strsplit(names(s), split = "braak")), 1))
#     })
#     names(labels) <- sampleIds[[dn]]
#     labels
#   })
# }
# 
# ### Braak labels 1 to 6 for each sample ###
# bsPerDonor <- revert.list(braakStages[c("braak1", "braak2", "braak3", "braak4", "braak5", "braak6")])
# braakLabels <- label.vector(bsPerDonor)
# ### Merged Braak labels 1-2, 3-4, and 4-6 for each sample ###
# bsPerDonor1 <- revert.list(braakStages[c("braak1-2", "braak3-4", "braak5-6")])
# braakLabels1 <- label.vector(bsPerDonor1)
# ### Merged Braak labels 1-3 and 4-6 for each sample ###
# bsPerDonor2 <- revert.list(braakStages[c("braak1-3", "braak4-6")])
# braakLabels2 <- label.vector(bsPerDonor2)
# 
# save(braakLabels, braakLabels1, braakLabels2, file = "resources/braakLabels.RData")
# 
# #Count non-zero labels, if >1 there are multiple labels assigned to a sample
# countLabels <- sapply(bsPerDonor, function(d){apply(d, 2, function(v){sum(v!=0)})})
# lapply(countLabels, function(d){which(d > 2)})

########################################################################

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
