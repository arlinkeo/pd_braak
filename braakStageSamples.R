#sample IDs of ROI for PD
setwd("C:/Users/dkeo/surfdrive/pd_braak")

source("PD/base_script.R")
load("resources/roiSamples.RData")

#Braak stages
braakRegions <- list(
  braak1 = c("myelencephalon"),
  braak2 = c("pontine tegmentum"),
  braak3 = c("substantia nigra", "basal nucleus of meynert, right", "basal nucleus of meynert, left", "CA2 field"),
  braak4 = c("amygdala", "occipito-temporal gyrus"), 
  braak5 = c("cingulate gyrus", "temporal lobe"),
  braak6 = c("frontal lobe", "parietal lobe")
)

# Braak stage samples
braakStages <- lapply(donorNames, function (d){
  roi <- roiSamples[[d]]
  sapply(braakRegions, function(b){
    m <- roi[, b, drop = FALSE]
    apply(m, 1, function(v){
      ifelse(Reduce("|", v), TRUE, FALSE)
    })
  })
})
t(sapply(braakStages, function(m)apply(m, 2, sum)))
save(braakStages, file = "resources/braakStages.RData")

# Function to get Braak labels
label.vector <- function(m){
  apply(m, 1, function(v){
    s <- which(v == 1)
    ifelse(length(s) == 0, 0, tail(unlist(strsplit(names(s), split = "braak")), 1))
  })
}

braakLabels <- lapply(braakStages, function(m) label.vector(m))
save(braakLabels, file = "resources/braakLabels.RData")

# 
# #########################################################################################################
# #All Braak stages
# allBraak <- lapply(braakStages, collapseMerge)
# sapply(allBraak, sum)
# 
# # Merged Braak stages
# braakMerged1 <- split(braakNames, ceiling(seq_along(braakNames)/2))
# names(braakMerged1) <- braakNamesMerged1
# braakMerged2 <- split(braakNames, ceiling(seq_along(braakNames)/3))
# names(braakMerged2) <- braakNamesMerged2
# braakMerged1 <- merge.per.donor(braakStages, braakMerged1)
# braakMerged2 <- merge.per.donor(braakStages, braakMerged2)
# sapply(braakMerged2, function(m)apply(m, 2 , sum))
# 
# #non-Braak regions
# nonBraakA1 <- lapply(allBraak, function(d){ifelse(!d, 1L, 0L)}) # All except braak
# nonBraakB1 <- lapply(allBraak, function(d){ d[d==0]= 1L; d}) # All
# nonBraakKnown <- c("basal part of pons", "red nucleus", "ventral tegmental area", "corpus callosum", 
#                 "midbrain reticular formation", "cerebellum")
# nonBraakC1 <- sapply(donorNames, function(d) {
#   m <- roiSamples[[d]]
#   mergeRegions(m, nonBraakKnown)
# })
# 
# # Function to exclude cerebellum
# cerebellum <- lapply(roiSamples, function(d) d[,"cerebellum"])
# excludeCb <- function(b){
#   lapply(donorNames, function(d){ # All except braak and cerebellum
#     b[[d]] - cerebellum[[d]]
#   })
# }
# 
# #non-Braak regions without cerebellum
# nonBraakA2 <- excludeCb(nonBraakA1) # All except braak and cerebellum
# nonBraakB2 <- excludeCb(nonBraakB1) # All except cerebellum
# nonBraakC2 <- excludeCb(nonBraakC1) # Known unaffected regions except cerebellum
# 
# nbnames <- c("nonBraakA1", "nonBraakA2", "nonBraakB1", "nonBraakB2", "nonBraakC1", "nonBraakC2")
# nonBraak <- mget(nbnames) # named list
# save(nonBraak, file = "resources/nonBraak.RData")