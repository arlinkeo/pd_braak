# Braak info
# Samples IDs, fixed colors
setwd("C:/Users/dkeo/surfdrive/pd_braak")

source("PD/base_script.R")
load("resources/roiSamples.RData")
library(RColorBrewer)

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
# Print table with sample sizes
t(sapply(braakStages, function(m)apply(m, 2, sum)))
# save(braakStages, file = "resources/braakStages.RData")

# Function to get Braak labels
label.vector <- function(m){
  apply(m, 1, function(v){
    s <- which(v == 1)
    ifelse(length(s) == 0, 0, tail(unlist(strsplit(names(s), split = "braak")), 1))
  })
}

braakLabels <- lapply(braakStages, function(m) label.vector(m))
# save(braakLabels, file = "resources/braakLabels.RData")

# Fixed colors for Braak related regions
braak.colors <- brewer.pal(6, "Set2")
names(braak.colors) <- braakNames

save(braakStages, braakLabels, file = "resources/braakInfo.RData")
