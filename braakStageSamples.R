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

# Function to get Braak labels
label.vector <- function(m){
  apply(m, 1, function(v){
    s <- which(v == 1)
    ifelse(length(s) == 0, 0, tail(unlist(strsplit(names(s), split = "braak")), 1))
  })
}
braakLabels <- lapply(braakStages, function(m) label.vector(m))

# Indexes to select Braak samples in order of Braak regions and anatomy
braak_idx <- lapply(donorNames, function(d){
  apply(braakStages[[d]], 2, function(b){
    samples <- which(b)
    graph_order <- sampleInfo[[d]][samples, "graph_order"]
    sample_order <- order(-graph_order)
    samples[sample_order]
  })
})
# braak_idx <- lapply(donorNames, function(d){
#   samples <- which(braakLabels[[d]] != 0)
#   labels <- braakLabels[[d]][samples]
#   graph_order <- sampleInfo[[d]][samples, "graph_order"]
#   sample_order <- order(labels, -graph_order)
#   samples[sample_order]
# })

# Fixed colors for Braak related regions
braakColors <- brewer.pal(6, "Set2")
names(braakColors) <- braakNames

save(braakStages, braakLabels, braakColors, braak_idx, file = "resources/braakInfo.RData")