#sample IDs of ROI for PD
setwd("C:/Users/dkeo/surfdrive/pd_braak")
library(RColorBrewer)
source("PD/base_script.R")
load("../ABA_Rdata/BrainExpr.RData")
source("PD/sample.ids.R")

# Regions of interest (roi)
roi <- c("myelencephalon", "pontine tegmentum", "substantia nigra", "amygdala", 
         "basal nucleus of meynert, right", "basal nucleus of meynert, left", "CA2 field", 
         "occipito-temporal gyrus", "cingulate gyrus", "temporal lobe", 
         "frontal lobe", "parietal lobe",
         #non-braak
         "cerebellum", "basal part of pons", "red nucleus", "ventral tegmental area",
         "corpus callosum", "midbrain reticular formation"
)

# Sample IDs of roi's
roiIDs <- sample.ids(roi)

# Function to get indices to select samples/columns per donor
samples.donor <- function(ids, d){ # list of IDs and donorname
  expr <- brainExpr[[d]]
  colnames <- colnames(expr)
  ids <- intersect(ids, colnames)
  cols <- which(colnames %in% ids)
  names(cols) <- colnames[cols]
  cols # column indices
}

roiSamples <- lapply(donorNames, function(d){
  lapply(roiIDs, function(ids){
    samples.donor(ids, d)
  })
})
sapply(roiSamples, function(x){lapply(x, length)})
save(roiSamples, file = "roiSamples.RData")

####################################################################

# Braak stages
braakRegions <- list(
  braak1 = c("myelencephalon"),
  braak2 = c("pontine tegmentum"),
  braak3 = c("substantia nigra", "basal nucleus of meynert, right", "basal nucleus of meynert, left", "CA2 field"),
  braak4 = c("amygdala", "occipito-temporal gyrus"), 
  braak5 = c("cingulate gyrus", "temporal lobe"),
  braak6 = c("frontal lobe", "parietal lobe")
)

# Braak stage samples
braak_idx <- lapply(donorNames, function (d){
  roi <- roiSamples[[d]]
  lapply(braakRegions, function(b){
    samples <- Reduce(c, roi[b])
    graph_order <- sampleInfo[[d]][samples, "graph_order"]
    sample_order <- order(-graph_order)
    samples[sample_order]
  })
})

# Print table with sample sizes
t(sapply(braak_idx, function(m) sapply(m, length)))

# Function to get Braak labels
# label.vector <- function(m){
#   apply(m, 1, function(v){
#     s <- which(v == 1)
#     ifelse(length(s) == 0, 0, tail(unlist(strsplit(names(s), split = "braak")), 1))
#   })
# }
# braakLabels <- lapply(braakStages, function(m) label.vector(m))

braakLabels <- lapply(donorNames, function(d){
  braak <- braak_idx[[d]]
  cols <- colnames(brainExpr[[d]])
  label <- rep("0", length(cols))
  label[braak$braak1] <- "1"
  label[braak$braak2] <- "2"
  label[braak$braak3] <- "3"
  label[braak$braak4] <- "4"
  label[braak$braak5] <- "5"
  label[braak$braak6] <- "6"
  names(label) <- cols
  label
})

# Fixed colors for Braak related regions
braakColors <- brewer.pal(6, "Set2")
names(braakColors) <- braakNames

save(braak_idx, braakLabels, braakColors, file = "resources/braakInfo.RData")
