#sample IDs of ROI for PD
setwd("C:/Users/dkeo/surfdrive/Parkinson")

load("../ABA_Rdata/BrainExprNorm.RData")
source("PD/base_script.R")
ontology <- read.csv("../ABA_human_processed/Ontology_edited.csv")

# Regions of interest (roi)
braakRoi <- c("myelencephalon", "pontine tegmentum", "substantia nigra", "amygdala", 
              "basal nucleus of meynert, right", "basal nucleus of meynert, left", "CA2 field", 
              "occipito-temporal gyrus", "cingulate gyrus", "temporal lobe", 
              "frontal lobe", "parietal lobe") 
nonBraakRoi  <-  c("cerebellum", "basal part of pons", "red nucleus", "ventral tegmental area", 
                   "corpus callosum", "midbrain reticular formation")
regions <- c(braakRoi, nonBraakRoi)
names(regions) <- regions

# Function to select region-specific samples in all six donors

selectIds <- function(r){ # for a single structure
  row <- match(r, ontology$name)
  id <- ontology$id[row]
  rows <- grep(id, ontology$structure_id_path)
  ontology$id[rows]
}

regionIDs <- lapply(regions, selectIds)

roiSamples <- lapply(donorNames, function(d){
  expr <- brainExprNorm[[d]]
  colnames <- colnames(expr)
  v <- sapply(regionIDs, function(ids){
    ids <- intersect(ids, colnames)
    cols <- colnames %in% ids
    as.integer(cols)
  })
  rownames(v) <- colnames
  v
})
sapply(roiSamples, function(df){apply(df, 2, sum)})
save(roiSamples, file = "resources/roiSamples.RData")