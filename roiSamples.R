#sample IDs of ROI for PD
setwd("C:/Users/dkeo/surfdrive/Parkinson")

load("../polyQ_coexpression/resources/BrainExpr.RData")
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
selectSamples <- function(r){# for a single structure
  row <- match(r, ontology$name)
  id <- ontology$id[row]
  rows <- grep(id, ontology$structure_id_path)
  selectIds <- ontology$id[rows]
  lapply(donorNames, function(d){
    expr <- brainExpr[[d]]
    ids <- intersect(selectIds, colnames(expr))
    cols <- colnames(expr) %in% ids
    as.integer(cols)
  })
}

roiSamples <- lapply(regions, function(s){
  selectSamples(s)
})
sapply(roiSamples, function(s){sapply(s, sum)})
save(roiSamples, file = "resources/roiSamples.RData")