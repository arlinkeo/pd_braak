#sample IDs of ROI for PD
setwd("C:/Users/dkeo/surfdrive/pd_braak")

load("../ABA_Rdata/BrainExprNorm.RData")
source("PD/base_script.R")
ontology <- read.csv("../ABA_human_processed/Ontology_edited.csv")

# Regions of interest (roi)
roi <- c("myelencephalon", "pontine tegmentum", "substantia nigra", "amygdala", 
              "basal nucleus of meynert, right", "basal nucleus of meynert, left", "CA2 field", 
              "occipito-temporal gyrus", "cingulate gyrus", "temporal lobe", 
              "frontal lobe", "parietal lobe",
         #non-braak
         "cerebellum", "basal part of pons", "red nucleus", "ventral tegmental area",
                   "corpus callosum", "midbrain reticular formation"
         )

# Function to select region-specific sample IDs in all six donors
selectIds <- function(r){ # for a single structure
  row <- match(r, ontology$name)
  id <- ontology$id[row]
  rows <- grep(id, ontology$structure_id_path)
  ontology$id[rows]
}

# Sample IDs of roi's
roiIDs <- sapply(roi, selectIds, simplify = FALSE)

# Function to get T/F vectors to select samples/columns per donor
samples.donor <- function(ids, d){ # list of IDs and donorname
    expr <- brainExprNorm[[d]]
    colnames <- colnames(expr)
    ids <- intersect(ids, colnames)
    cols <- colnames %in% ids
    names(cols) <- colnames
    cols # as.integer(cols)
}

roiSamples <- lapply(donorNames, function(d){
  sapply(roiIDs, function(ids){
    samples.donor(ids, d)
  })
})
sapply(roiSamples, function(df){apply(df, 2, sum)})
save(roiSamples, file = "resources/roiSamples.RData")
