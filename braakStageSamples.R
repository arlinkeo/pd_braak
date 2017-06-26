#sample IDs of ROI for PD
setwd("C:/Users/dkeo/surfdrive/Parkinson")
load("../polyQ_coexpression/resources/BrainExpr.RData")
source("../AHBA_ontology.R")
donorNames <- names(brainExpr)
names(donorNames) <- donorNames
ontology <- read.csv("../ABA_human_processed/Ontology_edited.csv")

structures <- c("myelencephalon", "dorsal motor nucleus of the vagus", "midbrain reticular formation", "intermediate zone, left", "intermediate zone, right",
                "pons", "locus ceruleus", "midbrain raphe nuclei", "pontine raphe nucleus", "mesencephalon", "substantia nigra", 
                "basal forebrain", "amygdala", "hippocampal formation", "substantia innominata", "basal nucleus of meynert, right",
                "basal nucleus of meynert, left", "CA2 field", "occipito-temporal gyrus", "cingulate gyrus", "temporal lobe", 
                "frontal lobe", "parietal lobe")
names(structures) <- structures

#Select anatomic region-specific samples
sampleIDs <- lapply(structures, function(s){
  row <- match(s, ontology$name)
  id <- ontology$id[row]
  rows <- grep(id, ontology$structure_id_path)
  selectIds <- ontology$id[rows]
  lapply(donorNames, function(d){
    expr <- brainExpr[[d]]
    ids <- intersect(selectIds, colnames(expr))
    cols <- colnames(expr) %in% ids
    as.integer(cols)
  })
})
sapply(sampleIDs, function(s){sapply(s, sum)})

# #Braak stage I: Medulla
# med <- ontology[grep("medulla", ontology$name), ]$id
# names(med) <- sapply(med, get.name)
# medIds <- unique(unlist(lapply(med, function(id){
#   rows <- grep(id, ontology$structure_id_path)
#   selectIds <- ontology$id[rows]
# })))
# medIds <- lapply(donorNames, function(d){
#   expr <- brainExpr[[d]]
#   ids <- intersect(medIds, colnames(expr))
#   cols <- colnames(expr) %in% ids
#   as.integer(cols)
# })
# sapply(medIds, sum)

#Braak stages
braakStages <- list(
  braak1 = c("myelencephalon"),
  braak2 = c("pons"),
  braak3 = c("substantia nigra", "basal nucleus of meynert, right", "basal nucleus of meynert, left", "CA2 field"),
  braak4 = c("amygdala", "occipito-temporal gyrus"), 
  braak5 = c("cingulate gyrus", "temporal lobe"),
  braak6 = c("frontal lobe", "parietal lobe")
)
braakStages <- lapply(braakStages, function(r){
   region <- sampleIDs[r]
   ids <- lapply(donorNames, function(d){
     tab <- lapply(region, function(s){s[[d]]})
     as.numeric(Reduce("|", tab))
   })
})
sapply(braakStages, function(s){sapply(s, sum)})

#All Braak stages
# braakStages$braak1 <- medIds
# braakStages <- braakStages[c(6, 1:5)]
allBraak <- lapply(donorNames, function(d){
  tab <- lapply(braakStages, function(s){s[[d]]})
  as.numeric(Reduce("|", tab))
})

#non-Braak regions
nonBraak <- lapply(allBraak, function(d){as.numeric(!d)})

#Braak and non-Braak
sumAllB <- sapply(allBraak, sum)
sumNonB <- sapply(nonBraak, sum)
total <- sumAllB + sumNonB
tabNonAll <- cbind(sumAllB, perc.All = sumAllB/total, sumNonB, perc.Non = sumNonB/total, total)


save(braakStages, nonBraak, file = "braakStages.RData")