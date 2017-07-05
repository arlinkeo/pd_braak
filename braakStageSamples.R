#sample IDs of ROI for PD
setwd("C:/Users/dkeo/surfdrive/Parkinson")
load("../polyQ_coexpression/resources/BrainExpr.RData")
source("../AHBA_ontology.R")
donorNames <- names(brainExpr)
names(donorNames) <- donorNames
sampleIds <- lapply(brainExpr, colnames)
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

#Braak regions 1-3 and 4-6 merged
braak1to3 <- lapply(donorNames, function(d){
  tab <- lapply(braakStages[1:3], function(s){s[[d]]})
  as.numeric(Reduce("|", tab))
})
sumB1to3 <- sapply(braak1to3, sum)
braak4to6 <- lapply(donorNames, function(d){
  tab <- lapply(braakStages[4:6], function(s){s[[d]]})
  as.numeric(Reduce("|", tab))
})
sumB4to6 <- sapply(braak4to6, sum)
cbind(sumB1to3, sumB4to6)

braakStages <- append(braakStages, list('braak1-3' = braak1to3, 'braak4-6' = braak4to6))

save(braakStages, nonBraak, file = "resources/braakStages.RData")

###### Functions to convert binary vectors to single vector with Braak labels

#Revert list in list
revert.list <- function(ll){
  lapply(donorNames, function(d){
    t(sapply(ll, function(bs){
      bs[[d]]
    }))
  })
}
#generate label vector for all samples
label.vector <- function(ll){
  lapply(donorNames, function(dn){
    d <- ll[[dn]]
    labels <- apply(d, 2, function(v){
      s <- which(v == 1)
      ifelse(length(s) == 0, 0, tail(unlist(strsplit(names(s), split = "braak")), 1))
    })
    names(labels) <- sampleIds[[dn]]
    labels
  })
}

### Braak labels 1 to 6 for each sample ###
bsPerDonor <- revert.list(braakStages[1:6])
braakLabels <- label.vector(bsPerDonor)

### Merged Braak labels 1-3 and 4-6 for each sample ###
mergedBsPerDonor <- revert.list(braakStages[7:8])
mergedBraakLabels <- label.vector(mergedBsPerDonor)

save(braakLabels, mergedBraakLabels, file = "resources/braakLabels.RData")

#Count non-zero labels, if >1 there are multiple labels assigned to a sample
countLabels <- sapply(bsPerDonor, function(d){apply(d, 2, function(v){sum(v!=0)})})
lapply(countLabels, function(d){which(d > 2)})