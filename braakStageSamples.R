# Sample IDs of ROI for PD
library(RColorBrewer)
library(ggplot2)

# Fixed colors for Braak related regions
braakColors <- brewer.pal(6, "Set2")
names(braakColors) <- braakRoi

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
sapply(roiSamples, function(x){sapply(x, length)})

####################################################################

# Braak stages
braakRegions <- list(
  R1 = c("myelencephalon"),
  R2 = c("pontine tegmentum"),
  R3 = c("substantia nigra", "basal nucleus of meynert, right", "basal nucleus of meynert, left", "CA2 field"),
  R4 = c("amygdala", "occipito-temporal gyrus"), 
  R5 = c("cingulate gyrus", "temporal lobe"),
  R6 = c("frontal lobe", "parietal lobe")
)

# Braak stage samples
braak_idx <- lapply(donorNames, function (d){
  roi <- roiSamples[[d]]
  lapply(braakRegions, function(b){
    samples <- Reduce(c, roi[b])
    graph_order <- sample_annot[[d]][samples, "graph_order"]
    sample_order <- order(-graph_order)
    samples[sample_order]
  })
})

# Sample sizes within Braak-related regions
numbers <- sapply(braak_idx, function(m) sapply(m, length))
colnames(numbers) <- sapply(colnames(numbers), function(x) gsub("donor", "Donor ", x))
anatomy <- sapply(braakRegions, function(x) {
  line <- paste(x, collapse = ", ")
  paste0(toupper(substring(line, 1, 1)), substring(line, 2))
})
numbers <- cbind('Braak stage-related regions' = rownames(numbers),
                 'Anatomical structures' = anatomy,
                 numbers,
                 'Total' = rowSums(numbers))
write.csv(numbers, file = "output/braakRoi_size.csv", row.names = FALSE)

# Vector of Braak labels for all samples for each donor
braakLabels <- lapply(donorNames, function(d){
  braak <- braak_idx[[d]]
  cols <- colnames(brainExpr[[d]])
  label <- rep("0", length(cols))
  label[braak$R1] <- "1"
  label[braak$R2] <- "2"
  label[braak$R3] <- "3"
  label[braak$R4] <- "4"
  label[braak$R5] <- "5"
  label[braak$R6] <- "6"
  names(label) <- cols
  label
})

# save(braak_idx, braakLabels, braakColors, file = "output/braakInfo.RData")