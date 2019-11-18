# Map probes to genes
# First download 'Complete normalized microarray datasets' for all six donors from http://human.brain-map.org/static/download
# Store the downloaded files in a AHBA download directory
library(WGCNA)

########## Load data ##########

# Set location of AHBA directories
ahba_download <- "C:/Users/dkeo/surfdrive/AHBA_download" # To store downloaded files
ahba_dir <- "C:/Users/dkeo/surfdrive/AHBA_Arlin" # To store processed files as done in this script

# Read ontology and probe info (files are the same for each donor)
ontology <- read.csv(paste0(ahba_download, "/normalized_microarray_donor9861/Ontology.csv")) # Same for each donor
probeInfo <- read.csv(paste0(ahba_download, "/normalized_microarray_donor9861/Probes.csv")) # Same for each donor

# Read expression data for each donor (probes x samples)
brainExpr <- lapply(donorNames, function(d){
  file1 <- paste0(ahba_download, "/normalized_microarray_", d, "/MicroarrayExpression.csv")
  e <- read.csv(file1, header = FALSE)
  rownames(e) <- e[,1]
  e <- e[,-1]
  file2 <- paste0(ahba_download, "/normalized_microarray_", d, "/SampleAnnot.csv")
  sample_annotation <- read.csv(file2)
  colnames(e) <- sample_annotation$structure_id
  e
})

########## Filter probes based on concatenated info from all donors ##########

# Probes with missing Entrez IDs
probes_missing_entrezID <- is.na(probe_info$entrez_id) # same in all donors
sum(probes_missing_entrezID)

# Probes with expression well above background in at least 1% of samples in all donors
pa_call <- lapply(donorNames, function(d){
  file <- paste0(ahba_download, "/normalized_microarray_", d, "/PACall.csv")  
  pa <- read.csv(file, header = FALSE)
  rownames(pa) <- pa[,1]
  pa[,-1]
})
pa_concat <- Reduce(cbind, pa_call)
sum <- rowSums(pa_concat)
low_presence <- sum < 0.01*ncol(pa_concat)
print(paste(sum(low_presence), "probes present in <1% of samples"))

# Concatenate data across all donors and filter probes by  ID
expr_concat <- Reduce(cbind, brainExpr)
expr_concat <- expr_concat[!(probes_missing_entrezID | low_presence),]
probeInfo <- probeInfo[!(probes_missing_entrezID | low_presence), ]

probe2gene <- collapseRows(expr_concat, 
                           rowGroup = probeInfo$entrez_id, 
                           rowID = probeInfo$probe_id,
                           method = "maxRowVariance",
                           connectivityBasedCollapsing = TRUE
                           )
selected_probes <- probe2gene$selectedRow
entrez_id <- probeInfo$entrez_id[selected_probes]
entrez_order <- order(entrez_id)

# Filter probes in original data and select probes based on probe2gene mapping
brainExpr <- lapply(brainExpr, function(e){
  e <- e[!(probes_missing_entrezID | low_presence),]
  e <- e[selected_probes, ]
  rownames(e) <- entrez_id
  e[entrez_order,]
})
saveRDS(brainExpr, paste0(ahba_dir, "/gene_expr.RDS"))

# Write csv-files
lapply(donorNames, function(d){
  e <- brainExpr[[d]]
  write.csv(e, file = paste0(ahba_dir, "/gene_expr_", d, "_2018-11-18.csv"))
})

probeInfo <- probeInfo[selected_probes, ]
probeInfo <- probeInfo[entrez_order,]
write.csv(probeInfo, file = paste0(ahba_dir, "/probe_info_2018-11-18.csv"), row.names = FALSE)

# Copy sample annotation info to ahba_dir
# Files do not change
lapply(donorNames, function(d){
  file <- paste0(ahba_download, "/normalized_microarray_", d, "/SampleAnnot.csv")
  sample_annotation <- read.csv(file)
  write.csv(sample_annotation, file = paste0(ahba_dir, "/sample_info_", d, "_2018-11-18.csv"), row.names = FALSE)
})

# Copy ontology info to ahba_dir
write.csv(ontology, file = paste0(ahba_dir, "/Ontology.csv"), row.names = FALSE)

# # Load AHBA data for analyses
# probeInfo <- read.csv(paste0(ahba_dir, "/probe_info_2018-11-18.csv"))
# brainExpr <- readRDS(paste0(ahba_dir, "/gene_expr.RDS"))
# ontology <- read.csv(paste0(ahba_dir, "/Ontology.csv"))
sample_annot <- lapply(donorNames, function(d){ # Sample info per donor
  t <- read.csv(paste0(ahba_dir, "/sample_info_", d, "_2018-11-18.csv"))
  info <- ontology[match(t$structure_id, ontology$id), ]
  info$color_hex_triplet <- paste0("#", info$color_hex_triplet)
  info
})
