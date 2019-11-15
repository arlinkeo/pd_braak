# Map probes to genes
setwd("C:/Users/dkeo/surfdrive/AHBA_Arlin")
library(WGCNA)

########## Load data ##########

# Brain donor names
donorNames <- c("donor9861", "donor10021", "donor12876", "donor14380", "donor15496", "donor15697")
names(donorNames) <- donorNames

# Read ontology and probe info (files are the same for each donor)
ontology <- read.csv("../AHBA_download/normalized_microarray_donor9861/Ontology.csv")
probe_info <- read.csv("../AHBA_download/normalized_microarray_donor9861/Probes.csv")

# Read expression data
expr <- lapply(donorNames, function(d){
  file1 <- paste0("../AHBA_download/normalized_microarray_", d, "/MicroarrayExpression.csv")
  e <- read.csv(file1, header = FALSE)
  rownames(e) <- e[,1]
  e <- e[,-1]
  file2 <- paste0("../AHBA_download/normalized_microarray_", d, "/SampleAnnot.csv")
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
  file <- paste0("../AHBA_download/normalized_microarray_", d, "/PACall.csv")  
  pa <- read.csv(file, header = FALSE)
  rownames(pa) <- pa[,1]
  pa[,-1]
})
pa_concat <- Reduce(cbind, pa_call)
sum <- rowSums(pa_concat)
low_presence <- sum < 0.01*ncol(pa_concat)
print(paste(sum(low_presence), "probes present in <1% of samples"))

# Concatenate data across all donors and filter probes by  ID
expr_concat <- Reduce(cbind, expr)
expr_concat <- expr_concat[!(probes_missing_entrezID | low_presence),]
probe_info <- probe_info[!(probes_missing_entrezID | low_presence), ]

probe2gene <- collapseRows(expr_concat, 
                           rowGroup = probe_info$entrez_id, 
                           rowID = probe_info$probe_id,
                           method = "maxRowVariance",
                           connectivityBasedCollapsing = TRUE
                           )
selected_probes <- probe2gene$selectedRow
entrez_id <- probe_info$entrez_id[selected_probes]
entrez_order <- order(entrez_id)

# Filter probes in original data and select probes based on probe2gene mapping
expr <- lapply(expr, function(e){
  e <- e[!(probes_missing_entrezID | low_presence),]
  e <- e[selected_probes, ]
  rownames(e) <- entrez_id
  e[entrez_order,]
})
saveRDS(expr, "gene_expr.RDS")

# Write csv-files
lapply(donorNames, function(d){
  e <- expr[[d]]
  write.csv(e, file = paste0("gene_expr_", d, "_2018-11-18.csv"))
})

probe_info <- probe_info[selected_probes, ]
probe_info <- probe_info[entrez_order,]
write.csv(probe_info, file = "probe_info_2018-11-18.csv", row.names = FALSE)

# Copy sample annotation info to dir
# Files do not change
lapply(donorNames, function(d){
  file <- paste0("../AHBA_download/normalized_microarray_", d, "/SampleAnnot.csv")
  sample_annotation <- read.csv(file)
  write.csv(sample_annotation, file = paste0("sample_info_", d, "_2018-11-18.csv"), row.names = FALSE)
})
