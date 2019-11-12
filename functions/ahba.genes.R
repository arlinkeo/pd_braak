# Function to get  AHBA genes

ahba.genes <- function(random = NULL){
  genes <- probeInfo$entrez_id
  genes <- if (!is.null(random)) sample(genes, random) else genes
  as.character(genes)
}