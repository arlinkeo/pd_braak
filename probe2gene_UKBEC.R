# Map probe to genes UKBEC datasets

setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(WGCNA)
load("../UKBEC/expr.maps.rda",verbose=T)

regions <- c("CRBL", "FCTX", "HIPP","MEDU", "OCTX", "PUTM", "SNIG", "TCTX", "THAL", "WHMT")
regionExpr <- sapply(regions, function(r){
  fName <- paste("../UKBEC/expr_", r, ".txt", sep = "")
  read.csv(fName, header = TRUE, sep = " ", row.names = 1)
}, simplify = FALSE)

#Concatenate data matrices
exprConcat <- Reduce(cbind, regionExpr)

# Probes to gene based on concatenated data
probeSelection <- collapseRows(exprConcat, expr.map$tID, rownames(exprConcat), method = "maxRowVariance", connectivityBasedCollapsing = TRUE)

# Select probes for each data matrix
probeNames <- names(probeSelection$selectedRow)[probeSelection$selectedRow] # exprID
regionExpr <- lapply(regionExpr, function(x) x[probeNames, ])
# regionExpr <- lapply(regionExpr, function(x) t(scale(t(x), center = TRUE)) ) # Normallize data

save(regionExpr, file = "../UKBEC_RData/regionExpr.RData")

##############################################################################################

expr.map2 <- cbind(expr.map, selected_probe = probeSelection$selectedRow)

# Function to split character string of genes by comma
genes.split <- function(x) {unlist(strsplit(x, split = ","))}

# Add names of AHBA genes in t.map and expr.map
AHBAgenes <- probeInfo$gene_symbol # all AHBA genes
mappedGenes<- apply(t.map[1:10,], 1, function(t){ # For each gene, check overlap in name symbols
  genes <- unique(c(genes.split(t["Gene"]), genes.split(t["Symbol.NA31"])))
  overlap <- intersect(genes, AHBAgenes)
  if (length(overlap) != 1) NA else overlap
})
t.map2 <- cbind(t.map, AHBA = mappedGenes)
expr.map2 <- cbind(expr.map, AHBA = t.map2$AHBA[match(expr.map$tID, t.map2$tID)], Gene = t.map2$Gene)
save(t.map2, expr.map2, file = "../UKBEC/expr.maps2.RData")