# Braak gene enrichment of modules
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
load("resources/modules.RData")
load("resources/braakGenes.RData")
load("resources/summaryLabelCorrEG.RData")

braakGenes <- unlist(braakGenes, use.names = FALSE)
labelCor <- do.call(rbind.data.frame, lapply(summaryLabelCorrEG, function(g) g["summary",]))
total <- length(ahba.genes())

b= "braak1-6"
m <- modules[[b]]
module_size <- sapply(m, length)

tab <- data.frame(module = paste0("M", names(m)), module_size)

# Add eigengene-braak correlations
tab$eigengene_r <- labelCor$r

# Add numbers of significant GO terms
correctedTerms <- sapply(names(m), function(l){
  file <- paste0("Functional_analyses/", b, "_modules/", l, "_goterms.txt")
  terms <- read.csv(file, header = TRUE, sep = "\t", colClasses = "character")
  rows <- which(as.numeric(terms$Benjamini) < 0.05)
  terms <- terms[rows, c("Category", "Term", "Bonferroni", "Benjamini")]
  terms
}, simplify = FALSE)
tab$number_goterms <- sapply(correctedTerms, nrow)

# Overlap with Braak genes and cell types
hyper.test <- function(a, b){
  genes <- intersect(a, b)
  overlap <- length(genes)
  ns1 <- length(a)
  ns2 <- length(b)
  p <- phyper(overlap - 1, ns1, total - ns1, ns2, lower.tail = FALSE)
  c(overlap = overlap, pvalue = p)
}

# Cell-type enrichment numbers
geneLists <- append(list(progression = braakGenes), celltype_genes)
overlap <- lapply(names(geneLists), function(n){
  type_genes <- geneLists[[n]]
  df <- as.data.frame(t(sapply(m, function(mod_genes){
    hyper.test(mod_genes, type_genes)
  })))
  df$pvalue <- p.adjust(df$pvalue, method = "BH") # corrected
  colnames(df) <- sapply(c("Overlap_", "pvalue_"), function(x) paste0(x, n))
  df
})
overlap <- Reduce(cbind, overlap)
overlap <- overlap[, order(!c(1:ncol(overlap))%%2)]
tab <- cbind(tab, overlap)

# Add lists of PD-implicated genes (and cell-types)
geneLists <- Reduce(append, list(pdGenesID, list(progression = braakGenes), celltype_genes))
pdOverlap <-  sapply(geneLists, function(pd){
  sapply(m, function(genes){
    res <- intersect(genes, pd)
    paste0(entrezId2Name(res), collapse = ",")
  })
})
tab <- cbind(tab, pdOverlap)

# Order genes based on summary correlation with Braak labels
orderEG <- rev(order(labelCor$r))
tab <- tab[orderEG, ]

# Print results in text-file
write.table(tab, file = paste0("module_braakgene_overlap_", b, ".txt"), sep ="\t", quote = FALSE, row.names = FALSE)