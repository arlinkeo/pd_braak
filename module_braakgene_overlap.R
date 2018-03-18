# Braak gene enrichment of modules
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
load("resources/modules.RData")
load("resources/braakGenes.RData")
load("resources/summaryLabelCorrEG.RData")

braakGenes <- unlist(braakGenes)
labelCor <- do.call(rbind.data.frame, lapply(summaryLabelCorrEG, function(g) g["summary",]))
total <- length(ahba.genes())

# Overlap with Braak genes
b= "braak1-6"
m <- modules[[b]]
tab <- as.data.frame(t(sapply(m, function(module){ # for each module
  genes <- intersect(module, braakGenes)
  overlap <- length(genes)
  ns1 <- length(module)
  ns2 <- length(braakGenes)
  p <- phyper(overlap - 1, ns1, total - ns1, ns2, lower.tail = FALSE)
  c(overlap = overlap, module.size = ns1, pvalue = p)
})))
tab$benjamini_hochberg <- p.adjust(tab$pvalue)

# Add numbers of significant GO terms
m <- modules[[b]]
correctedTerms <- sapply(names(m), function(l){
  file <- paste0("Functional_analyses/", b, "_modules/", l, "_goterms.txt")
  terms <- read.csv(file, header = TRUE, sep = "\t", colClasses = "character")
  rows <- which(as.numeric(terms$Benjamini) < 0.05)
  terms <- terms[rows, c("Category", "Term", "Bonferroni", "Benjamini")]
  terms
}, simplify = FALSE)
tab$numer_goterms <- sapply(correctedTerms, nrow)

# Add eigengene-braak correlations
tab$r_braak <- labelCor$r

# Add enrichment of PD-implicated genes
pdGenesID$progression <- braakGenes
pdOverlap <-  sapply(pdGenesID, function(pd){
  sapply(m, function(genes){
    res <- intersect(genes, pd)
    paste0(entrezId2Name(res), collapse = ",")
  })
})
tab <- cbind(tab, pdOverlap)

orderEG <- rev(order(labelCor$r))
tab <- tab[orderEG, ]

# Print results in text-file
tab$pvalue <- NULL
tab <- cbind(module = rownames(tab), tab)
write.table(tab, file = paste0("module_braakgene_overlap_", b, ".txt"), sep ="\t", quote = FALSE, row.names = FALSE)