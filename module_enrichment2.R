setwd("/tudelft.net/staff-bulk/ewi/insy/DBL/Arlin/pd_braak/PD")

load("../resources/signif_modules.RData")
load("../resources/genelists.RData")

# hypergeometric test
total <- 19992
hyper.test <- function(a, b){
  genes <- intersect(a, b)
  overlap <- length(genes)
  ns1 <- length(a)
  ns2 <- length(b)
  p <- phyper(overlap - 1, ns1, total - ns1, ns2, lower.tail = FALSE)
  p
}

modEnrich <- sapply(names(genelists), function(n1){
  t <- genelists[[n1]]
  t(sapply(names(t), function(n2){
    print(paste0(n1, "; ", which(names(t)==n2), ": ", n2))
    set <- t[[n2]]
    # Overlap with each module
    pvalue <- sapply(signif_modules, function(mod_genes){
      hyper.test(mod_genes, set)
    }, simplify = FALSE)
    pvalue <- p.adjust(pvalue, method = "BH") # corrected
    pvalue
  }))
}, simplify = FALSE)
save(modEnrich, file = "../resources/modEnrich.RData")