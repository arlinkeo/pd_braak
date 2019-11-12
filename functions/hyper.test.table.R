# Test for two lists of gene sets and correct P for cell-types tested
hyper.test.table <- function(l1, l2){ # two lists of gene sets
  pvalue <- sapply(names(l1), function(n){
    print(paste0(which(names(l1)==n), ": ", n))
    set <- l1[[n]]
    # Overlap with each module
    sapply(l2, function(mod_genes){
      hyper.test(mod_genes, set, length(ahba.genes()))
    })
  })
  apply(pvalue, 2, function(x) p.adjust(x, method = "BH"))
}
