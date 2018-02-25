# Braak gene enrichment of modules
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
load("resources/modules.RData")
load("resources/braakGenes.RData")

braakGenes <- unlist(braakGenes)

total = 19992
module_overlap <- lapply(modules, function(m){
  tab <- as.data.frame(t(sapply(m, function(module){
    genes <- intersect(module, braakGenes)
    overlap <- length(genes)
    ns1 <- length(module)
    ns2 <- length(braakGenes)
    p <- phyper(overlap - 1, ns1, total - ns1, ns2, lower.tail = FALSE)
    c(overlap = overlap, module.size = ns1, pvalue = p)
    # row$genes <- list(genes)
    # row
  })))
  tab$benjamini_hochberg <- p.adjust(tab$pvalue)
  tab <- tab[tab$benjamini_hochberg < 0.05, ]
  tab[order(tab$benjamini_hochberg), ]
})

lapply(module_overlap, function(b){
  # sum(b$module.size)
  b[,-4]
})

lapply(modules, function(m){
  present <- sapply(m, function(module){
    genes <- intersect(module, braakGenes)
    any(genes == "6622")
  })
  ll <- m[present]
  names(ll)
  entrezId2Name(unlist(ll))
})

lapply(names(module_overlap), function(b){
  tab <-  module_overlap[[b]]
  tab$module <- rownames(tab)
  write.table(tab, file = paste0("module_braakgene_overlap_", b, ".txt"), sep ="\t", quote = FALSE, row.names = FALSE)
})
