# Braak gene enrichment of modules
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
load("resources/modules_absCor.RData")
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
  })))
  tab$benjamini_hochberg <- p.adjust(tab$pvalue)
  tab <- tab[tab$benjamini_hochberg < 0.05, ]
  tab[order(tab$benjamini_hochberg), ]
})

lapply(names(module_overlap), function(b){
  tab <-  module_overlap[[b]]
  tab$module <- rownames(tab)
  write.table(tab, file = paste0("module_braakgene_overlap_", b, "_absCor.txt"), sep ="\t", quote = FALSE, row.names = FALSE)
})
modules_braak <- sapply(c("braak1", "braak6", "braak1-6"), function(b){
  modNames <- rownames(module_overlap[[b]])
  modules[[b]][modNames]
}, simplify = FALSE)
save(modules_braak, file = "resources/modules_braak_absCor.RData")

sapply(modules_braak, function(m){
  present <- sapply(names(m), function(n){
    module <- m[[n]]
    present <- sapply(pdGenes$lysosome, function(g){
      any(module == name2EntrezId(g))
    })
    present <- present[!is.na(present)]
    paste0(names(present)[present], collapse = ",")
  })
  
})

# overlap of genes between modules in B1 and B6
sapply(modules_braak$braak1, function(m1){
  sapply(modules_braak$braak6, function(m6){
    round(length(intersect(m1, m6)) / min(length(m1), length(m6)) *100, digits = 0)
  })
})
