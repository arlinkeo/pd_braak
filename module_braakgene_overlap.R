# Braak gene enrichment of modules
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
load("resources/modules.RData")
load("resources/braakGenes.RData")
load("resources/summaryLabelCorrEG.RData")

braakGenes <- unlist(braakGenes)
labelCor <- do.call(rbind.data.frame, lapply(summaryLabelCorrEG, function(g) g["summary",]))
total <- length(ahba.genes())
orderEG <- rev(order(labelCor$r))

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
# tab <- tab[tab$benjamini_hochberg < 0.05, ] # select only significant modules
# tab[order(tab$benjamini_hochberg), ]
tab$r_braak <- labelCor$r

# modules_braak <- sapply(c("braak1", "braak6", "braak1-6"), function(b){
#   modNames <- rownames(module_overlap[[b]])
#   modules[[b]][modNames]
# }, simplify = FALSE)
# save(modules_braak, file = "resources/modules_braak.RData")

# Add enrichment of PD genes
pdGenesID$progression <- braakGenes
pdOverlap <-  sapply(pdGenesID, function(pd){
  sapply(m, function(genes){
    res <- intersect(genes, pd)
    paste0(entrezId2Name(res), collapse = ",")
  })
})

tab <- cbind(tab, pdOverlap)
tab <- tab[orderEG, ]


# Print results in text-file
tab$pvalue <- NULL
tab <- cbind(module = rownames(tab), tab)
write.table(tab, file = paste0("module_braakgene_overlap_", b, ".txt"), sep ="\t", quote = FALSE, row.names = FALSE)
