# Check for presence of PD genes among BRGs and within modules

# Variant associated PD genes
pdGenes <- list(hiImpact = c("SNCA", "LRRK2", "GBA", "VPS35", "PARK2", "PINK1", "PARK7", "ATP13A2", "PLA2G6", "FBXO7", "DNAJC6", "SYNJ1",
                             "EIF4G1", "DNAJC13", "CHCHD2", "C20orf30", "RIC3", "LRP10"), #TMEM230 is C20orf30
                jansen2017 = c("INPP5F", "TMEM175", "ASH1L", "MAPT", "RIT1", "C14orf83", "STK39", "GPNMB", "BST1",
                               "SIPA1L2", "DLG2", "NUCKS1", "GCH1", "MCCC1", "FAM47E", "BCKDK", "TMPRSS9", "UBOX5",
                               "CCDC62", "SYNJ1", "EIF4G1", "FBXO7", "C20orf30", "POLG", "VPS13C", "PLA2G6"),
                hla = c("HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DQB1"),
                'Chang et al. 2017' = read.table("../chang2017_riskgenes.txt", comment.char = "#", sep = "\n", row.names = NULL, stringsAsFactors = FALSE)[, 1],
                'Nalls et al. 2014' = read.table("../nalls2014_riskgenes.txt", comment.char = "#", sep = "\n", row.names = NULL, stringsAsFactors = FALSE)[, 1],
                parkingenes = read.table("../parkin_genes_hgnc.txt", sep = "\t", header = TRUE)$Approved.Symbol
)
pdGenesID <- lapply(pdGenes, name2EntrezId)
pdGenesID <- lapply(pdGenesID, function(x) x[!is.na(x)])

# # Boxplot for each PD-implicated gene
# plot.pdf("output/boxplot_high_impact_genes.pdf", pdGenesID$hiImpact)
# plot.pdf("output/boxplot_susceptible_genes.pdf", pdGenesID$jansen2017)
# plot.pdf("output/boxplot_HLA_genes.pdf", pdGenesID$hla)
# plot.pdf("boxplot_hemegenes.pdf", name2EntrezId(c("HBD", "HBB", "HBA1", "HBA2", "OASL")))

# pd_brgs <- read.table("pdgenes_stats.txt", sep ="\t", header = TRUE)
# plot.pdf("boxplot_pd_brgs.pdf", as.character(pd_brgs$entrez_id))


##### Find PD-mplicated genes #####

# Presence PD genes in all modules
sapply(modules[unlist(braakModules)], function(m){
  paste0(entrezId2Name(intersect(unlist(pdGenesID),m)), collapse = ", ")
})

########## Presence of PD-implicated genes as BRGs and module member ##########
tab <- lapply(names(pdGenesID), function(n){
  x <- pdGenesID[[n]]
  g <- intersect(braakGenes$entrez_id, x)
  cbind(braakGenes[braakGenes$entrez_id %in% g,], study = rep(n, length(g)))
})
tab <- Reduce(rbind, tab)
tab$study <- sapply((tab[, "entrez_id"]), function(g) paste0(tab[tab$entrez_id==g, "study"], collapse = ", "))  # Merge studies if there are duplicate genes
tab <- tab[!duplicated(tab$entrez_id), ]
tab$module <- sapply(tab$entrez_id, function(g) {
  presence <- sapply(modules, function(m){g %in% m})
  ifelse(any(presence), names(which(presence)), "-")
})
tab <- tab[order(tab$r), c(1:6,8,7)]
tab[, c(3,5)] <- round(tab[, c(3,5)], digits = 2)
tab[, c(4,6)] <- format(tab[, c(4,6)], digits = 2, scientific = TRUE)
colnames(tab) <- c("Gene symbol", "Entrez ID", "Correlation with Braak (r)", "P-value (BH-corrected)", "Fold-change", "P-value (BH-corrected)",
                   "Module member", "Reference")
write.table(tab, file = "output/pdgenes_stats.txt", sep ="\t", quote = FALSE, row.names = FALSE)
