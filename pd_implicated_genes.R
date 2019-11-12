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