setwd("C:/Users/dkeo/surfdrive/Parkinson")
library("gplots")
options(stringsAsFactors = FALSE)

make.italic <- function(x) {as.expression(lapply(x, function(x) bquote(italic(.(x)))))}

load("../ABA_Rdata/BrainExprNorm.RData")
donorNames <- names(brainExprNorm)
names(donorNames) <- donorNames
entrezId2Name <- function (x) { row <- which(probeInfo$entrez_id == x); probeInfo[row, 4]} #Input is single element

probeInfo <- read.csv("../ABA_human_processed/probe_info_2014-11-11.csv")
pdGenes1 <- c("GBA", "LRRK2", "PINK1", "PARK7", "SNCA", "VPS35", "DNAJC13", "CHCHD2")
pdGenes2 <- c("INPP5F", "TMEM175", "ASH1L", "MAPT", "RIT1", "C14orf83", "STK39", "GPNMB", "BST1", 
              "SIPA1L2", "DLG2", "NUCKS1", "GCH1", "MCCC1", "FAM47E", "BCKDK", "TMPRSS9", "UBOX5", 
              "CCDC62", "SYNJ1", "EIF4G1", "FBXO7", "C20orf30", "POLG", "VSP13C", "PLA2G6")
pdGenes3 <- c("HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DQB1")
pdGenes4 <- unlist(read.table("lysosome_geneset.txt", header = FALSE, comment.char = "#", sep = "\n", row.names = NULL))
names(pdGenes4) <- NULL
pdCol <- c(sapply(pdGenes1, function(x) "darkgreen"), 
           sapply(pdGenes2, function(x) "orange"), 
           sapply(pdGenes3, function(x) "purple"),
           sapply(pdGenes4, function(x) "pink"))
geneInfo <- probeInfo[probeInfo$gene_symbol %in% c(pdGenes1, pdGenes2, pdGenes3, pdGenes4), ]
pdEntrezIDs <- as.character(geneInfo[ , 6])
pdGenes <- geneInfo[ , 4]
pdCol <- pdCol[pdGenes]

ontology <- read.csv("../ABA_human_processed/Ontology_edited.csv")
structures <- c("brain", "locus ceruleus", "substantia nigra", "amygdala", "hippocampal formation", "temporal lobe", 
                "cingulate gyrus", "frontal lobe", "parietal lobe", "anterior olfactory nucleus", "dorsal motor nucleus of the vagus")
structureIDs <- ontology[ontology$name %in% structures, ][ , c(1:3)]
structureIDs$name <- gsub(" ", "_", structureIDs$name)
rownames(structureIDs) <- structureIDs$acronym

#Select anatomic region-specific samples
sampleIDs <- apply(structureIDs, 1, function(id){
  print(id)
  structName <- id[2]
  ontologyRows <- grep(id[1], ontology$structure_id_path)
  selectIds <- as.character(ontology$id[ontologyRows])
  lapply(donorNames, function(d){
    expr <- brainExprNorm[[d]]
    ids <- intersect(selectIds, colnames(expr))
    cols <- colnames(expr) %in% ids
    as.numeric(cols)
  })
})
sampleSize <- sapply(sampleIDs, function(s){
  sapply(s, function(d){
    length(which(unlist(d) == 1))
  })
})
sampleSize

sampleIDs[c("Br", "AO", "LC", "10")] <- NULL

#Expr in all regional samples and donors
regionalExpr <- lapply(sampleIDs, function(s){
  res <- lapply(donorNames, function(d){
    expr <- brainExprNorm[[d]]
    expr2 <- expr[pdEntrezIDs, as.logical(s[[d]])]
  })
})

# Average expression of a polyQ gene in a structure across donors and samples.
avgExpr <- sapply(regionalExpr, function(s){
  res <- sapply(donorNames, function(d){
    apply(s[[d]], 1, mean) # Avg across region-specific samples per PQ per donor
  })
  apply(res, 1, mean) # Avg across region-specific samples and donors per PQ
})
rownames(avgExpr) <- sapply(rownames(avgExpr), entrezId2Name)

# #Co-expression of PD genes
# corMat <- lapply(regionalExpr, function(s){
#   res <- lapply(donorNames, function(d){ # co-expression per donor
#     expr <- s[[d]]
#     mat <- cor(t(expr))
#     diag(mat) <- 0
#     mat
#   })
#   meanCor <- apply(simplify2array(res), 1:2, mean)
#   rownames(meanCor) <- sapply(rownames(meanCor), entrezId2Name)
#   colnames(meanCor) <- sapply(colnames(meanCor), entrezId2Name)
#   meanCor
# })

colpalette <- colorRampPalette(c("blue", "white", "red"))(n = 200)

pdf(file = "pdGenes_expr3.pdf", 8, 12)
heatmap.2(avgExpr, hclustfun = function(y){hclust(y, method = "average")}, trace = "none", 
          col = colpalette, cellnote = round(avgExpr, digits = 2), notecol = "black", notecex = 0.6, colRow = pdCol)
# plot(geneTree, xlab="", sub="", main = paste("Mean correlation of clustered modules in ", id[3], sep =""), 
#      hang = 0.04, font = 3, axes = FALSE, ylab = "");
# labeledHeatmap(avgExpr, xLabels = colnames(avgExpr), yLabels = make.italic(rownames(avgExpr)),
#                colors = blueWhiteRed(200), textMatrix = round(avgExpr, digits = 2), zlim = c(0,14), 
#                main = "Mean expression of PD genes")
dev.off()