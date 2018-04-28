# Heatmap of modules enriched for different gene sets
# Braak gene enrichment of modules
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(ggplot2)
library(reshape2)
library(biomaRt)
load("resources/modules.RData")
load("resources/braakGenes.RData")
load("resources/summaryLabelCorrEG.RData")

labelCor <- do.call(rbind.data.frame, lapply(summaryLabelCorrEG, function(g) g["summary",]))
labelCor$pvalue <- p.adjust(labelCor$pvalue, method = "BH")
rownames(labelCor) <- paste0("M", rownames(labelCor))
orderEG <- order(labelCor$r)
total <- length(ahba.genes())

# # significance stars
# star <- function(v){
#   sapply(v, function(x) if (x<0.001) "***" else if (x<0.01) "**" else if (x<0.05) "*" else "")
# }

# hypergeometric test
hyper.test <- function(a, b){
  genes <- intersect(a, b)
  overlap <- length(genes)
  ns1 <- length(a)
  ns2 <- length(b)
  p <- phyper(overlap - 1, ns1, total - ns1, ns2, lower.tail = FALSE)
  p
}

##########################

# Braak label correlation plot
m <- modules[["braak1-6"]]
names(m) <- rownames(labelCor)
# module_size <- sapply(m, length)
# tab <- data.frame(module = names(m), module_size)
# tab$r <- labelCor$r
# tab$pvalue <- labelCor$pvalue
# tab <- tab[orderEG, ]
# # tab$type<- rep("Braak label correlation", nrow(tab))
# tab$star <- star(tab$pvalue)
# tab$module <- factor(tab$module, levels = unique(tab$module))
# offset <- 0.1
# tab$y <- sapply(tab$r, function(x) if (x>0) x+offset else x-offset)
# 
# colPal <- c("darkblue", "white", "darkred")
# rampcols <- colorRampPalette(colors = colPal, space="Lab")(201)
# tab$color <- rampcols[as.numeric(cut(tab$r, breaks = 201))]
# tab$color<- factor(tab$color, levels = unique(tab$color))
# 
# p1 <- ggplot(tab) + 
#   geom_col(aes(x=module, y = r, fill=r)) +
#   geom_text(aes(x=module, y=y, label=star), size = 4, vjust = 0.75) +
#   scale_y_continuous(limits = c(min(tab$r)-0.1, max(tab$r)+0.1), position = "top") +
#   scale_fill_gradientn(colours = rampcols) +
#   theme(
#     axis.text = element_text(size = 8),
#     axis.ticks.x = (element_blank()),
#     axis.line.y = element_line(),
#     panel.background = element_blank(),
#     legend.position="none"
#   ) + coord_flip()
# p1
##########################

# Genes associated with GO-terms
goterms <- as.matrix(read.table("goterms_of_interest.txt", header = FALSE, comment.char = "#", sep = "~", 
                                row.names = 2))[,1]
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=91)
go_genes <- lapply(goterms, function(go){
  as.character(unlist(getBM(c('entrezgene'), filters='go', mart=ensembl, values=go)))
})
names(go_genes) <- sapply(names(go_genes), function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x))))

# Genes from studies
studies <- list('Jansen et al. 2017' = c("INPP5F", "TMEM175", "ASH1L", "MAPT", "RIT1", "C14orf83", "STK39", "GPNMB", "BST1", 
                                "SIPA1L2", "DLG2", "NUCKS1", "GCH1", "MCCC1", "FAM47E", "BCKDK", "TMPRSS9", "UBOX5", 
                                "CCDC62", "SYNJ1", "EIF4G1", "FBXO7", "C20orf30", "POLG", "VPS13C", "PLA2G6"),
                'Chang et al. 2017' = read.table("chang2017_riskgenes.txt", comment.char = "#", sep = "\n", row.names = NULL, stringsAsFactors = FALSE)[, 1], 
                'Nalls et al. 2014' = read.table("nalls2014_riskgenes.txt", comment.char = "#", sep = "\n", row.names = NULL, stringsAsFactors = FALSE)[, 1]#,
                # liscovitch2014 = read.table("ifn_signaling_genes.txt", comment.char = "#", sep = "\n", row.names = NULL, stringsAsFactors = FALSE)[, 1]
)
studies <- lapply(studies, name2EntrezId)
studies <- lapply(studies, function(x) x[!is.na(x)])

# Disease-associated genes
disease_table <- read.delim("gene_associated_with_disease.txt")
diseases <- c("Parkinson Disease", "Young onset Parkinson disease", 
              # "Autosomal Dominant Parkinsonism", "Dementia in Parkinson's disease", "Parkinsonian Disorders",   
              # "Autosomal dominant late onset Parkinson disease", "Autosomal Recessive Parkinsonism", 
              # "Parkinson Disease 6, Autosomal Recessive Early-Onset", "Parkinsonism, Juvenile", "Parkinson Disease, Familial, Type 1",
              # "Parkinsonism-Dystonia, Infantile", "Secondary Parkinson Disease", "Wolff-Parkinson-White Syndrome",
              "Alzheimer Disease, Late Onset", "Alzheimer Disease, Early Onset", "Alzheimer's Disease", 
              "Dementia", "Lewy Body Disease", "Schizophrenia", "Lysosomal Storage Diseases", "Amyotrophic Lateral Sclerosis",
              "Multiple System Atrophy", "Multiple Sclerosis", "Amyotrophic Lateral Sclerosis", 
              "Huntington Disease", "Gastrointestinal Diseases", "Anemia", 
              "Autoimmune Diseases", "Inflammation", "Prion Diseases",
              "Muscle Rigidity", "Dystonia", "Bradykinesia", "Mental Depression")
disease_genes <- sapply(diseases, function(x){
  disease_table[disease_table$diseaseName %in% x, "geneId"]
}, simplify = FALSE)

# Cell-type genes
celltype_genes <- sapply(c("Astrocytes", "Endothelial_cells", "Microglia", "Neurons", "Oligodendrocytes"), function(type){
  file = paste0("brainscope_celltypes/", type, ".txt")
  as.character(read.csv(file, header = TRUE)$entrez_id)
}, simplify = FALSE)
names(celltype_genes) <- paste0("BS_", names(celltype_genes))

celltype_genes2 <- as.list(read.delim("celltype_markers_Ritchie2018.txt", na.strings = ""))
celltype_genes2 <- lapply(celltype_genes2, toupper)
celltype_genes2 <- lapply(celltype_genes2, name2EntrezId)
celltype_genes2 <- lapply(celltype_genes2, function(x) x[!is.na(x)])

# Combine all gene lists
genelists <- Reduce(append, list(list('Progression genes' = unlist(braakGenes, use.names = FALSE)), 
                                 go_genes, celltype_genes, studies, disease_genes))
data.frame(sapply(genelists, length))

# Filter for genes present in AHBA
genelists <- lapply(genelists, function(x){
  x[x %in% ahba.genes()]
})
data.frame(sapply(genelists, length))

names(genelists) <- paste0(names(genelists), " (", sapply(genelists, length), ")")


# Overlap with hypergeometrix test
overlap <- sapply(names(genelists), function(n){
  type_genes <- genelists[[n]]
  pvalue <- sapply(m, function(mod_genes){
    hyper.test(mod_genes, type_genes)
  })
  pvalue <- p.adjust(pvalue, method = "BH") # corrected
  pvalue
})
colnames(overlap) 
overlap <- overlap[orderEG, ]
overlap <- overlap[labelCor$pvalue[orderEG] < 0.001, ]

#Heatmap
overlap <- ifelse(overlap < 0.05, "p<0.05", "p>=0.5")
overlap <- melt(overlap)
overlap$r <- sapply(overlap$Var1, function(m){
  cor <- labelCor[as.character(m), "r"]
  if (cor > 0) "r>0" else "r<0"
})
overlap$Var2 <- factor(overlap$Var2, levels = rev(unique(overlap$Var2)))
overlap$type <- sapply(overlap$Var2, function(x){
  x <- unlist(strsplit(as.character(x), split = " \\("))[1]
  if (x == "Progression genes") 1
  else if (x %in% names(celltype_genes)) 2
  else if (x %in% names(go_genes)) 3
  else if (x %in% names(studies)) 4
  else if (x %in% names(disease_genes)) 5
  else 0
})

p2 <- ggplot(overlap) +
  geom_tile(aes(x=Var1, y=Var2, fill=value), colour = "black") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0, size = 8)) +
  facet_grid(type~r, scales = "free", space = "free")
p2

pdf("module_enrichment.pdf", 9, 6)
# p1
p2
dev.off()