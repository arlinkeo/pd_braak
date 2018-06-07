# Heatmap of modules enriched for different gene sets
# Braak gene enrichment of modules
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(ggplot2)
library(reshape2)
library(biomaRt)
library(GO.db)
load("resources/modules.RData")
load("resources/braakGenes.RData")
load("resources/summaryLabelCorrEG.RData")

labelCor <- do.call(rbind.data.frame, lapply(summaryLabelCorrEG, function(g) g["summary",]))
labelCor$pvalue <- p.adjust(labelCor$pvalue, method = "BH")
rownames(labelCor) <- paste0("M", rownames(labelCor))
orderEG <- order(labelCor$r)
labelCor <- labelCor[orderEG, ]
total <- length(ahba.genes())

# Sorted, significant modules
m <- modules[["braak1-6"]][orderEG]
names(m) <- paste0("M", names(m))
signif_modules <- m[labelCor$pvalue < 0.001]
save(signif_modules, file = "resources/signif_modules.RData")

# # significance stars
# star <- function(v){
#   sapply(v, function(x) if (x<0.001) "***" else if (x<0.01) "**" else if (x<0.05) "*" else "")
# }

# Braak genes
braak <- lapply(c(positive = "pos", negative = "neg"), function(x){
  braakGenes$entrez_id[braakGenes$dir %in% x]
})

# Genes associated with GO-terms
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=92)
goterms <- unlist(getBM(c('go_id'), mart = ensembl, values = 'go'))
names(goterms) <- goterms
go_genes <- lapply(goterms, function(go){
  as.character(unlist(getBM(c('entrezgene'), filters='go', mart=ensembl, values=go)))
})
save(go_genes, file = "resources/go_genes.RData")
names(go_genes) <- Term(names(go_genes))

# Disease-associated genes
disease_table <- read.delim("gene_associated_with_disease.txt")
diseases <- unique(disease_table$diseaseName)
disease_genes <- sapply(diseases, function(x){
  disease_table[disease_table$diseaseName %in% x, "geneId"]
}, simplify = FALSE)

# Cell-type genes
celltype_genes <- sapply(c("Astrocytes", "Endothelial_cells", "Microglia", "Neurons", "Oligodendrocytes"), function(type){
  file = paste0("brainscope_celltypes/", type, ".txt")
  as.character(read.csv(file, header = TRUE)$entrez_id)
}, simplify = FALSE)

# Combine all gene lists
# genelists <- Reduce(append, list(list('Progression genes' = unlist(braakGenes, use.names = FALSE)), 
                                 # go_genes, celltype_genes, studies, disease_genes))
genelists <- list('BSGs' = braak, 'celltype' = celltype_genes, 'GO' = go_genes, 'disease' = disease_genes)

# Filter for genes present in AHBA and for sets with at least 10 genes
genelists <- lapply(genelists, function(t){
  sets <- lapply(t, function(x){
    x[x %in% ahba.genes()]
  })
  sets <- sets[sapply(sets, length) >= 10]
  # Add size of gene sets
  # names(sets) <- paste0(names(sets), " (", sapply(sets, length), ")")
  sets
})
save(genelists, file = "resources/genelists.RData")

#########################################################################

# Ran on server: module_enrichment2.R

# Overlap with hypergeometric test
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

#########################################################################

load("resources/modEnrich.RData")

# Filter significant gene sets for GO terms and diseases only
modEnrich$GO <- modEnrich$GO[apply(modEnrich$GO, 1, function(x){
  any(x < 0.05)
}), ]
modEnrich$disease <- modEnrich$disease[apply(modEnrich$disease, 1, function(x){
  any(x < 0.05)
}), ]

# Heatmaps
mod_names <- list(
  'r<0' = names(m)[labelCor$pvalue < 0.001 & labelCor$r < 0], 
  'r>0' = names(m)[labelCor$pvalue < 0.001 & labelCor$r > 0]
)

prepare.table <- function(l) {
  t <- lapply(l, function(n){
    m <- modEnrich[[n]]
    m <- ifelse(m < 0.05, "<0.05", ">=0.5")
    rownames(m) <- paste0(rownames(m), " (", sapply(genelists[[n]], length)[rownames(m)], ")") # Add gene set size
    rownames(m) <- paste0(toupper(substring(rownames(m),1,1)), substring(rownames(m),2))
    t <- lapply(names(mod_names), function(dir){
      t <- melt(m[, mod_names[[dir]]])
      t$r <- dir
      t
    })
    t <- Reduce(rbind, t)
    colnames(t) <- c("geneset", "module", "pvalue", "r")
    t$module <- factor(t$module, levels = unique(t$module))
    t$geneset <- factor(t$geneset, levels = rev(unique(t$geneset)))
    t$type <- n
    t
  })
  t <- Reduce(rbind, t)
  t$type <- factor(t$type, levels = unique(t$type))
  t
}

heat.plot <- function(t) {
  ggplot(t) +
  geom_tile(aes(x=module, y=geneset, fill=pvalue), colour = "black") +
  scale_fill_discrete(name = "P-value") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_blank(),
        legend.position = "top", legend.text = element_text(size = 10), legend.title = element_text(size = 10),
        panel.background = element_blank()
  ) + 
  facet_grid(type~r, scales = "free", space = "free")
}

# Supplementary heatmap with all GO-terms and diseases
t1 <- prepare.table(c("GO", "disease"))
levels(t1$type) <- c("GO-terms", "Diseases")
pdf("module_enrichment_supplement.pdf", 14, 150)
heat.plot(t1)
dev.off()

# Heatmap with Braak genes, cell-types, and selected GO terms and diseases
selected_GO <- c("lysosome", "synapse", "cell junction", "nervous system development", "plasma membrane",
                 "extracellular region", "immune system process", "immune response", "DNA binding", "angiogenesis",
                 "lipid metabolic process", "nucleosome", "mitochondrion", 
                 "apoptotic process")
modEnrich$GO <- modEnrich$GO[selected_GO, ]
selected_diseases <- c(
                "Alzheimer Disease, Late Onset", "Alzheimer's Disease",
                "Dementia", "Lewy Body Disease", "Schizophrenia", "Epilepsy",
                "Multiple Sclerosis", "Anemia", "Autoimmune Diseases", "Inflammation", "Mental Depression"
                )
modEnrich$disease <- modEnrich$disease[selected_diseases, ]
t2 <- prepare.table(names(modEnrich))
levels(t2$type) <- c("BSGs", "Cell-types", "GO-terms", "Diseases")

pdf("module_enrichment.pdf", 10, 7)
heat.plot(t2)
dev.off()

# Genes from studies (count genes)
studies <- list('Jansen et al. 2017' = c("INPP5F", "TMEM175", "ASH1L", "MAPT", "RIT1", "C14orf83", "STK39", "GPNMB", "BST1", 
                                         "SIPA1L2", "DLG2", "NUCKS1", "GCH1", "MCCC1", "FAM47E", "BCKDK", "TMPRSS9", "UBOX5", 
                                         "CCDC62", "SYNJ1", "EIF4G1", "FBXO7", "C20orf30", "POLG", "VPS13C", "PLA2G6"),
                'Chang et al. 2017' = read.table("chang2017_riskgenes.txt", comment.char = "#", sep = "\n", row.names = NULL, stringsAsFactors = FALSE)[, 1], 
                'Nalls et al. 2014' = read.table("nalls2014_riskgenes.txt", comment.char = "#", sep = "\n", row.names = NULL, stringsAsFactors = FALSE)[, 1]#,
                # liscovitch2014 = read.table("ifn_signaling_genes.txt", comment.char = "#", sep = "\n", row.names = NULL, stringsAsFactors = FALSE)[, 1]
)
studies <- lapply(studies, name2EntrezId)
studies <- lapply(studies, function(x) x[!is.na(x)])

pdf("module_enrichment.pdf", 9, 6)
# p1
p2
dev.off()