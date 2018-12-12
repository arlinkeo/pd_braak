# Heatmap of modules enriched for different gene sets
# Braak gene enrichment of modules
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(ggplot2)
library(reshape2)
library(biomaRt)
# library(GO.db)
load("resources/modules.RData")
load("resources/braakGenes.RData")
load("resources/braakModules.RData")
# load("resources/braakModules2.RData")

# dget("PD/diff.expr.lm.R")

########## Prepare list of Braak genes, cell types, GO-term, and diseases ##########

# Braak genes
braak <- list( 
  downregulated = braakGenes$entrez_id[braakGenes$r < 0],
  upregulated = braakGenes$entrez_id[braakGenes$r > 0])

# Genes associated with GO-terms
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=92)
goterms <- unlist(getBM(c('go_id'), mart = ensembl, values = 'go'))
names(goterms) <- goterms
go_genes <- lapply(goterms, function(go){
  as.character(unlist(getBM(c('entrezgene'), filters='go', mart=ensembl, values=go)))
})
save(go_genes, file = "resources/go_genes.RData")
names(go_genes) <- Term(names(go_genes))
go_genes <- go_genes[which(!is.na(names(go_genes)))]

# Disease-associated genes
disease_table <- read.delim("gene_associated_with_disease.txt")
diseases <- unique(disease_table$diseaseName)
disease_genes <- sapply(diseases, function(x){
  disease_table$geneId[disease_table$diseaseName %in% x]
}, simplify = FALSE)

# Cell-type genes
celltype_genes <- sapply(c("Astrocytes", "Endothelial_cells", "Microglia", "Neurons", "Oligodendrocytes"), function(type){
  file = paste0("brainscope_celltypes/", type, ".txt")
  as.character(read.csv(file, header = TRUE)$entrez_id)
}, simplify = FALSE)

# Combine all gene lists
genelists <- list('BSGs' = braak, 'celltype' = celltype_genes, 'GO' = go_genes, 'disease' = disease_genes)

# Filter for genes present in AHBA and for sets with at least 10 genes
genelists <- lapply(genelists, function(t){
  sets <- lapply(t, function(x){
    x[x %in% ahba.genes()]
  })
  sets <- sets[sapply(sets, length) >= 10]
  sets
})
save(genelists, file = "resources/genelists.RData")
load("resources/genelists.RData")
    
########## Module enrichment functions ##########

# hypergeometric test
hyper.test <- function(a, b, total){
  genes <- intersect(a, b)
  overlap <- length(genes)
  ns1 <- length(a)
  ns2 <- length(b)
  p <- phyper(overlap - 1, ns1, total - ns1, ns2, lower.tail = FALSE)
  p
}

# Test for each category of gene set and correct P for modules and gene sets
hyper.test.table <- function(l1, l2){ # two lists of gene sets
  pvalue <- t(sapply(names(l1), function(n){
    print(paste0(which(names(l1)==n), ": ", n))
    set <- l1[[n]]
    # Overlap with each module
    sapply(l2, function(mod_genes){
      hyper.test(mod_genes, set, length(ahba.genes()))
    })
  }))
  pvalue <- melt(pvalue)
  pvalue$value <- p.adjust(pvalue$value, method = "BH") # corrected for significant modules and tested gene sets
  m <- dcast(pvalue, Var1 ~ Var2)
  rownames(m) <- m$Var1
  m[,-1]
}

##### Plotting functions #####

module_size <- sapply(modules, length)

prepare.data <- function(m){ # input matrix
  m <- ifelse(m < 0.05, "<0.05", ">=0.05")
  # m <- ifelse(m<0.05, ifelse(m<0.01, "<0.01", "0.01<P<0.05"), ">=0.05")
  rowOrder <- unique(unlist(apply(m, 2, function(x) which(x == "<0.05"))))
  m <- m[rowOrder,]
  df <- lapply(braakModules, function(x) m[,x])
  df <- melt(df)
  colnames(df) <- c("geneset", "module", "value", "dir")
  df$module <- paste0(df$module, " (", module_size[df$module], ")")
  df$module <- factor(df$module, levels = unique(df$module))
  df$geneset <- factor(df$geneset, levels = rev(unique(df$geneset)))
  # df$value <- factor(df$value, levels = unique(df$value)[c(1,3,2)])
  df
}

heat.plot <- function(t) {
  ggplot(t) +
    geom_tile(aes(x=module, y=geneset, fill=value), colour = "black") +
    scale_fill_manual(name = "P-value", values = c("chocolate", "white")) +
    scale_x_discrete(position = "top") +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, size = 10),
          axis.text.y = element_text(size = 10),
          axis.title = element_blank(),
          legend.text = element_text(size = 10), legend.title = element_text(size = 10),
          panel.background = element_blank()
    ) + 
      facet_grid(category~dir, scales = "free", space = "free")
}

########## Module enrichment and plotting table ##########

# Apply hypergeometrix test between gene sets
modEnrich <- lapply(genelists, function(l){ # For each category
  t <- hyper.test.table(l, modules[unlist(braakModules)]) # Apply hypergeometric test
  rows <- apply(t, 1, function(x) any(x < 0.05)) # Filter significant gene sets
  t <- t[rows, ]
  size_row <- sapply(l, length)
  rownames(t) <- paste0(rownames(t), " (", size_row[rownames(t)], ")") # Add gene set size
  rownames(t) <- paste0(toupper(substring(rownames(t),1,1)), substring(rownames(t),2)) # Capital first character
  t
})
save(modEnrich, file = "resources/modEnrich.RData")

# Prepare data for plotting
t1 <- lapply(modEnrich, prepare.data)
t1 <- melt(t1)
colnames(t1)[5] <- "category"
t1$category <- factor(t1$category, levels = unique(t1$category))

# Supplementary heatmap with all GO-terms and diseases
pdf("module_enrichment.pdf", 11, 14)
heat.plot(t1) + theme(legend.position = "top")
dev.off()

##### Find PD-mplicated genes #####

# Genes from studies (count genes)
studies <- list(
  'Jansen et al. 2017' = c("INPP5F", "TMEM175", "ASH1L", "MAPT", "RIT1", "C14orf83", "STK39", "GPNMB", "BST1",
                                         "SIPA1L2", "DLG2", "NUCKS1", "GCH1", "MCCC1", "FAM47E", "BCKDK", "TMPRSS9", "UBOX5",
                                         "CCDC62", "SYNJ1", "EIF4G1", "FBXO7", "C20orf30", "POLG", "VPS13C", "PLA2G6"),
  'Chang et al. 2017' = read.table("chang2017_riskgenes.txt", comment.char = "#", sep = "\n", row.names = NULL, stringsAsFactors = FALSE)[, 1],
   'Nalls et al. 2014' = read.table("nalls2014_riskgenes.txt", comment.char = "#", sep = "\n", row.names = NULL, stringsAsFactors = FALSE)[, 1]#,
                # liscovitch2014 = read.table("ifn_signaling_genes.txt", comment.char = "#", sep = "\n", row.names = NULL, stringsAsFactors = FALSE)[, 1]
)
studies <- lapply(studies, name2EntrezId)
studies <- lapply(studies, function(x) x[!is.na(x)])
studies <- unique(Reduce(c, studies))

# Presence PD genes
t=as.data.frame(sapply(modules[unlist(braakModules)], function(m){
  paste0(entrezId2Name(intersect(studies,m)), collapse = ", ")
}))

##### Write table with overlap and p-value of cell-type enrichment #####

cor <- round(labelCor[names(signif_modules), "r"], digits = 2)
braak_pval <- t(modEnrich$BSGs)
braak_overlap <- sapply(braak, function(set){
  sapply(signif_modules, function(mod_genes){
    length(intersect(mod_genes, set))
  })
})
cell_pval <- t(modEnrich$celltype)
cell_overlap <- sapply(celltype_genes, function(set){
  sapply(signif_modules, function(mod_genes){
    length(intersect(mod_genes, set))
  })
})

table <- cbind('Module' = names(signif_modules), 'r' = cor, 
               braak_overlap, braak_pval,
               cell_overlap, cell_pval)
write.table(table, file = "module_enrichment.txt", sep = "\t", row.names = FALSE)