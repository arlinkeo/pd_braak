# Heatmap of modules enriched for different gene sets
# Braak gene enrichment of modules
library(GO.db)

########## Prepare list of Braak genes, cell types, GO-term, and diseases ##########

# Genes associated with GO-terms
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=92)
goterms <- unlist(getBM(c('go_id'), mart = ensembl, values = 'go'))
names(goterms) <- goterms
go_genes <- lapply(goterms, function(go){
  as.character(unlist(getBM(c('entrezgene'), filters='go', mart=ensembl, values=go)))
})
saveRDS(go_genes, file = "output/go_genes.rds")
names(go_genes) <- Term(names(go_genes))
go_genes <- go_genes[which(!is.na(names(go_genes)))]

# Disease-associated genes
disease_table <- read.delim("../gene_associated_with_disease.txt")
diseases <- unique(disease_table$diseaseName)
disease_genes <- sapply(diseases, function(x){
  disease_table$geneId[disease_table$diseaseName %in% x]
}, simplify = FALSE)

# Cell-type genes
celltype_genes <- sapply(c("Astrocytes", "Endothelial_cells", "Microglia", "Neurons", "Oligodendrocytes"), function(type){
  file = paste0("../brainscope_celltypes/", type, ".txt")
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
saveRDS(genelists, file = "output/genelists.rds")
    
########## Module enrichment functions ##########

# Test for each category of gene set and correct P for modules AND gene sets
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
  rowOrder <- unique(unlist(apply(m, 2, function(x) which(x == "<0.05"))))
  m <- m[rowOrder,]
  df <- lapply(braakModules, function(x) {
    sub_m <- m[,x]
    colnames(sub_m) <- paste0(colnames(sub_m), " (", module_size[colnames(sub_m)], ")")
    sub_m
  })
  df <- melt(df)
  colnames(df) <- c("geneset", "module", "value", "dir")
  df$module <- factor(df$module, levels = unique(df$module))
  df$geneset <- factor(df$geneset, levels = rev(unique(df$geneset)))
  df
}

heat.plot <- function(t) {
  ggplot(t) +
    geom_tile(aes(x=module, y=geneset, fill=value), colour = "black") +
    scale_fill_manual(name = "P-value", values = c("chocolate", "white")) +
    scale_x_discrete(position = "top") +
    # scale_y_discrete(labels = function(x) paste(strwrap(x, width = 100), collapse = "\n")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, size = 10),
          axis.text.y = element_text(size = 10),
          axis.title = element_blank(),
          legend.text = element_text(size = 10), legend.title = element_text(size = 10),
          panel.background = element_blank()
    ) + 
      facet_grid(category~dir, scales = "free", space = "free")
}

########## Module enrichment and plotting table ##########

# Apply hypergeometric test between gene sets
modEnrich <- lapply(genelists, function(l){ # For each category
  t <- hyper.test.table(l, modules[unlist(braakModules)]) # Apply hypergeometric test
  rows <- apply(t, 1, function(x) any(x < 0.05)) # Filter significant gene sets
  t <- t[rows, ]
  size_row <- sapply(l, length)
  rownames(t) <- paste0(rownames(t), " (", size_row[rownames(t)], ")") # Add gene set size
  rownames(t) <- paste0(toupper(substring(rownames(t),1,1)), substring(rownames(t),2)) # Capital first character
  t
})
saveRDS(modEnrich, file = "output/modEnrich.rds")

# Poster version
rowOrder <- unique(unlist(apply(modEnrich$GO, 2, function(x) which(x< 0.05))))
modEnrich$GO <- modEnrich$GO[rowOrder, ]
modEnrich$GO <- modEnrich$GO[-c(4:13,17,19:20,25:28,30,35,36,38,40,43:47,50), ]
rowOrder <- unique(unlist(apply(modEnrich$disease, 2, function(x) which(x< 0.05))))
modEnrich$disease <- modEnrich$disease[rowOrder, ]
modEnrich$disease <- modEnrich$disease[-c(4,9:11:15,20,21,24), ]

# Prepare data for plotting
t <- lapply(modEnrich, prepare.data)
t <- melt(t)
colnames(t)[5] <- "category"
t$category <- factor(t$category, levels = unique(t$category))

# Supplementary heatmap with all GO-terms and diseases
pdf("module_enrichment.pdf", 11, 14)
# pdf("module_enrichment_posterversion.pdf", 7.5, 9)
heat.plot(t) + theme(legend.position = "top")
dev.off()

# ##### Write table with overlap and p-value of cell-type enrichment #####
# 
# cor <- round(labelCor[names(signif_modules), "r"], digits = 2)
# braak_pval <- t(modEnrich$BSGs)
# braak_overlap <- sapply(braak, function(set){
#   sapply(signif_modules, function(mod_genes){
#     length(intersect(mod_genes, set))
#   })
# })
# cell_pval <- t(modEnrich$celltype)
# cell_overlap <- sapply(celltype_genes, function(set){
#   sapply(signif_modules, function(mod_genes){
#     length(intersect(mod_genes, set))
#   })
# })
# 
# table <- cbind('Module' = names(signif_modules), 'r' = cor, 
#                braak_overlap, braak_pval,
#                cell_overlap, cell_pval)
# write.table(table, file = "module_enrichment.txt", sep = "\t", row.names = FALSE)