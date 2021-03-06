# Module eigengene expression correlation with Braak stages
library(ggrepel)
library(ComplexHeatmap)

# eigen gene function for data matrix (samples x genes)
eigen.gene <- function(x){
  eg <- prcomp(x, scale. = TRUE)$x[, 1]# 1st PC (eigen gene expr)
  mean <- apply(x, 1, mean)
  if (cor(eg, mean) > 0) eg else -eg # flip sign of eigen gene based on the data
}

##### PCA eigen gene #####

# PCA first component of subselection expr. matrices
eigenExpr <- lapply(donorNames, function(d){
  s <- unlist(braak_idx[[d]])
  expr <- brainExpr[[d]][, s]
  as.data.frame(t(sapply(modules, function(genes){ # For each module with genes (grouped gene rows)
    eigen.gene(t(expr[genes, ]))
  })))
})
saveRDS(eigenExpr, file = "output/eigenExpr.rds")

# Summary Braak correlation
labels <- lapply(donorNames, function(d){
  s <- unlist(braak_idx[[d]])
  braakLabels[[d]][s]
})
summaryLabelCorEG <- summary.braak.cor(eigenExpr, labels)
saveRDS(summaryLabelCorEG, file = "output/summaryLabelCorEG.rds")

labelCor <- do.call(rbind.data.frame, lapply(summaryLabelCorEG, function(g) g["summary",]))
labelCor$BH <- p.adjust(labelCor$pvalue, method = "BH")
orderEG <- order(labelCor$r)
labelCor <- labelCor[orderEG, ]

mod_neg <- rownames(labelCor)[labelCor$BH < 0.0001 & labelCor$r < 0] # significant correlated modules
mod_pos <- rownames(labelCor)[labelCor$BH < 0.0001 & labelCor$r > 0] # significant correlated modules

braakModules <- list('r < 0' = mod_neg, 'r > 0' = mod_pos)

#####  Volcano plot #####

theme <- theme(legend.position = "none",
               panel.background = element_blank(),
               axis.line = element_line(colour = "black"),
               axis.title =  element_text(size = 16),
               plot.title = element_text(size = 16),
               axis.text = element_text(size = 16)
)

tab <- labelCor
tab$'logp' <- -log10(tab$BH)
# eg <- rownames(tab)[tab$BH < 0.001 & tab$r < 0] # significant correlated modules
tab$info <- sapply(rownames(tab), function(m){
  if (m %in% mod_neg) 1
  else if (m %in% mod_pos) 2
  else 0
})
tab$label <- rownames(tab)
tab$label[!tab$label %in% c(mod_neg, mod_pos)] <- ""

# Plotting order of data points 
tab$info <- as.factor(tab$info)
order <- order(tab$info)
tab <- tab[order, ]

xmax <- max(tab$r)+.2
xmin <- min(tab$r)-.2
ymax <-  ceiling(max(tab$'logp'))+.5

p <- ggplot(tab, aes(r, logp, colour = info)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_text_repel(aes(label=label), colour = "black", size = 4, nudge_x = 0) +
  scale_colour_manual(values = c("0"="grey", "1"="blue", "2"="red")) +
  labs(x = bquote("Correlation with Braak "*italic(r)), y = bquote("-log"[10]*" "*italic(P)*"-value")) +
  scale_x_continuous(limits = c(xmin, xmax), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, ymax), expand = c(0,0)) +
  theme

pdf("output/eigengene_r_volcanoplot.pdf", 5, 5)
p
dev.off()

##### Heatmap Expression of eigen gene #####

# Order of modules
cols <- rownames(labelCor)

# Heatmap colors
colPal <- c("darkblue", "white", "darkred")
rampcols <- colorRampPalette(colors = colPal, space="Lab")(201)
colColor <- rampcols[as.numeric(cut(labelCor$r, breaks = 201))]
colsep <- tail(which(labelCor$r < 0), 1)

pdf("output/heatmap_expr_eigengenes.pdf", 5, 12)
lapply(donorNames, function(d){
  samples <- unlist(braak_idx[[d]])
  df <- sample_annot[[d]][samples,]
  labels <- braakLabels[[d]][samples]
  exprMat <- t(as.matrix(eigenExpr[[d]][cols, ]))
  rowsep <- match(c(2:6), labels) -1# separate Braak regions
  rowColor <- df$color_hex_triplet
  heatmap.2(exprMat, col = rampcols, 
            labRow = df$name, 
            rowsep = rowsep, colsep = colsep, sepcolor = "black",
            Rowv=FALSE, Colv=FALSE, 
            cexCol = 1, cexRow = .1,
            scale = "none", trace = "none", dendrogram = "none", key = FALSE, 
            RowSideColors = rowColor, ColSideColors = colColor,
            main = d,
            margins = c(20,5))
})
dev.off()

##### Print table with module info #####

module_info <- data.frame(
  # Module = rownames(labelCor),
  Size = sapply(modules, length)[rownames(labelCor)], 
  r = round(labelCor$r, digits = 2),
  "BH-corrected P" = format(labelCor$pvalue, digits = 3, scientific = TRUE),
  genes = sapply(modules, function(m) paste0(entrezId2Name(m), collapse = ","))[rownames(labelCor)]
)
write.table(module_info, file = "output/module_info.txt", sep = "\t", quote = FALSE)

##### Heatmap of eigengene co-expression #####

pdf("output/heatmap_coexpr_modules.pdf", 3.2, 2.6)
lapply(donorNames, function(d){
  lapply(names(braakModules), function(r){
    m <- braakModules[[r]]
    e <- eigenExpr[[d]][m, ]
    cor <- cor(t(e)) # correlation between eigengenes (diagonal matrix)
    Heatmap(cor, name = "r",
            column_title = paste0(gsub("donor", "Donor ", d), ", ", r),
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            row_names_gp = gpar(fontsize = 8),
            column_names_gp = gpar(fontsize = 8),
            column_names_rot = 45,
            width = unit(ncol(cor)*.8, "lines"),
            height = unit(nrow(cor)*.8, "lines")
    )
  })
})
dev.off()


# # Module membership (correlation with eigengene) to identify hub genes
# module_membership <- sapply(unname(unlist(braakModules)), function(m){
#   mm <- sapply(donorNames, function(d){
#     eg <- unlist(eigenExpr[[d]][m, ]) # eigengene expression of module
#     s <- unlist(braak_idx[[d]]) # sample indices for all braak samples
#     module_genes <- modules[[m]] # genes in module
#     gene_expr <- brainExpr[[d]][module_genes, s] # gene expression in all Braak regions
#     cor(t(gene_expr), eg)[,1] # module membership
#     # sort(mm)
#   })
#   avg_mm <- apply(mm, 1, mean) # module membership averaged across donors
#   rev(sort(avg_mm))
# }, simplify = FALSE)
# 
# hubgenes <- sapply(module_membership, function(m){
#   hg <- names(m)[1:10]
#   pd_hg <- lapply(pdGenesID, function(x)intersect(hg, x))
#   sapply(pd_hg, length)
# })
