# Line plot of dopaminergic genes

# Plot function
plot.matrix <- function(expr, title){
  df <- melt(as.matrix(expr))
  colnames(df) <- c("gene", "sample", "expr")
  df$gene <- factor(df$gene, levels = unique(df$gene))
  df$sample <- factor(df$sample, levels = unique(df$sample))
  
  ggplot(df, aes(x=sample, y=expr, group=gene, color=gene)) + 
    geom_point() +
    geom_line() +
    # geom_smooth(aes(group=gene), span = 0.1) +
    # geom_vline(xintercept = vline) +
    ggtitle(title) +
    theme_classic()
}

# Function get expression of genes and modules
get.expr <- function(g = NULL, m = NULL){
  lapply(donorNames, function(d){
    idx <- braak_idx[[d]]
    labels <- braakLabels[[d]][unlist(idx)]
    idx_eg <- sapply(names(idx), function(i){ #index of sorted labels
      which(labels == gsub("R","", i))
    }, simplify = FALSE) 
    e <- sapply(names(idx), function(i){
      if (is.null(g)) e <- integer(0) else e <- brainExpr[[d]][g, idx[[i]]] # expr of genes
      if (is.null(m)) em <- integer(0) else  em <- eigenExpr[[d]][m, idx_eg[[i]]]# expr of module
      e <- rbind(em, e)
      apply(e, 1, median)
    })
    rownames(e) <- c(m, entrezId2Name(g))
    e
  })
}

# Data of dopaminergic genes & module
genes <- c("SNCA", "GCH1", "TH", "SLC6A3", "SLC18A2")
genes_id <- name2EntrezId(genes)
module <- c("M127")
expr_dopa <- get.expr(genes_id, module)

pdf("lineplot_dopaminergic_genes.pdf", 6, 4)
# # Per donor
# lapply(donorNames, function(d){
#   e <- expr[[d]]
#   plot.matrix(e, d)
# })
# Average median across donor
mean_median <- apply(simplify2array(expr_dopa), c(1,2), mean)
plot.matrix(mean_median, "mean median across donors")
# median_median <- apply(simplify2array(expr_dopa), c(1,2), median)
# plot.matrix(median_median, "median median across donors")
# max_median <- apply(simplify2array(expr_dopa), c(1,2), max)
# plot.matrix(max_median, "max median across donors")
# min_median <- apply(simplify2array(expr_dopa), c(1,2), min)
# plot.matrix(min_median, "min median across donors")
dev.off()

# Data of blood oxygen genes
genes <- c("HBD", "HBB", "HBA1", "HBA2", "OASL")
genes_id <- name2EntrezId(genes)
module <- c("M47")
expr_bloodoxy <- get.expr(genes_id, module)

pdf("output/lineplot_bloodoxygen_genes.pdf", 6, 4)
mean_median <- apply(simplify2array(expr_bloodoxy), c(1,2), mean)
plot.matrix(mean_median, "mean median across donors")
dev.off()

##### Line plot of module eigengenes #####

pdf("output/lineplot_braakmodules.pdf", 6, 4)
expr_braakmodules <- get.expr(m = unlist(braakModules))
mean_median <- apply(simplify2array(expr_braakmodules), c(1,2), mean)
plot.matrix(mean_median, "mean median across donors")
dev.off()

# lapply(donorNames[1], function(d){
#   e <- as.matrix(eigenExpr[[d]][unlist(braakModules), ])
#   df <- melt(e)
#   colnames(df) <- c("module", "sample", "expr")
#   df$module <- factor(df$module, levels = unique(df$module))
#   df$sample <- factor(df$sample, levels = unique(df$sample))
#   ggplot(df, aes(x=sample, y=expr, group=module, color=module)) +
#     geom_smooth(aes(group=module), span = 0.1, size = 0.5) +
#     theme_classic()
# })