# Differential expression of eigengenes
# Eigen gene differential expression
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(metafor)
library(reshape2)
library(ggplot2)
library(abind)
library(ggrepel)
load("../ABA_Rdata/BrainExpr.RData")
load("resources/modules.RData")
load("resources/braakInfo.RData")

#####  Functions #####
source("PD/diff.expr.lm.R")

# eigen gene function for data matrix (samples x genes)
eigen.gene <- function(x){
  eg <- prcomp(x, center = FALSE)$x[, 1]# 1st PC (eigen gene expr)
  mean <- rowMeans(x)
  # weight = matrix(abs(cor(x, eg, use = 'p')), nrow(x), ncol(x), byrow = TRUE);
  # eg <- scale(rowMeans(x * weight))
  if (cor(eg, mean) > 0) eg else -eg # flip sign of eigen gene based on the data
}

##### PCA eigen gene #####

# PCA first component of subselection expr. matrices
eigenExpr <- lapply(donorNames, function(d){
  idx <- unlist(braak_idx[[d]]) # idx for samples
  expr <- brainExpr[[d]][, idx]
  as.data.frame(t(sapply(modules, function(genes){ # For each module with genes (grouped gene rows)
    eigen.gene(t(expr[genes, ]))
  })))
})
# save(eigenExpr, file = "resources/eigenExpr.RData")

########### Cell-type reference signals #############

# Cell-type genes
celltypes <- sapply(c("Neurons", "Astrocytes", "Oligodendrocytes", "Microglia", "Endothelial_cells"), function(type){
  file = paste0("brainscope_celltypes/", type, ".txt")
  as.character(read.csv(file, header = TRUE)$entrez_id)
}, simplify = FALSE)

# Cell-type mean gene expression for all samples (whole brain)
mean_celltype <- lapply(donorNames, function(d){
  t(sapply(celltypes, function(ct){
    x <- brainExpr[[d]][ct, ]
    apply(x, 2, mean)
  }))
})

############## Fit linear model & get coefficients ##################
diffExpr <- lapply(donorNames, function(d){
  idx <- unlist(braak_idx[[d]]) # idx for samples
  braak <- paste0("braak", braakLabels[[d]][idx]) # braak labels
  ct <- t(mean_celltype[[d]][, idx]) # samples x cell-types
  expr <- t(eigenExpr[[d]]) # samples x genes
  diff.expr.lm(expr, ct, braak)
})
diffExpr <- simplify2array(diffExpr) # 4D-array: braak x measures x genes x donors

apply(diffExpr, c(1,4), function(b){
  sum(b["BH", ] < 0.05 & abs(b["Estimate",]) > 1)
})

summaryCoef <- apply(diffExpr, 1, function(b){ # For each Braak region
  apply(b, 2, function(g){ # For each gene
    gene <- as.data.frame(t(g))
    summary <- rma(yi = gene$Estimate, vi = gene$`Std. Error`^2, method = "DL", test = "t") # Summary effect size
    gene$weight <- weights(summary)
    rbind(gene, 'summary' = list(summary$beta, summary$se, summary$pval, NA, sum(gene$weight)))
  })
})

# Barplot
summTables <- lapply(summaryCoef, function(b){
  t <- do.call(rbind.data.frame, lapply(b, function(g) g["summary",]))
  t$BH <- p.adjust(t$`Pr(>|t|)`, method = "BH")
  t
})

sapply(summTables, function(t){
  sum(abs(t$Estimate) > 8 & t$BH < 0.01)
})

tab <- summTables$braak6
tab$logp <- -log10(tab$BH)
mod_neg <- rownames(tab)[tab$BH < 0.01 & tab$Estimate < -7] # significant correlated modules
mod_pos <- rownames(tab)[tab$BH < 0.01 & tab$Estimate > 7] # significant correlated modules
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

# xmax <- max(tab$Estimate)+.2
# xmin <- min(tab$Estimate)-.2
# ymax <-  ceiling(max(tab$'logp'))
p <- ggplot(tab, aes(Estimate, logp, colour = info)) +
  geom_point(size = 1, alpha = 0.5) +
  geom_text_repel(aes(label=label), colour = "black", size = 3, nudge_x = 0) +
  scale_colour_manual(values = c("0"="grey", "1"="blue", "2"="red")) +
  labs(y = "-log10 p-value") +
  # scale_x_continuous(limits = c(xmin, xmax), expand = c(0,0)) +
  # scale_y_continuous(limits = c(0, ymax), expand = c(0,0)) +
  theme(legend.position = "none",
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black")
                 # axis.title =  element_text(size = 16),
                 # plot.title = element_text(size = 16),
                 # axis.text = element_text(size = 16)
  )
p
pdf("eigengene_diff_expr_lm.pdf", 6, 3)
print(p)
dev.off()
