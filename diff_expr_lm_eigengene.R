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
  eg <- prcomp(x)$x[, 1]# 1st PC (eigen gene expr)
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

# Cell-type eigengene expression across all samples (whole brain)
ct_ref <- lapply(donorNames, function(d){
  t(sapply(celltypes, function(ct){
    x <- t(brainExpr[[d]][ct, ])
    eg <- prcomp(x)$x[, 1]# 1st PC (eigen gene expr)
    mean <- apply(x, 1, mean)
    if (cor(eg, mean) > 0) eg else -eg # flip sign of eigen gene based on the data
  }))
})

############## Fit linear model & get coefficients ##################
diffExpr <- lapply(donorNames, function(d){
  idx <- unlist(braak_idx[[d]]) # idx for samples
  braak <- paste0("braak", braakLabels[[d]][idx]) # braak labels
  ct <- t(ct_ref[[d]][, idx]) # samples x cell-types
  expr <- t(eigenExpr[[d]]) # samples x genes
  diff.expr.lm(expr, ct, braak)
})
diffExpr <- simplify2array(diffExpr) # 4D-array: braak x measures x genes x donors

apply(diffExpr, c(1,4), function(b){
  sum(b["BH", ] < 0.05 & abs(b["Estimate",]) > 1)
})

summaryCoef <- aaply(diffExpr, c(1,3), function(g){ # For each Braak region and gene
    gene <- as.data.frame(t(g))
    summary <- rma(yi = gene$Estimate, vi = gene$`Std. Error`^2, method = "DL", test = "t") # Summary effect size
    gene$weight <- weights(summary)
    t <- rbind(gene, 'summary' = list(summary$beta, summary$se, summary$pval, NA, sum(gene$weight)))
    as.matrix(t)
})

diffGenes <- apply(summaryCoef2, 1, function(t){
  # rownames(t)[which(abs(t[, "Estimate"])>1 & t[,"BH"] < 0.05)]
  down <-  rownames(t)[which(t[, "Estimate"] < -1 & t[,"BH"] < 0.05)]
  up <- rownames(t)[which(t[, "Estimate"] > 1 & t[,"BH"] < 0.05)]
  list(down = down, up = up)
})
sapply(diffGenes, function(x){
  deg <- sapply(x, length)
  c(deg, sum = sum(deg))
})

# Correct summary P-values
summaryCoef2 <- summaryCoef[,,"summary", -which(dimnames(summaryCoef)[[4]] %in% c("BH", "weight"))] # braak region x genes x measures
summaryCoef2 <- aaply(summaryCoef2, 1, function(t){ # P-value corrected for genes
  b <- p.adjust(t[, "Pr(>|t|)"], method = "BH")
  cbind(t, BH = b)
})

diffGenes <- apply(summaryCoef2, 1, function(t){
  down <-  rownames(t)[which(t[, "Estimate"] < -1 & t[,"BH"] < 0.05)]
  up <- rownames(t)[which(t[, "Estimate"] > 1 & t[,"BH"] < 0.05)]
  list(down = down, up = up)
})
sapply(diffGenes, function(x){
  deg <- sapply(x, length)
  c(deg, sum = sum(deg))
})
tab <- as.data.frame(summaryCoef2["braak6", , ])
tab$logp <- -log10(tab$BH)
mod_down <- rownames(tab)[tab$BH < 0.05 & tab$Estimate < -1] # significant correlated modules
mod_up <- rownames(tab)[tab$BH < 0.05 & tab$Estimate > 1] # significant correlated modules

braakModules2 <- list(down = mod_down, up = mod_up)
save(braakModules2, file = "resources/braakModules2.RData")

tab$info <- sapply(rownames(tab), function(m){
  if (m %in% mod_down) 1
  else if (m %in% mod_up) 2
  else 0
})
tab$label <- rownames(tab)
tab$label[!tab$label %in% c(mod_down, mod_up)] <- ""

# Plotting order of data points 
tab$info <- as.factor(tab$info)
order <- order(tab$info)
tab <- tab[order, ]

ymax <- max(tab$Estimate)+.5
ymin <- min(tab$Estimate)-.5
xmax <-  ceiling(max(tab$'logp'))

p <- ggplot(tab, aes(logp, Estimate, colour = info)) +
  geom_point(size = 2, alpha = 0.5) +
  # geom_text(aes(label=label)) +
  geom_text_repel(aes(label=label), colour = "black", size = 4, nudge_x = 0) +
  scale_colour_manual(values = c("0"="grey", "1"="blue", "2"="red")) +
  labs(x = "-log10 p-value") +
  scale_y_continuous(limits = c(ymin, ymax), expand = c(0,0)) +
  scale_x_continuous(limits = c(0, xmax), expand = c(0,0)) +
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title =  element_text(size = 16),
        plot.title = element_text(size = 16),
        axis.text = element_text(size = 16)
  )
p
pdf("eigengene_diff_expr_lm.pdf", 5, 4.5)
print(p)
dev.off()
