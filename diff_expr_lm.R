#  Differential expression based on linear regression
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(metafor)
library(reshape2)
library(ggplot2)
library(abind)
# library(plyr)
load("resources/braakInfo.RData")
load("../ABA_Rdata/BrainExpr.RData")

# Load function
source("PD/diff.expr.lm.R")

# Cell-type genes
celltypes <- sapply(c("Neurons", "Astrocytes", "Oligodendrocytes", "Microglia", "Endothelial_cells"), function(type){
  file = paste0("brainscope_celltypes/", type, ".txt")
  as.character(read.csv(file, header = TRUE)$entrez_id)
}, simplify = FALSE)

# Cell-type mean gene expression for all samples (whole brain)
pca_celltype <- lapply(donorNames, function(d){
  t(sapply(celltypes, function(ct){
    x <- t(brainExpr[[d]][ct, ])
    eg <- prcomp(x)$x[, 1]# 1st PC (eigen gene expr)
    mean <- apply(x, 1, mean)
    if (cor(eg, mean) > 0) eg else -eg # flip sign of eigen gene based on the data
  }))
})

# Fit linear model & get coefficients
diffExpr <- lapply(donorNames, function(d){
  idx <- unlist(braak_idx[[d]]) # idx for samples
  braak <- paste0("braak", braakLabels[[d]][idx]) # braak labels
  ct <- t(pca_celltype[[d]][, idx]) # samples x cell-types
  expr <- t(brainExpr[[d]][, idx]) # samples x genes
  diff.expr.lm(expr, ct, braak)
})
diffExpr <- simplify2array(diffExpr) # 4D-array: braak x measures x genes x donors

apply(diffExpr, c(1,4), function(b){
  sum(b["BH", ] < 0.05 & abs(b["Estimate",]) > 1)
})

# Meta-analysis of mean expression difference across donors
summaryCoef <- apply(diffExpr, 1, function(b){ # For each Braak region
  apply(b, 2, function(g){ # For each gene
    gene <- as.data.frame(t(g))
    summary <- rma(yi = gene$Estimate, vi = gene$`Std. Error`^2, method = "DL", test = "t") # Summary effect size
    gene$weight <- weights(summary)
    rbind(gene, 'summary' = list(summary$beta, summary$se, summary$pval, NA, sum(gene$weight)))
  })
})
saveRDS(summaryCoef, file = "resources/summaryCoef_pca.rds")
  
# Barplot
summaryCoef_pca <- readRDS("resources/summaryCoef_pca.rds")
summaryCoef_mean <- readRDS("resources/summaryCoef_mean.rds")
summaryCoef_svd <- readRDS("resources/summaryCoef_svd.rds")

summ.table <- function(x){
  lapply(x, function(b){
    t <- do.call(rbind.data.frame, lapply(b, function(g) g["summary",]))
    t$BH <- p.adjust(t$`Pr(>|t|)`, method = "BH")
    t
  })
}

summTables_pca <- summ.table(summaryCoef_pca)
summTables_mean <- summ.table(summaryCoef_mean)
summTables_svd <- summ.table(summaryCoef_svd)
mean_b6 <- summTables_mean$braak6
pca_b6 <- summTables_pca$braak6
svd_b6 <- summTables_svd$braak6

b <- "braak6"

# venn diagram overlap between methods
load("resources/not_celltype_corrected/braakGenes.RData")
brg_neg <- braakGenes$entrez_id[braakGenes$braak_r<0]
brg_pos <- braakGenes$entrez_id[braakGenes$braak_r>0]

find_diffgenes <- function(x, pos = TRUE){
  if (pos) rownames(x)[x$Estimate > 1 & x$BH < 0.05]
  else rownames(x)[x$Estimate < -1 & x$BH < 0.05]
}
mean_neg <- find_diffgenes(mean_b6, pos = FALSE)
mean_pos <- find_diffgenes(mean_b6)
pca_neg <- find_diffgenes(pca_b6, pos = FALSE)
pca_pos <- find_diffgenes(pca_b6)
svd_neg <- find_diffgenes(svd_b6, pos = FALSE)
svd_pos <- find_diffgenes(svd_b6)

library(venn)
ll1 <- list(brg = brg_neg, mean = mean_neg, pca = pca_neg, svd = svd_neg)
ll2 <- list(brg = brg_pos, mean = mean_pos, pca = pca_pos, svd = svd_pos)
pdf("overlap_results_brg_eq2.pdf", 2, 2)
venn(ll1[c(1,2)]); title("negative", line = -1)
venn(ll1[c(1,3)]); title("negative", line = -1)
venn(ll1[c(1,4)]); title("negative", line = -1)
venn(ll2[c(1,2)]); title("positive", line = -1)
venn(ll2[c(1,3)]); title("positive", line = -1)
venn(ll2[c(1,4)]); title("positive", line = -1)
dev.off()


# Fold-change of marker genes
arr <- simplify2array(lapply(summTables, as.matrix))
fc <- arr[celltypes$Microglia, "Estimate",]
fc <- fc[order(fc[,5]),]
df <- melt(fc)
colnames(df) <- c("gene", "braak", "estimate")
df$gene <- factor(df$gene, levels = unique(df$gene))
lim <- max(abs(df$estimate))
ggplot(df, aes(x=braak, y=gene, fill = estimate)) + 
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                                    midpoint = 0, limit = c(-lim, lim), space = "Lab", 
                                    name="estimate") 

degLists <- t(sapply(summTables, function(t){
  negative_r <- -sum(t$Estimate < -1 & t$BH < 0.05)
  positive_r <- sum(t$Estimate > 1 & t$BH < 0.05)
  c('negative' = negative_r, 'positive' = positive_r)
}))

df <- melt(degLists)
colnames(df) <- c("region", "dir", "deg")
df$region <- factor(df$region, levels = rev(unique(df$region)))
df$y <- ifelse(df$dir == "positive", df$deg+500, df$deg-500)

p <- ggplot(df) + 
  geom_col(aes(x=region, y = deg, fill=dir), size = 0.5, colour = "black") + 
  geom_text(aes(x=region, y= y, label=format(abs(df$deg), big.mark=","))) + 
  scale_fill_manual(values = c("blue", "red")) +
  scale_y_continuous(expand = c(0.1,0.1)) +
  coord_flip() +
  labs(x = "", y = "Number of differentially expressed genes") +
  theme(
    axis.text = element_text(size = 11),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    legend.title = element_blank()
  )
p
pdf("diff_expr_lm_eg_barplot.pdf", 6, 3)
print(p)
dev.off()