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
mean_celltype <- lapply(donorNames, function(d){
  t(sapply(celltypes, function(ct){
    x <- brainExpr[[d]][ct, ]
    apply(x, 2, mean)
  }))
})

# Fit linear model & get coefficients
diffExpr <- lapply(donorNames, function(d){
  idx <- unlist(braak_idx[[d]]) # idx for samples
  braak <- paste0("braak", braakLabels[[d]][idx]) # braak labels
  ct <- t(mean_celltype[[d]][, idx]) # samples x cell-types
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
saveRDS(summaryCoef, file = "resources/summaryCoef_mean.rds")
  
# Barplot

summTables <- lapply(summaryCoef, function(b){
  t <- do.call(rbind.data.frame, lapply(b, function(g) g["summary",]))
  t$BH <- p.adjust(t$`Pr(>|t|)`, method = "BH")
  t
})

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