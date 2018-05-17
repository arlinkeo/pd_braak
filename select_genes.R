# select significant genes based on significant summary estimate
setwd("C:/Users/dkeo/surfdrive/pd_braak")
library("metap")
library(gplots)
library(ggplot2)
library(reshape2)
source("PD/base_script.R")
load("resources/summaryDiffExpr.RData")
load("resources/summaryLabelCorr.RData")

# Select region pair
diffExpr <- summaryDiffExpr$`braak1-braak6`

#Filter for summary effect
diffExpr <- do.call(rbind.data.frame, lapply(diffExpr, function(g) g["summary",]))
labelCor <- do.call(rbind.data.frame, lapply(summaryLabelCorr, function(g) g["summary",]))

# number of genes in top 10%
top10 <- floor(nrow(labelCor)*0.1)

# Get top 10% correlated genes
order <- rev(order(abs(labelCor$r))) # order absolute corr.
corrGenes <- rownames(labelCor)[order[1:top10]] # top 10% genes
corrGenes_neg <- corrGenes[labelCor[corrGenes, "r"] < 0]
corrGenes_pos <- corrGenes[labelCor[corrGenes, "r"] > 0]

# Top 10% mean Diff. expressed genes braak 1 vs. braak 6
order <- rev(order(abs(diffExpr$meanDiff))) # order absolute corr.
diffGenes2 <- rownames(diffExpr)[order[1:top10]] # top 10% genes
diffGenes2_up <- diffGenes2[diffExpr[diffGenes2, "meanDiff"] > 0]
diffGenes2_down <- diffGenes2[diffExpr[diffGenes2, "meanDiff"] < 0]

# Top 10% significant Diff. expressed genes braak 1 vs. braak 6
diffExpr$BH <- p.adjust(diffExpr$pvalue, method = "BH")
order <- order(diffExpr$BH) # order absolute corr.
diffGenes1 <- rownames(diffExpr)[order[1:top10]] # top 10% genes
diffGenes1_up <- diffGenes1[diffExpr[diffGenes1, "meanDiff"] > 0]
diffGenes1_down <- diffGenes1[diffExpr[diffGenes1, "meanDiff"] < 0]

# Venn
criteria <- c("r", "fc", "pval_fc")
ll <- list(corrGenes, diffGenes2, diffGenes1)
ll1 <- list(corrGenes_neg, diffGenes2_up, diffGenes1_up)
ll2 <- list(corrGenes_pos, diffGenes2_down, diffGenes1_down)
names(ll) <- criteria
names(ll1) <- criteria
names(ll2) <- criteria
venn <- venn(ll)
venn1 <- venn(ll1)
venn2 <- venn(ll2)

braakGenes <- attr(venn, "intersections")[[paste0(criteria, collapse = ":")]]
names <- entrezId2Name(braakGenes)
braakGenes <- data.frame(entrez_id = braakGenes, gene_symbol = names)
braakGenes$braak_r <- labelCor[braakGenes$entrez_id, "r"]
braakGenes$FC <- diffExpr[braakGenes$entrez_id, "meanDiff"]
braakGenes$BH <- diffExpr[braakGenes$entrez_id, "BH"]
braakGenes$dir <- ifelse(braakGenes$braak_r>0, "pos", "neg")
braakGenes <- braakGenes[order(braakGenes$braak_r), ]
write.table(braakGenes, file = "braakGenes.txt", row.names = FALSE, quote = FALSE, sep = "\t")

# positive_r <- braakGenes[labelCor[braakGenes, "r"] > 0]
# negative_r <- braakGenes[labelCor[braakGenes, "r"] < 0]
# braakGenes <- list(positive_r = positive_r, negative_r = negative_r)
save(braakGenes, file = "resources/braakGenes.RData")

#thresholds
t1 <- paste("|r| >", floor(min(abs(braakGenes$braak_r))*100)/100)
t2 <- paste("|FC| >", floor(min(abs(braakGenes$FC))*100)/100)
t3 <- paste("Pfc <", floor(max(braakGenes$BH)*1e5)/1e5)
labels <- c(t1, t2, t3)

# Bar plot of selected genes
tab <- cbind(negative = -sapply(ll1, length), positive = sapply(ll2, length))
rownames(tab) <- labels
tab <- rbind(tab, Overlap = c(-sum(braakGenes$dir=="neg"), sum(braakGenes$dir=="pos")))

tab <- melt(tab)
colnames(tab) <- c("type", "dir", "size")
tab$type <- factor(tab$type, levels = rev(unique(tab$type)))
tab$y <- ifelse(tab$dir == "positive", tab$size+200, tab$size-300)

p <- ggplot(tab) + 
  geom_col(aes(x=type, y = size, fill=dir), size = 0.5, colour = "black") + 
  geom_text(aes(x=type, y=y, label=format(abs(tab$size), big.mark=","))) + 
  scale_fill_manual(values = c("blue", "red")) +
  scale_y_continuous(expand = c(0.1,0.1)) +
  coord_flip() +
  labs(x = "", y = "Number of selected genes") +
  theme(
    axis.text = element_text(size = 11),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    legend.title = element_blank()
  )
p
pdf("braakgenes_barplot.pdf", 5, 2)
print(p)
dev.off()
# 
# #Presence of PD-implicated genes
# pd.genes <- function(x){
#   lapply(pdGenesID, function(l){
#    intersect(l, x)
#   })
# }
# 
# pd <- pd.genes(unlist(braakGenes))
# lapply(pd, entrezId2Name)
# 
# # Print braak genes with correlations, mean diff, and p-values
# gene.stats <- function(g){
#   r <- labelCor[g, "r", drop = FALSE] 
#   geneorder <- rownames(r)[order(r)]
#   r <- r[geneorder, ]
#   r <- round(r, digits = 2)
#   de <- diffExpr[geneorder, c("meanDiff", "benjamini_hochberg")]
#   de$meanDiff <- round(de$meanDiff, digits = 2)
#   de$benjamini_hochberg <- format(de$benjamini_hochberg, digits = 3, scientific = TRUE)
#   id <- geneorder
#   name <- entrezId2Name(geneorder)
#   df <- data.frame(id, name, r, de)
#   df
# }
# 
# tab <- Reduce(rbind, lapply(pd, gene.stats))
# write.table(tab, file = "pdgenes_stats.txt", sep ="\t", quote = FALSE, row.names = FALSE)