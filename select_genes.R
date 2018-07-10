# select significant genes based on significant summary estimate
setwd("C:/Users/dkeo/surfdrive/pd_braak")
library("metap")
library(gplots)
library(ggplot2)
library(reshape2)
library(biomaRt)
source("PD/base_script.R")
load("resources/summaryDiffExpr.RData")
load("resources/summaryLabelCorr.RData")

# Select region pair
diffExpr <- summaryDiffExpr$`braak1-braak6`

#Filter for summary effect
diffExpr <- do.call(rbind.data.frame, lapply(diffExpr, function(g) g["summary",]))
labelCor <- do.call(rbind.data.frame, lapply(summaryLabelCorr, function(g) g["summary",]))

### Selecting Braak-related genes ###

top10 <- floor(nrow(labelCor)*0.1) # number of genes in top 10%

# Get top 10% correlated genes
order <- rev(order(abs(labelCor$r))) # order absolute corr.
corrGenes <- rownames(labelCor)[order[1:top10]] # top 10% genes
corrGenes_neg <- corrGenes[labelCor[corrGenes, "r"] < 0]
corrGenes_pos <- corrGenes[labelCor[corrGenes, "r"] > 0]

# Top 10% mean Diff. expressed genes braak 1 vs. braak 6
order <- rev(order(abs(diffExpr$meanDiff))) # order absolute corr.
diffGenes1 <- rownames(diffExpr)[order[1:top10]] # top 10% genes
diffGenes1_pos <- diffGenes1[diffExpr[diffGenes1, "meanDiff"] > 0]
diffGenes1_neg <- diffGenes1[diffExpr[diffGenes1, "meanDiff"] < 0]

# Top 10% significant Diff. expressed genes braak 1 vs. braak 6
diffExpr$BH <- p.adjust(diffExpr$pvalue, method = "BH")
order <- order(diffExpr$BH) # order absolute corr.
diffGenes2 <- rownames(diffExpr)[order[1:top10]] # top 10% genes
diffGenes2_pos <- diffGenes2[diffExpr[diffGenes2, "meanDiff"] > 0]
diffGenes2_neg <- diffGenes2[diffExpr[diffGenes2, "meanDiff"] < 0]

# Venn diagram to visualize overlap of selections
criteria <- c("r", "fc", "pval_fc")
ll <- list(corrGenes, diffGenes1, diffGenes2)
ll1 <- list(corrGenes_neg, diffGenes2_neg, diffGenes1_neg)
ll2 <- list(corrGenes_pos, diffGenes2_pos, diffGenes1_pos)
names(ll) <- criteria
names(ll1) <- criteria
names(ll2) <- criteria
venn <- venn(ll)
venn1 <- venn(ll1)
venn2 <- venn(ll2)

# Selection based on 3 criteria and collect info of selected genes
braakGenes <- attr(venn, "intersections")[[paste0(criteria, collapse = ":")]]
braakGenes <- data.frame(entrez_id = braakGenes, gene_symbol = entrezId2Name(braakGenes))
braakGenes$braak_r <- round(labelCor[braakGenes$entrez_id, "r"], digits = 2)
braakGenes$FC <- round(diffExpr[braakGenes$entrez_id, "meanDiff"], digits = 2)
braakGenes$BH <- format(diffExpr[braakGenes$entrez_id, "BH"], digits = 2)
braakGenes$dir <- ifelse(braakGenes$braak_r>0, "pos", "neg")
write.table(braakGenes, file = "braakGenes.txt", row.names = FALSE, quote = FALSE, sep = "\t")
save(braakGenes, file = "resources/braakGenes.RData")

# Bar plot of selected genes
t1 <- paste("|r| >", floor(min(abs(braakGenes$braak_r))*100)/100)
t2 <- paste("|FC| >", floor(min(abs(braakGenes$FC))*100)/100)
t3 <- paste("Pfc <", floor(max(braakGenes$BH)*1e5)/1e5)
labels <- c(t1, t2, t3)
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
 
#Presence of PD-implicated genes
parkingenes <- read.table("parkin_genes_hgnc.txt", sep = "\t", header = TRUE)$Approved.Symbol

pdGenes <- list(hiImpact = c("SNCA", "LRRK2", "GBA", "VPS35", "PARK2", "UCHL1", "PINK1", "PARK7", "ATP13A2", "PLA2G6", "FBXO7", "DNAJC6", "SYNJ1", 
                             "EIF4G1", "DNAJC13", "CHCHD2", "C20orf30", "RIC3", "LRP10"), #TMEM230 is C20orf30
                jansen2017 = c("INPP5F", "TMEM175", "ASH1L", "MAPT", "RIT1", "C14orf83", "STK39", "GPNMB", "BST1", 
                               "SIPA1L2", "DLG2", "NUCKS1", "GCH1", "MCCC1", "FAM47E", "BCKDK", "TMPRSS9", "UBOX5", 
                               "CCDC62", "SYNJ1", "EIF4G1", "FBXO7", "C20orf30", "POLG", "VPS13C", "PLA2G6"),
                hla = c("HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DQB1"),
                'Chang et al. 2017' = read.table("chang2017_riskgenes.txt", comment.char = "#", sep = "\n", row.names = NULL, stringsAsFactors = FALSE)[, 1], 
                'Nalls et al. 2014' = read.table("nalls2014_riskgenes.txt", comment.char = "#", sep = "\n", row.names = NULL, stringsAsFactors = FALSE)[, 1]
)
pdGenesID <- lapply(pdGenes, name2EntrezId)
tab <- lapply(names(pdGenesID), function(n){
  x <- pdGenesID[[n]]
  g <- intersect(braakGenes$entrez_id, x)
  cbind(study = rep(n, length(g)), braakGenes[braakGenes$entrez_id %in% g,])
})
tab <- Reduce(rbind, tab)
write.table(tab, file = "pdgenes_stats.txt", sep ="\t", quote = FALSE, row.names = FALSE)