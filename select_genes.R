# select significant genes based on significant summary estimate
setwd("C:/Users/dkeo/surfdrive/pd_braak")
library("metap")
library(gplots)
library(ggplot2)
library(reshape2)
source("PD/base_script.R")
load("resources/summaryDiffExpr.RData")
load("resources/summaryLabelCor.RData")

########## Prepare data ##########

# Select region pair
diffExpr <- summaryDiffExpr$`braak1-braak6`

#Filter for summary effect
diffExpr <- do.call(rbind.data.frame, lapply(diffExpr, function(g) g["summary",]))
diffExpr$BH <- p.adjust(diffExpr$pvalue, method = "BH")
labelCor <- do.call(rbind.data.frame, lapply(summaryLabelCor, function(g) g["summary",]))
labelCor$BH <- p.adjust(labelCor$pvalue, method = "BH")

########## Selecting Braak-related genes ##########

# top10 <- floor(nrow(labelCor)*0.1) # number of genes in top 10%
# 
# # Get top 10% correlated genes
# order <- rev(order(abs(labelCor$r))) # order absolute corr.
# corrGenes <- rownames(labelCor)[order[1:top10]] # top 10% genes
corrGenes <- rownames(labelCor)[labelCor$BH < 0.05 & abs(labelCor$r) > 0.6]
corrGenes_neg <- corrGenes[labelCor[corrGenes, "r"] < 0]
corrGenes_pos <- corrGenes[labelCor[corrGenes, "r"] > 0]
max(abs(labelCor[corrGenes, "r"]))
min(abs(labelCor[corrGenes, "r"]))

# Top 10% mean Diff. expressed genes braak 1 vs. braak 6
# order <- rev(order(abs(diffExpr$meanDiff))) # order absolute corr.
# diffGenes1 <- rownames(diffExpr)[order[1:top10]] # top 10% genes
# diffGenes1_pos <- diffGenes1[diffExpr[diffGenes1, "meanDiff"] > 0]
# diffGenes1_neg <- diffGenes1[diffExpr[diffGenes1, "meanDiff"] < 0]
# max(abs(diffExpr[diffGenes1, "meanDiff"]))
# min(abs(diffExpr[diffGenes1, "meanDiff"]))
# 
# # Top 10% significant Diff. expressed genes braak 1 vs. braak 6
# diffExpr$BH <- p.adjust(diffExpr$pvalue, method = "BH")
# order <- order(diffExpr$BH) # order absolute corr.
# diffGenes2 <- rownames(diffExpr)[order[1:top10]] # top 10% genes
# diffGenes2_pos <- diffGenes2[diffExpr[diffGenes2, "meanDiff"] > 0]
# diffGenes2_neg <- diffGenes2[diffExpr[diffGenes2, "meanDiff"] < 0]

diffGenes <- rownames(diffExpr)[diffExpr$BH < 0.05 & abs(diffExpr$meanDiff) > 1] # top 10% genes
diffGenes_pos <- diffGenes[diffExpr[diffGenes, "meanDiff"] > 0]
diffGenes_neg <- diffGenes[diffExpr[diffGenes, "meanDiff"] < 0]

# Venn diagram to visualize overlap of selections
criteria <- c("r", "fc")#, "pval_fc")
ll <- list(corrGenes, diffGenes)#diffGenes1, diffGenes2)
ll1 <- list(corrGenes_neg, diffGenes_neg)
ll2 <- list(corrGenes_pos, diffGenes_pos)
# ll1 <- list(corrGenes_neg, diffGenes2_neg, diffGenes1_neg)
# ll2 <- list(corrGenes_pos, diffGenes2_pos, diffGenes1_pos)
names(ll) <- criteria
names(ll1) <- criteria
names(ll2) <- criteria
venn <- venn(ll)
venn1 <- venn(ll1)
venn2 <- venn(ll2)

# Selection based on 2 criteria and collect info of selected genes
braakGenes <- attr(venn, "intersections")[[paste0(criteria, collapse = ":")]]
braakGenes <- data.frame(entrez_id = braakGenes, gene_symbol = entrezId2Name(braakGenes))
braakGenes$r <- round(labelCor[braakGenes$entrez_id, "r"], digits = 2)
braakGenes$r_BH <- format(labelCor[braakGenes$entrez_id, "BH"], digits = 2)
braakGenes$FC <- round(diffExpr[braakGenes$entrez_id, "meanDiff"], digits = 2)
braakGenes$FC_BH <- format(diffExpr[braakGenes$entrez_id, "BH"], digits = 2)
# braakGenes$dir <- ifelse(braakGenes$braak_r>0, "pos", "neg")
write.table(braakGenes, file = "braakGenes.txt", row.names = FALSE, quote = FALSE, sep = "\t")
save(braakGenes, file = "resources/braakGenes.RData")

########## Bar plot of selected genes ##########

# row labels
t1 <- paste("|r| >", floor(min(abs(braakGenes$r))*100)/100, "& P < 0.05")
t2 <- paste("|FC| >", floor(min(abs(braakGenes$FC))*100)/100, "& P < 0.05")
# t3 <- paste("Pfc <", floor(as.numeric(max(braakGenes$FC_BH))*1e5)/1e5)
labels <- c(t1, t2)#, t3)
tab <- cbind(negative = -sapply(ll1, length), positive = sapply(ll2, length))
rownames(tab) <- labels
tab <- rbind(tab, Overlap = c(-sum(braakGenes$r<0), sum(braakGenes$r>0)))

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
 
########## Presence of PD-implicated genes ##########
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

########## Volcano plot for differential expression analysis ##########
# library(gridExtra)
# library(ggrepel)

# Extract summary statistics, BH-correct and log-transformed p-values
summaryDiffExpr <- lapply(summaryDiffExpr, function(rp){
  tab <- do.call(rbind.data.frame, lapply(rp, function(g) g["summary",]))
  tab$BH <- p.adjust(tab$pvalue, method = "BH")
  tab$'logp' <- -log10(tab$BH)
  tab
})

braak_pos <- braakGenes$entrez_id[braakGenes$r>0]
braak_neg <- braakGenes$entrez_id[braakGenes$r<0]

# Plot theme
theme <- theme(legend.position = "none",
               panel.background = element_blank(),
               axis.line = element_line(colour = "black"),
               axis.title =  element_text(size = 12),
               plot.title = element_text(size = 12, face = "bold"))

# axis limits for all plots
ctab <- Reduce(rbind, summaryDiffExpr)
xmax <- max(ctab$meanDiff)
xmin <- min(ctab$meanDiff)
ymax <- ceiling(max(ctab$'logp'[is.finite(ctab$'logp')]))

plotll <- lapply(names(summaryDiffExpr), function(rp){
  tab <- summaryDiffExpr[[rp]]
  # tab$info <- as.numeric(rownames(tab) %in% braak_pos)
  tab$info <- sapply(rownames(tab), function(x){
    if (x %in% braak_pos) 1
    else if (x %in% braak_neg) 2
    else 0
  })
  
  # tab$labels <- entrezId2Name(rownames(tab))
  # pd <- c("DNAJC13","SNCA","GCH1","INPP5F", "ASH1L", "ITPKB", "ELOVL7", "ZNF184", "SCARB2", "DDRGK1")
  # tab$labels[!tab$labels %in% pd] <- ""
  # tab[name2EntrezId(pd), "info"] <- 2
  
  # Plotting order of data points 
  tab$info <- as.factor(tab$info)
  order <- order(tab$info)
  tab <- tab[order, ]
  
  p <- ggplot(tab, aes(meanDiff, logp, colour = info)) +
    geom_point(size = 0.25, alpha = 0.3) +
    scale_colour_manual(values = c("0"="grey", "1"="red", "2"="blue")) +
    # geom_text_repel(label = tab$labels, fontface = "italic", colour = "black", size = 3, nudge_x = 0.2) +
    labs(x = "Fold-change", y = "-log10 P-value") +
    scale_x_continuous(limits = c(xmin, xmax), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, ymax), expand = c(0,0)) +
    ggtitle(paste("Braak stage ", gsub("braak", "", gsub("-", " vs ", rp)))) +
    theme
  
  name <- paste0("DiffExpr_braak/volcanoplot_", rp, ".pdf")
  pdf(file = name, 4, 3)
  print(p)
  dev.off()
})