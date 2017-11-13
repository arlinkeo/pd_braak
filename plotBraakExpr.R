# Plot braak stage expresssion
options(stringsAsFactors = FALSE)
setwd("C:/Users/dkeo/surfdrive/Parkinson")
library(ggplot2)
library(gridExtra)
library(grid)
library(reshape2)
library(plyr)
library(ggrepel)

load("../ABA_Rdata/BrainExprNorm.RData")
donorNames <- names(brainExprNorm)
names(donorNames) <- donorNames
load("resources/geneLabelCor.RData")
load("resources/diffGenesBraak.RData") # for looking up p-value (2 diff. data structures)
load("resources/braakLabels.RData") # Braak stage label vectors
load("resources/summaryCorr.RData")

susceptibilityGenes <- as.character(name2entrezId(c("INPP5F", "TMEM175", "ASH1L", "MAPT", "RIT1", "C14orf83", "STK39", "GPNMB", "BST1", 
              "SIPA1L2", "DLG2", "NUCKS1", "GCH1", "MCCC1", "FAM47E", "BCKDK", "TMPRSS9", "UBOX5", 
              "CCDC62", "SYNJ1", "EIF4G1", "FBXO7", "C20orf30", "POLG", "VSP13C", "PLA2G6")))
susceptibilityGenes <- susceptibilityGenes[!is.na(susceptibilityGenes)]
hlaGenes <- as.character(name2entrezId(c("HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DQB1")))

# Mapping entrez IDs to gene symbols and vice versa
probeInfo <- read.csv("../ABA_human_processed/probe_info_2014-11-11.csv")
entrezId2Name <- function (x) {probeInfo$gene_symbol[match(x, probeInfo$entrez_id)]} #Input is vector
name2entrezId <- function (x) {probeInfo$entrez_id[match(x, probeInfo$gene_symbol)]} #Input is vector

###### p-value/fold-change Functions
# Get diff. expr. p-value for a gene and donor
pval.up <- function(gene, donor){
  sapply(diffGenesList, function(braak){
    format(braak[[donor]][gene, "pval_up"], digits = 2)
  })
}

pval.down <- function(gene, donor){
  sapply(diffGenesList, function(braak){
    format(braak[[donor]][gene, "pval_down"], digits = 2)
  })
}

# pval.2tail <- function(gene, donor){
#   sapply(diffGenesList, function(braak){
#     format(braak[[donor]][gene, "pval_2tail"], digits = 2)
#   })
# }
# 
# fold.change <- function(gene, donor){
#   sapply(diffGenesList, function(braak){
#     format(braak[[donor]][gene, "fold-change"], digits = 2)
#   })
# }

######### Box plot functions

#braak stage plot order 
braakNames <- sapply(names(diffGenesList), function(n){tail(unlist(strsplit(n, split = "braak")), 1)})
names(braakNames) <- braakNames
braakOrder <- c("0", braakNames)
names(braakOrder) <- NULL

# single boxplot for braak stages
boxplot.braakstage <- function(tab, title){
  tab$braakstage <- factor(tab$braakstage, levels = braakOrder)
    nSamples <- function(x){c(y = max(tab[["expr"]]) + 0.5, label = length(x))}# Number of samples
  p <- ggplot(tab, aes(x = braakstage, y = expr)) + geom_boxplot(aes(fill = color)) +
    stat_summary(fun.data = nSamples, geom = "text", size = 2.5) +
    labs(x = "Braak stage", y = "Expression") +
    guides(fill = FALSE) + 
    ggtitle(title)
  p
}

boxplot.with.pval <- function(tab, title, pvalUp.df, pvalDown.df){
  p <- boxplot.braakstage(tab, title)
  p + geom_text(data = pvalUp.df, aes(y = min(tab[["expr"]]) - .5, label = pvalue), size = 2.5) +
    geom_text(data = pvalDown.df, aes(y = min(tab[["expr"]]) - 1, label = pvalue), size = 2.5)
}

# multiple plots on a grid
multi.boxplot <- function(p, main){
  main = textGrob(main, gp=gpar(fontface="bold"))
  grid.arrange(grobs = p, top = main) 
}

# Create data frame (boxplot function input) given expression matrix in a donor
expr.df <- function(expr, d){
  labels <- braakLabels[[d]]
  mergedLabels <- mergedBraakLabels[[d]]
  tab1 <- data.frame(expr = expr, braakstage = labels, color = mergedLabels)
  tab2 <- data.frame(expr = expr, braakstage = mergedLabels, color = mergedLabels)
  rows <- which(tab2$braakstage != 0)#otherwise expr data of Non-braak is double
  tab2 <- tab2[rows, ]
  tab <- rbind(tab1, tab2)
}
  
# get gene(s) expression, braak labels and p-value (diff. expression) for each donor
gene.expr.lab <- function(gene, d){
  geneExpr <- unlist(brainExprNorm[[d]][gene, ])
  expr.df(geneExpr, d)
}

#Boxplot with p-values given gene and donor
plot.gene <- function(gene, d){
  tab <- gene.expr.lab(gene,d)
  pvaluesUp <- c('0' = "", pval.up(gene, d))
  pvalUp.df <- data.frame(braakstage = names(pvaluesUp), pvalue = pvaluesUp)
  pvaluesDown <- c('0' = "", pval.down(gene, d))
  pvalDown.df <- data.frame(braakstage = names(pvaluesDown), pvalue = pvaluesDown)
  corr <- format(geneLabelCor[gene, d], digits = 2)
  main <- paste(d, ", r=", corr, sep = "")
  boxplot.with.pval(tab, main, pvalUp.df, pvalDown.df)
}
############

#Average correlation across brains
# avgCor <- apply(geneLabelCor, 1, mean)
sumCor <- sapply(summaryCorr, function(g)g["summary", "z"])
#Sort genes based on avgCor
sort.genes <- function(x){
  x <- sumCor[x]
  names(x) <- entrezId2Name(names(x))
  sort(x)
}


# Box plot min and max genes
maxGene <- names(which(avgCor == max(avgCor)))
exprMaxGene <- lapply(donorNames, function(d){ plot.gene(maxGene,d) })
minGene <- names(which(avgCor == min(avgCor)))
exprMinGene <- lapply(donorNames, function(d){ plot.gene(minGene,d) })

#Expression of lysosome genes across Braak region
lysosomeGenes <- unlist(read.table("lysosome_geneset.txt", header = FALSE, comment.char = "#", sep = "\n"))
lysosomeGenes <- as.character(name2entrezId(lysosomeGenes))
lysosomeGenes<- lysosomeGenes[!is.na(lysosomeGenes)]
avgCorLyso <- avgCor[lysosomeGenes]
names(avgCorLyso) <- entrezId2Name(names(avgCorLyso))

# average across genes
exprLyso1 <- lapply(donorNames, function(d){
  geneExpr <- brainExprNorm[[d]][lysosomeGenes, ]
  avgSampleExpr <- apply(geneExpr, 2, mean)
  tab <- expr.df(avgSampleExpr, d)
  boxplot.braakstage(tab, d)
})
#Average within region
exprLyso2 <- lapply(donorNames, function(d){
  geneExpr <- brainExprNorm[[d]][lysosomeGenes, ]
  labels <- braakLabels[[d]]
  mergedLabels <- mergedBraakLabels[[d]]
  avgGeneExpr <- sapply(braakOrder, function(s){
    cols <- unique(c(which(labels == s), which(mergedLabels == s)))
    stageExpr <- geneExpr[ , cols]
    apply(stageExpr, 1, mean)
  })
  tab <- melt(avgGeneExpr)
  colnames(tab) <- c("gene", "braakstage", "expr")
  color <- mapvalues(tab$braakstage, from = c(1:6), to = c(rep("1-3",3), rep("4-6",3)))
  tab$color <- color
  tab
  boxplot.braakstage(tab, d)
})

#Plot each lysosome in donor9861
sortLyso <- sort.genes(lysosomeGenes)
topLyso <- c(head(sortLyso,3), tail(sortLyso,3))
topLyso <- as.character(sapply(names(topLyso), name2entrezId))
exprLyso3a <- lapply(topLyso, function(gene){
  d <- donorNames[1]
  plot.gene(gene,d)
})
nonCorLyso <- as.character(sapply(names(head(sort(abs(avgCorLyso)),6)), name2entrezId))
exprLyso3b <- lapply(nonCorLyso, function(gene){
  d <- donorNames[1]
  plot.gene(gene,d)
})

# Presence of high impact genes
hiGenes <- c("GBA", "LRRK2", "PINK1", "PARK7", "SNCA", "VPS35", "DNAJC13", "CHCHD2") # High impact genes
# presence_hiGenes <- is.present(hiGenes)
hiGenes <- as.character(sapply(hiGenes, name2entrezId))
exprHiGenes1 <- lapply(hiGenes[1:4], function(gene){
  d <- donorNames[1]
  plot.gene(gene,d)
})
exprHiGenes2 <- lapply(hiGenes[5:8], function(gene){
  d <- donorNames[1]
  plot.gene(gene,d)
})

############## Plot boxplots and histograms
pdf(file = "corr_boxplot.pdf",8,8)

hist(avgCor, main = "Correlation between gene expression and braak regions")
multi.boxplot(exprMaxGene, main = paste("Gene with max. positive corrrelation: ", entrezId2Name(maxGene), sep = ""))
multi.boxplot(exprMinGene, main = paste("Gene with min. negative corrrelation: ", entrezId2Name(minGene), sep = ""))
hist(avgCorLyso)
multi.boxplot(exprLyso1, main = "Expression of lysosome genes in Braak regions \n (averaged across genes)")
multi.boxplot(exprLyso2, main = "Expression of lysosome genes in Braak regions \n (averaged across samples)")
multi.boxplot(exprLyso3a, main = "Expression of top 3 correlated lysosome genes in Braak regions")
multi.boxplot(exprLyso3b, main = "Expression of low correlated lysosome genes in Braak regions")
multi.boxplot(exprHiGenes1, main = "Expression of high impact genes in Braak regions")
multi.boxplot(exprHiGenes2, main = "Expression of high impact genes in Braak regions")

dev.off()

####### Heatmap functions

genes.expr.braak <- function(genes, d){
  expr <- brainExprNorm[[d]] # Normalized expression data
  labels <- braakLabels[[d]]
  samples <- labels != 0
  expr <- expr[genes, samples]
  labels <- as.integer(labels[samples])
  rownames(expr) <- entrezId2Name(rownames(expr))
  expr <- cbind(sample= colnames(expr), braakstage = labels, as.data.frame(t(expr)))
  expr
}

plot.heatmap <- function(tab, main, zlim = NULL){
  tab.m <- melt(tab, id.vars = c("sample", "braakstage"))
  tab.m$braakstage <- factor(tab.m$braakstage, levels = c(1:6)) # keep order
  ggplot(tab.m, aes(sample, variable, group = braakstage)) +
    geom_tile(aes(fill = value)) +
    facet_grid(~braakstage, scales = "free", space = "free") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", name = "Expression", limits = zlim) +
    scale_x_discrete(position = "top") +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(face = "italic", size = 2.5),
          axis.ticks = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black", size = 0.6)#,
          # plot.margin = unit(c(8,0,1,0), unit = "cm")
    ) +
    labs(x = "Braak stage", y = "Gene") +
    ggtitle(main)
}

exprLimit <- function(dfList){
  dfL <- lapply(dfList, function(df){
    df[ , !(colnames(df) == c("sample", "braakstage"))]
  })
  allValues <- Reduce(rbind, dfL)
  range(allValues)
}

expr.heatmaps <- function(genes) {
  tabs <- lapply(donorNames, function(d){
    genes.expr.braak(genes, d)
  })
  lim <- exprLimit(tabs)
  plots <- lapply(donorNames, function(d){
    t <- tabs[[d]]
    plot.heatmap(t, d, lim)
  })
}

cor.heatmap <- function(genes){
  cor <- geneLabelCor[genes, ]
  rownames(cor) <- entrezId2Name(rownames(cor))
  cor.m <- melt(cor)
  cor.m$Var2 <- factor(cor.m$Var2, levels = donorNames) # keep order of cols
  
  p <- ggplot(cor.m, aes(Var2, Var1)) +
    geom_tile(aes(fill = value), colour = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", name = "Correlation") +
    scale_x_discrete(position = "top") +
    theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 12),
          axis.text.y = element_text(face = "italic")#,
          #plot.margin = unit(c(10,0,10,0), unit = "cm")
          ) +
    labs(x = "Donor", y = "Gene") +
    coord_fixed(1)
}
######

# Heat maps of expression of lysosome genes
sortedLysoIds <- as.character(name2entrezId(names(sortLyso)))
tabLyso <- expr.heatmaps(sortedLysoIds)
#Heatmap of correlations across donors of lysosome genes
corLyso <- cor.heatmap(sortedLysoIds)

# Heat maps of expression of high impact genes
sortHiGenes <- sort.genes(hiGenes)
sortHiGenes <- as.character(name2entrezId(names(sortHiGenes)))
tabHiGenes <- expr.heatmaps(sortHiGenes)
#Heatmap of correlations across donors of high impact genes
corHiGenes <- cor.heatmap(sortHiGenes)

######## Plot heatmaps

pdf(file = "expr_corr_heatmap2.pdf",12,8)
multi.boxplot(tabLyso, main = "Expression of lysosome genes in Braak stages")
corLyso
multi.boxplot(tabHiGenes, main = "Expression of high impact genes in Braak stages")
corHiGenes

dev.off()
####################

### Volcano plots of diff. expr. genes
fullBraakNames <- names(diffGenesList)
braakOrder2 <- fullBraakNames[c(1:3, 7, 4:6, 8)]
braak1to6 <- fullBraakNames[1:6]
braakMerged <- fullBraakNames[7:8]

volcano.plot <- function(tab, bs, labels){
  ggplot(tab, aes(fold_change, pval_2tail, colour = info)) +
    geom_point(alpha = 1, size=1.5) +
    scale_colour_manual(values = c("0"="grey", "1"="red", "2"="black", "3" = "blue",
                                   "4" = "green", "5" = "pink")) +
    # geom_text(label = labels, colour = "black", size = 3, nudge_x = 0.2) +
    geom_text_repel(label = labels, colour = "black", size = 3, nudge_x = 0.2) +
    theme(legend.position = "none",
          panel.background = element_rect(fill = "white"),
          axis.line = element_line(colour = "black"),
          axis.title =  element_text(size = 12),
          plot.title = element_text(size = 12, face = "bold")
    ) +
    geom_vline(xintercept = 0, colour = "black") +
    geom_hline(yintercept = -log10(0.05), colour = "black") +
    labs(x = "log2 fold-change", y = "-log10 p-value") +
    ggtitle(paste("Braak stage", tail(unlist(strsplit(bs, split = "braak")), 1)))
}

v.plot.list <- function(braakList){
  lapply(braakList, function(bs){
    tab <- diffGenesList[[bs]][[d]]
    tab <- tab[, c("fold-change", "pval_2tail")]
    colnames(tab) <- c("fold_change", "pval_2tail")
    signifGenes <- rownames(tab)[tab$pval_2tail < 0.05 & abs(tab$fold_change) > 1]
    tab$info <- as.numeric(tab$pval_2tail < 0.05 & abs(tab$fold_change) > 1)
    tab[hiGenes, "info"] <- 2
    # tab[lysosomeGenes, "info"] <- 3
    # tab[susceptibilityGenes, "info"] <- 4
    # tab[hlaGenes, "info"] <- 5
    tab <- tab[order(tab$info),]
    tab$info <- as.factor(tab$info)
    tab$pval_2tail <- -log10(tab$pval_2tail)
    labels <- entrezId2Name(rownames(tab))
    labels[!(labels %in% entrezId2Name(c(hiGenes, signifGenes)))] <- ""
    volcano.plot(tab, bs, labels)
  })
}

d <- donorNames[1]

png(file = "diff_expr_volcanoplot1.png",900,600)
# for (d in donorNames) {
  vPlots1to6 <- v.plot.list(braak1to6)
  main = textGrob(d, gp=gpar(fontface="bold"))
  grid.arrange(grobs = vPlots1to6, top = main, nrow = 2, ncol = 3) 
# }
dev.off()

png(file = "diff_expr_volcanoplot2.png",600,300)
  d <- donorNames[1]
  main = textGrob(d, gp=gpar(fontface="bold"))
  vPlotsMerged <- v.plot.list(braakMerged)
  grid.arrange(grobs = vPlotsMerged, top = main, nrow = 1, ncol = 2) 
dev.off()
