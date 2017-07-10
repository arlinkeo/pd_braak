# Plot braak stage expresssion
options(stringsAsFactors = FALSE)
setwd("C:/Users/dkeo/surfdrive/Parkinson")
library(ggplot2)
library(gridExtra)
library(grid)
library(reshape2)
library(plyr)

load("../ABA_Rdata/BrainExpr.RData")
donorNames <- names(brainExpr)
names(donorNames) <- donorNames
load("../ABA_Rdata/BrainExprNorm.RData")
load("resources/geneLabelCor.RData")
load("resources/diffGenesBraak.RData")
load("resources/braakLabels.RData")

# Mapping entrez IDs to gene symbols and vice versa
probeInfo <- read.csv("../ABA_human_processed/probe_info_2014-11-11.csv")
entrezId2Name <- function (x) {probeInfo$gene_symbol[match(x, probeInfo$entrez_id)]} #Input is vector
name2entrezId <- function (x) {probeInfo$entrez_id[match(x, probeInfo$gene_symbol)]} #Input is vector

#Average correlation across brains
avgCor <- apply(geneLabelCor, 1, mean)

#braak stage plot order 
braakNames <- sapply(names(diffGenesBraak1), function(n){tail(unlist(strsplit(n, split = "braak")), 1)})
names(diffGenesBraak1) <- braakNames
names(braakNames) <- braakNames
braakOrder <- c("0", braakNames)
names(braakOrder) <- NULL

###### Functions
# Get diff. expr. p-value for a gene and donor
pval.up <- function(gene, donor){
  sapply(diffGenesBraak1, function(braak){
    format(braak[[donor]][gene, "corr_pvalUp"], digits = 2)
  })
}

pval.down <- function(gene, donor){
  sapply(diffGenesBraak1, function(braak){
    format(braak[[donor]][gene, "corr_pvalDown"], digits = 2)
  })
}

# is.present <- function(x){
#   res <- sapply(diffGenesBraak2, function(bs){
#     upGenes <- bs[["upregulated"]]
#     downGenes <- bs[["downregulated"]]
#     upVec <- x %in% upGenes[ , "gene_symbol"]
#     downVec <- x %in% downGenes[ , "gene_symbol"]
#     upAndDown <- intersect(which(upVec), which(downVec)) # which gene (index) is up and downregulated? (should not be possible)
#     if(length(upAndDown) > 0) paste("Genes both up- and downregulated: ", paste(x[upAndDown], collapse = ", "), sep = "")
#     upVec <- as.numeric(upVec)
#     downVec <- -as.numeric(downVec)
#     apply(rbind(upVec, downVec), 2, sum)
#   })
#   res <- as.data.frame((res))
#   res <- if(length(x) == 1) t(res) else res
#   rownames(res) <- x
#   colnames(res) <- names(diffGenesBraak2)
#   res
# }

######### Box plot functions

# Number of samples
nSamples <- function(x){
  c(y = 9, label = length(x))
}

# single plot for braak stages
boxplot.braakstage <- function(tab, title){
  tab$braakstage <- factor(tab$braakstage, levels = braakOrder)
  p <- ggplot(tab, aes(x = braakstage, y = expr)) + geom_boxplot(aes(fill = color)) +
    stat_summary(fun.data = nSamples, geom = "text", size = 2.5) +
    labs(x = "Braak stage", y = "Expression") +
    guides(fill = FALSE) + 
    ggtitle(title)
  p
}

boxplot.with.pval <- function(tab, title, pvalUp.df, pvalDown.df){
  p <- boxplot.braakstage(tab, title)
  p + geom_text(data = pvalUp.df, aes(y = 0, label = pvalue), size = 2.5) +
    geom_text(data = pvalDown.df, aes(y = -1, label = pvalue), size = 2.5)
}

# multiple plots given multiple tables
multi.boxplot <- function(p, main){
  main = textGrob(main, gp=gpar(fontface="bold"))
  grid.arrange(grobs = p, top = main) 
}

# Create data frame
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
  geneExpr <- unlist(brainExpr[[d]][gene, ])
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

####### Heatmap functions


#######

# Box plot min and max genes
maxGene <- names(which(avgCor == max(avgCor)))
minGene <- names(which(avgCor == min(avgCor)))
exprMaxGene <- lapply(donorNames, function(d){ plot.gene(maxGene,d) })
exprMinGene <- lapply(donorNames, function(d){ plot.gene(minGene,d) })

#Expression of lysosome genes across Braak region
lysosomeGenes <- unlist(read.table("lysosome_geneset.txt", header = FALSE, comment.char = "#", sep = "\n"))
lysosomeGenes <- as.character(name2entrezId(lysosomeGenes))
lysosomeGenes<- lysosomeGenes[!is.na(lysosomeGenes)]
avgCorLyso <- avgCor[lysosomeGenes]
names(avgCorLyso) <- entrezId2Name(names(avgCorLyso))

# average across genes
exprLyso1 <- lapply(donorNames, function(d){
  geneExpr <- brainExpr[[d]][lysosomeGenes, ]
  avgSampleExpr <- apply(geneExpr, 2, mean)
  tab <- expr.df(avgSampleExpr, d)
  boxplot.braakstage(tab, d)
})
#Average within region
exprLyso2 <- lapply(donorNames, function(d){
  geneExpr <- brainExpr[[d]][lysosomeGenes, ]
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
sortLyso <- sort(avgCorLyso)
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

################3 Heat maps of expression of lysosome genes
sortedLysoIds <- as.character(name2entrezId(names(sortLyso)))

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

heatmap.expr <- function(tab, main){
  tab.m <- melt(tab, id.vars = c("sample", "braakstage"))
  tab.m$braakstage <- factor(tab.m$braakstage, levels = c(1:6)) # keep order
  ggplot(tab.m, aes(sample, variable, group = braakstage)) +
    geom_tile(aes(fill = value)) +
    facet_grid(~braakstage, scales = "free", space = "free") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", name = "Expression") +
    scale_x_discrete(position = "top") +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(face = "italic", size = 6),
          axis.ticks = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black", size = 0.6)#,
          # plot.margin = unit(c(8,0,1,0), unit = "cm")
    ) +
    labs(x = "Braak stage", y = "Lysosome gene") +
    ggtitle(main)
}

tabLyso <- lapply(donorNames, function(d){
  tab <- genes.expr.braak(sortedLysoIds, d)  
  heatmap.expr(tab, d)
})


#Heatmap of correlations across donors of lysosome genes
corLyso <- geneLabelCor[sortedLysoIds, ]
rownames(corLyso) <- entrezId2Name(rownames(corLyso))
corLyso.m <- melt(corLyso)
corLyso.m$Var2 <- factor(corLyso.m$Var2, levels = donorNames) # keep order of cols

corLyso.heatmap <- ggplot(corLyso.m, aes(Var2, Var1)) +
  geom_tile(aes(fill = value), colour = "black") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", name = "Correlation") +
  scale_x_discrete(position = "top") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 12),
        axis.text.y = element_text(face = "italic")) +
  labs(x = "Braak stage", y = "Lysosome gene") +
  coord_fixed(0.4)

############## Plot everything
pdf(file = "expr_corr_boxplot.pdf",8,8)

hist(avgCor, main = "Correlation between gene expression and braak regions")
multi.boxplot(exprMaxGene, main = paste("Gene with max. positive corrrelation: ", entrezId2Name(maxGene), sep = ""))
multi.boxplot(exprMinGene, main = paste("Gene with min. negative corrrelation: ", entrezId2Name(minGene), sep = ""))
hist(corLyso)
multi.boxplot(exprLyso1, main = "Expression of lysosome genes in Braak regions \n (averaged across genes)")
multi.boxplot(exprLyso2, main = "Expression of lysosome genes in Braak regions \n (averaged across samples)")
multi.boxplot(exprLyso3a, main = "Expression of top 3 correlated lysosome genes in Braak regions")
multi.boxplot(exprLyso3b, main = "Expression of low correlated lysosome genes in Braak regions")
multi.boxplot(exprHiGenes1, main = "Expression of high impact genes in Braak regions")
multi.boxplot(exprHiGenes2, main = "Expression of high impact genes in Braak regions")

dev.off()


pdf(file = "expr_corr_heatmap.pdf",8,9)
multi.boxplot(tabLyso, main = "Expression of lysosome genes in Braak stages")

exprLyso.heatmap
corLyso.heatmap

dev.off()
####################

