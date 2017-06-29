#Correlate expression  with braakstage labels

options(stringsAsFactors = FALSE)
setwd("C:/Users/dkeo/surfdrive/Parkinson")
library(ggplot2)
library(gridExtra)
library(grid)
library(reshape2)
library(plyr)

load("resources/braakStages.RData")
braakNames <- names(braakStages)
names(braakNames) <- braakNames
load("../polyQ_coexpression/resources/BrainExpr.RData")
donorNames <- names(brainExpr)
names(donorNames) <- donorNames
sampleIds <- lapply(brainExpr, colnames)
probeInfo <- read.csv("../ABA_human_processed/probe_info_2014-11-11.csv")
entrezId2Name <- function (x) {probeInfo$gene_symbol[match(x, probeInfo$entrez_id)]} #Input is vector
name2entrezId <- function (x) {probeInfo$entrez_id[match(x, probeInfo$gene_symbol)]} #Input is vector

### Functions
#Revert list in list
revert.list <- function(ll){
  lapply(donorNames, function(d){
    t(sapply(ll, function(bs){
      bs[[d]]
    }))
  })
}
#generate label vector for all samples
label.vector <- function(ll){
  lapply(donorNames, function(dn){
    d <- ll[[dn]]
    labels <- apply(d, 2, function(v){
      s <- which(v == 1)
      ifelse(length(s) == 0, 0, tail(unlist(strsplit(names(s), split = "braak")), 1))
    })
    names(labels) <- sampleIds[[dn]]
    labels
  })
}

### Braak labels 1 to 6 for each sample ###
bsPerDonor <- revert.list(braakStages[1:6])
braakLabels <- label.vector(bsPerDonor)

# #Count non-zero labels, if >1 there are multiple labels assigned to a sample
# countLabels <- sapply(bsPerDonor, function(d){apply(d, 2, function(v){sum(v!=0)})})
# lapply(countLabels, function(d){which(d > 2)})

# #Sorted labels
# sortedLabels <- lapply(braakLabels, sort)
# sortedSampleIds <- lapply(sortedLabels, names)

#Correlate gene expression to Braak labels for each brain (without non-Braak region)
geneLabelCor <- sapply(donorNames, function(d){
  labels <- as.numeric(braakLabels[[d]])
  nbRows <- as.logical(nonBraak[[d]])
  labels <- labels[!nbRows]
  expr <- brainExpr[[d]][, !nbRows]
  apply(expr, 1, function(g){# corr. each gene
    g <- unlist(g)
    cor(g, labels) # Pearson's r by default
  })
})

### Merged Braak labels 1-3 and 4-6 for each sample ###
mergedBsPerDonor <- revert.list(braakStages[7:8])
mergedBraakLabels <- label.vector(mergedBsPerDonor)

#braak stage plot order 
braakOrder <- c("0", sapply(braakNames, function(n){tail(unlist(strsplit(n, split = "braak")), 1)}))
names(braakOrder) <- NULL

#### Bar plot functions
# single plot for braak stages
boxplot.braakstage <- function(tab, title){
  tab$braakstage <- factor(tab$braakstage, levels = braakOrder)
  p <- ggplot(tab, aes(x = braakstage, y = expr)) + geom_boxplot(aes(fill = color)) +
    labs(x = "Braak stage", y = "Expression") +
    guides(fill = FALSE) + 
    ggtitle(title)
  p
}

# multiple plots given multiple tables
multi.boxplot <- function(tabList, main){
  tNames <- names(tabList)
  p <- lapply(tNames, function(tn){
    boxplot.braakstage(tabList[[tn]], tn)
  })
  main = textGrob(main, gp=gpar(fontface="bold"))
  grid.arrange(grobs = p, top = main) 
}

# get gene(s) expression and braak labels for each donor
gene.expr.lab <- function(gene){
  lapply(donorNames, function(d){
    geneExpr <- unlist(brainExpr[[d]][gene, ])
    labels <- braakLabels[[d]]
    mergedLabels <- mergedBraakLabels[[d]]
    tab1 <- data.frame(expr = geneExpr, braakstage = labels, color = mergedLabels)
    tab2 <- data.frame(expr = geneExpr, braakstage = mergedLabels, color = mergedLabels)
    tab <- rbind(tab1, tab2)
  })
}
####

pdf(file = "expr_corr.pdf")

#Average correlation across brains
avgCor <- apply(geneLabelCor, 1, mean)
hist(avgCor, main = "Correlation between gene expression and braak regions")

# Bar plot min and max genes
maxGene <- names(which(avgCor == max(avgCor)))
minGene <- names(which(avgCor == min(avgCor)))
exprMaxGene <- gene.expr.lab(maxGene)
multi.boxplot(exprMaxGene, main = "Gene with max. positive corrrelation")
exprMinGene <- gene.expr.lab(minGene)
multi.boxplot(exprMinGene, main = "Gene with min. negative corrrelation")

#Expression of lysosome genes across Braak region
lysosomeGenes <- unlist(read.table("lysosome_geneset.txt", header = FALSE, comment.char = "#", sep = "\n"))
lysosomeGenes <- as.character(name2entrezId(lysosomeGenes))
lysosomeGenes<- lysosomeGenes[!is.na(lysosomeGenes)]
corLyso <- avgCor[lysosomeGenes]
names(corLyso) <- entrezId2Name(names(corLyso))
# average across genes
exprLyso1 <- lapply(donorNames, function(d){
  geneExpr <- brainExpr[[d]][lysosomeGenes, ]
  avgSampleExpr <- apply(geneExpr, 2, mean)
  labels <- braakLabels[[d]]
  mergedLabels <- mergedBraakLabels[[d]]
  tab1 <- data.frame(expr = avgSampleExpr, braakstage = labels, color = mergedLabels)
  tab2 <- data.frame(expr = avgSampleExpr, braakstage = mergedLabels, color = mergedLabels)
  tab <- rbind(tab1, tab2)
})
multi.boxplot(exprLyso1, main = "Expression of lysosome genes in Braak regions \n (averaged across genes)")
#Average within region
exprLyso2 <- lapply(donorNames, function(d){
  geneExpr <- brainExpr[[d]][lysosomeGenes, ]
  labels <- braakLabels[[d]]
  mergedLabels <- mergedBraakLabels[[d]]
  # uniqueLabels <- unique(c(labels, mergedLabels))
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
})
multi.boxplot(exprLyso2, main = "Expression of lysosome genes in Braak regions \n (averaged across samples)")

dev.off()



