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

#Correlate gene expression to Braak labels for each brain (without non-Braak region)
# geneLabelCor <- sapply(donorNames, function(d){
#   labels <- as.numeric(braakLabels[[d]])
#   nbRows <- as.logical(nonBraak[[d]])
#   labels <- labels[!nbRows]
#   expr <- brainExpr[[d]][, !nbRows]
#   apply(expr, 1, function(g){# corr. each gene
#     g <- unlist(g)
#     cor(g, labels) # Pearson's r by default
#   })
# })
# save(geneLabelCor, file = "resources/geneLabelCor.RData")
load("resources/geneLabelCor.RData")

### Merged Braak labels 1-3 and 4-6 for each sample ###
mergedBsPerDonor <- revert.list(braakStages[7:8])
mergedBraakLabels <- label.vector(mergedBsPerDonor)

#braak stage plot order 
braakOrder <- c("0", sapply(braakNames, function(n){tail(unlist(strsplit(n, split = "braak")), 1)}))
names(braakOrder) <- NULL

######### load Differential expression info
load("resources/diffGenesBraak.RData")
names(diffGenesBraak1) <- braakOrder[-1]

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

# Presence of high impact genes
hiGenes <- c("GBA", "LRRK2", "PINK1", "PARK7", "SCNA", "VPS35", "DNAJC13", "CHCHD2") # High impact genes
# presence_hiGenes <- is.present(hiGenes)
#########

#### Box plot functions
nSamples <- function(x){
  c(y = 9, label = length(x))
}

# single plot for braak stages
boxplot.braakstage <- function(tab, title){
  tab$braakstage <- factor(tab$braakstage, levels = braakOrder)
  p <- ggplot(tab, aes(x = braakstage, y = expr)) + geom_boxplot(aes(fill = color)) +
    stat_summary(fun.data = nSamples, geom = "text", size = 2) +
    labs(x = "Braak stage", y = "Expression") +
    guides(fill = FALSE) + 
    ggtitle(title)
  p
}

boxplot.with.pval <- function(tab, title, pvalUp.df, pvalDown.df){
  p <- boxplot.braakstage(tab, title)
  p + geom_text(data = pvalUp.df, aes(y = 0, label = pvalue), size = 2) +
    geom_text(data = pvalDown.df, aes(y = -1, label = pvalue), size = 2)
}

# multiple plots given multiple tables
multi.boxplot <- function(p, main){
  main = textGrob(main, gp=gpar(fontface="bold"))
  grid.arrange(grobs = p, top = main) 
}

# get gene(s) expression, braak labels and p-value (diff. expression) for each donor
gene.expr.lab <- function(gene, d){
    geneExpr <- unlist(brainExpr[[d]][gene, ])
    labels <- braakLabels[[d]]
    mergedLabels <- mergedBraakLabels[[d]]
    tab1 <- data.frame(expr = geneExpr, braakstage = labels, color = mergedLabels)
    tab2 <- data.frame(expr = geneExpr, braakstage = mergedLabels, color = mergedLabels)
    rows <- which(tab2$braakstage != 0)#otherwise expr data of Non-braak is double
    tab2 <- tab2[rows, ]
    tab <- rbind(tab1, tab2)
}

####

plot.gene <- function(gene, d){
  tab <- gene.expr.lab(gene,d)
  pvaluesUp <- c('0' = "", pval.up(gene, d))
  pvalUp.df <- data.frame(braakstage = names(pvaluesUp), pvalue = pvaluesUp)
  pvaluesDown <- c('0' = "", pval.down(gene, d))
  pvalDown.df <- data.frame(braakstage = names(pvaluesDown), pvalue = pvaluesDown)
  boxplot.with.pval(tab, d, pvalUp.df, pvalDown.df)
}

pdf(file = "expr_corr.pdf")

#Average correlation across brains
avgCor <- apply(geneLabelCor, 1, mean)
hist(avgCor, main = "Correlation between gene expression and braak regions")

# Bar plot min and max genes
maxGene <- names(which(avgCor == max(avgCor)))
minGene <- names(which(avgCor == min(avgCor)))
exprMaxGene <- lapply(donorNames, function(d){ plot.gene(maxGene,d) })
multi.boxplot(exprMaxGene, main = "Gene with max. positive corrrelation")
exprMinGene <- lapply(donorNames, function(d){ plot.gene(minGene,d) })
multi.boxplot(exprMinGene, main = "Gene with min. negative corrrelation")

#Expression of lysosome genes across Braak region
lysosomeGenes <- unlist(read.table("lysosome_geneset.txt", header = FALSE, comment.char = "#", sep = "\n"))
lysosomeGenes <- as.character(name2entrezId(lysosomeGenes))
lysosomeGenes<- lysosomeGenes[!is.na(lysosomeGenes)]
corLyso <- avgCor[lysosomeGenes]
names(corLyso) <- entrezId2Name(names(corLyso))
hist(corLyso)
# average across genes
exprLyso1 <- lapply(donorNames, function(d){
  geneExpr <- brainExpr[[d]][lysosomeGenes, ]
  avgSampleExpr <- apply(geneExpr, 2, mean)
  labels <- braakLabels[[d]]
  mergedLabels <- mergedBraakLabels[[d]]
  tab1 <- data.frame(expr = avgSampleExpr, braakstage = labels, color = mergedLabels)
  tab2 <- data.frame(expr = avgSampleExpr, braakstage = mergedLabels, color = mergedLabels)
  rows <- which(tab2$braakstage != 0)#otherwise expr data of Non-braak is double
  tab2 <- tab2[rows, ]
  tab <- rbind(tab1, tab2)
  boxplot.braakstage(tab, d)
})
multi.boxplot(exprLyso1, main = "Expression of lysosome genes in Braak regions \n (averaged across genes)")
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
multi.boxplot(exprLyso2, main = "Expression of lysosome genes in Braak regions \n (averaged across samples)")

dev.off()