# Differential expression between regions in UKBEC
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(biomaRt)
library(WGCNA)
library(plyr)
library(reshape2)
library(ggplot2)
load("../UKBEC/expr.maps.rda",verbose=T)
load("resources/braakGenes.RData")
load("resources/braakInfo.RData") # Braak colors
source("PD/t.test.table.R")

##############################################################################################

# Map affy ID to entrez IDs
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version = 92)
affyID <- expr.map$exprID
system.time({
  ukbecGeneID <- getBM(c('ensembl_gene_id', 'entrezgene', 'hgnc_symbol', "affy_huex_1_0_st_v2"), 
                       filters=c("affy_huex_1_0_st_v2"), mart=ensembl, values=affyID)
})
save(ukbecGeneID, file = "resources/ukbecGeneID.RData")
write.table(ukbecGeneID, file = "ukbecGeneID.txt", row.names = FALSE, quote = FALSE, sep = "\t")
load("resources/ukbecGeneID.RData")

##############################################################################################
# Probes to genes

ukbecGeneID$entrezgene <- as.character(ukbecGeneID$entrezgene)
ukbecGeneID$affy_huex_1_0_st_v2 <- as.character(ukbecGeneID$affy_huex_1_0_st_v2)

# Read expression data for all brain regions
regions <- c("CRBL", "FCTX", "HIPP","MEDU", "OCTX", "PUTM", "SNIG", "TCTX", "THAL", "WHMT")
regionExpr <- sapply(regions, function(r){
  fName <- paste0("../UKBEC/expr_", r, ".txt")
  x=read.csv(fName, header = TRUE, sep = " ", row.names = 1)
}, simplify = FALSE)
exprConcat <- Reduce(cbind, regionExpr) #Concatenate brain regions

# Mapping probes to genes
probes <- unique(ukbecGeneID$affy_huex_1_0_st_v2)
genes <- ukbecGeneID$entrezgene[match(probes, ukbecGeneID$affy_huex_1_0_st_v2)] # entrez IDs
nas <- is.na(genes)
probes <- probes[!nas]
genes <- genes[!nas]
exprConcat <- exprConcat[probes,]
probeSelection <- collapseRows(exprConcat, rowGroup = genes, rowID = probes, 
                                method = "maxRowVariance", connectivityBasedCollapsing = TRUE)

# Select probes for each data matrix
probes <- probes[probeSelection$selectedRow]
genes <- genes[probeSelection$selectedRow]
regionExpr <- lapply(regionExpr, function(x) {
  x <- x[probes, ]
  rownames(x) <- genes
  x
})
save(regionExpr, file = "../UKBEC/regionExpr.RData")

##############################################################################################
# T-test in UKBEC

roi <- c('1' = "MEDU", '3' = "SNIG", '5' = "TCTX", '6' = "FCTX")
regionpairs <- combn(roi, 2)
colnames(regionpairs) <- apply(regionpairs, 2, function(x) paste0(x[1], "-",  x[2]))
ttest <- alply(regionpairs, 2, function(x){
  df1 <- regionExpr[[x[1]]]
  df2 <- regionExpr[[x[2]]]
  t.test.table(df1,df2)
}, .dims = TRUE)
ttest <- simplify2array(ttest) # 3D array: genes x measures x region pairs

# Number of diff. genes
apply(ttest, c(3), function(x){
  sum(x[, "BH"] < 0.05 & abs(x[, "meanDiff"]) > 1)
})
save(ttest, file = "resources/ttest_ukbec.RData")

# Number of diff. genes
sum(abs(ttest_ukbec$fc) > 1 & ttest_ukbec$BH < 0.05)

##############################################################################################
# Plotting functions

prepare.data <- function(g){ # prepare ggplot dataframe for single genes
  expr <- sapply(roi, function(r)  unlist(regionExpr[[r]][g, ] ))
  df <- melt(expr)
  colnames(df) <- c("sample", "region", "expr")
  df$region <- factor(df$region, levels = unique(df$region))
  df
}

names(braakColors) <- gsub("braak", "", names(braakColors))
theme <- theme(panel.background = element_blank(), panel.grid = element_blank(), 
               axis.line = element_line(colour = "black"))

box.plot <- function(df, title){
  ggplot(df) + 
    geom_boxplot(aes(y = expr, x = region, fill = region)) +
    labs(x = "Brain region", y = "Expression (log2-transformed)") +
    ggtitle(title) +
    scale_x_discrete(expand=c(0.2,0)) +
    scale_fill_manual(values = braakColors, guide = FALSE) +
    theme
}

plot.pdf <- function(name, genes){ # For plots of single genes
  pdf(name, 2, 3)
  lapply(genes, function(g){
    title <- paste0(entrezId2Name(g))
    df <- prepare.data(g)
    p <- box.plot(df, title)
    print(p)
  })
  dev.off()
}

##############################################################################################
# Box plots

# Expression of Braak genes
braak_neg <- braakGenes$entrez_id[braakGenes$r < 0]
braak_pos <- braakGenes$entrez_id[braakGenes$r > 0]
braak_neg <- intersect(genes, braak_neg)
braak_pos <- intersect(genes, braak_pos)
braak <- list('r<0' = braak_neg, 'r>0' = braak_pos)

# Mean expression across Braak genes within regions
meanExpr <- lapply(roi, function(r){
  t <- sapply(braak, function(g){
    expr <- regionExpr[[r]][g, ]
    apply(expr, 2, mean)
  })
  melt(t)
})
meanExpr <- melt(meanExpr)
colnames(meanExpr) <- c("sample", "r", "variable", "expr", "region")
meanExpr$region <- factor(meanExpr$region, levels = unique(meanExpr$region))
p1 <- box.plot(meanExpr, "Expression of Braak genes in UKBEC") +
  facet_grid(.~r, scales = 'free', space = 'free', switch = "y")
pdf("boxplot_UKBEC.pdf", 4, 3)
p1
dev.off()

#boxplot of PD genes
plot.pdf("boxplot_UKBEC_PD_variant_genes.pdf", 
         name2EntrezId(c("DNAJC13", "SNCA", "GCH1", "INPP5F", "ASH1L", "ZNF184", "DDRGK1", "ITPKB", "ELOVL7", "SCARB2")))
