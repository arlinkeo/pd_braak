  # Differential expression MEDU (1) vs. FCTX (6)
# Map probe to genes UKBEC datasets
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(WGCNA)
library(reshape)
library(ggplot2)
load("../UKBEC/expr.maps.rda",verbose=T)
load("resources/ukbecGeneID.RData")
load("resources/braakGenes.RData")
load("resources/braakInfo.RData") # Braak colors

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
# regionExpr <- lapply(regionExpr, function(x) t(scale(t(x), center = TRUE)) ) # Normallize data
save(regionExpr, file = "../UKBEC/regionExpr.RData")

##############################################################################################

# T-test in UKBEC
ttest_ukbec <- as.data.frame(t(sapply(genes, function(g){
  r1 <- regionExpr$MEDU[g, ]
  r2 <- regionExpr$FCTX[g, ]
  t <- t.test(r1, r2)
  fc <- t$estimate[1]-t$estimate[2]
  pvalue <- t$p.value
  c(fc=unname(fc), pvalue=unname(pvalue))
})))
ttest_ukbec$BH <- p.adjust(ttest_ukbec$pvalue, method = "BH")
save(ttest_ukbec, file = "resources/ttest_ukbec.RData")

# Number of diff. genes
sum(abs(ttest_ukbec$fc) > 1 & ttest_ukbec$BH < 0.05)

# Braak genes
braak_neg <- braakGenes$entrez_id[braakGenes$braak_r < 0]
braak_pos <- braakGenes$entrez_id[braakGenes$braak_r > 0]
braak_neg <- intersect(genes, braak_neg)
braak_pos <- intersect(genes, braak_pos)
braak <- list('r<0' = braak_neg, 'r>0' = braak_pos)

# Mean expression across Braak genes within regions
roi <- c('1' = "MEDU", '3' = "SNIG", '5' = "TCTX", '6' = "FCTX")
meanExpr <- lapply(roi, function(r){
  t <- sapply(braak, function(g){
    expr <- regionExpr[[r]][g, ]
    apply(expr, 2, mean)
  })
  melt(t)
})
meanExpr <- melt.list(meanExpr)
colnames(meanExpr) <- c("sample", "r", "variable", "expr", "region")
meanExpr$region <- factor(meanExpr$region, levels = unique(meanExpr$region))

theme <- theme(panel.background = element_blank(), panel.grid = element_blank(), 
               axis.line = element_line(colour = "black"))

names(braakColors) <- gsub("braak", "", names(braakColors))

box.plot <- function(df, title){
  ggplot(df) + 
    geom_boxplot(aes(y = expr, x = region, fill = region)) +
    labs(x = "Braak stage", y = "Expression (log2-transformed)") +
    ggtitle(title) +
    scale_x_discrete(expand=c(0.2,0)) +
    scale_fill_manual(values = braakColors, guide = FALSE) +
    theme
}
p1 <- box.plot(meanExpr, "Expression of Braak genes in UKBEC") +
  facet_grid(.~r, scales = 'free', space = 'free', switch = "y")
p1

#boxplot of SNCA
expr_snca <- sapply(roi, function(r)  unlist(regionExpr[[r]]["6622", ] ))
expr_snca <- melt(expr_snca)
colnames(expr_snca) <- c("sample", "region", "expr")
expr_snca$region <- factor(expr_snca$region, levels = unique(expr_snca$region))

p2 <- box.plot(expr_snca, "Expression of SNCA in UKBEC")
p2
pdf("boxplot_UKBEC.pdf", 4, 3)
print(p1)
print(p2)
dev.off()