# select significant genes based on significant summary estimate
setwd("C:/Users/dkeo/surfdrive/pd_braak")
library("metap")
library(gplots)
library(ggplot2)
library(reshape2)
library("RDAVIDWebService")
source("PD/base_script.R")
load("../ABA_Rdata/BrainExpr.RData")
load("resources/braakInfo.RData")
load("resources/not_celltype_corrected/summaryDiffExpr.RData")
load("resources/not_celltype_corrected/summaryLabelCor.RData")

########## Prepare data ##########

# Extract summary statistics, BH-correct
summaryDiffExpr <- lapply(summaryDiffExpr, function(rp){
  tab <- do.call(rbind.data.frame, lapply(rp, function(g) g["summary",]))
  tab$BH <- p.adjust(tab$pvalue, method = "BH")
  tab
})
# Select region pair
diffExpr <- summaryDiffExpr$`braak1-braak6`

labelCor <- do.call(rbind.data.frame, lapply(summaryLabelCor, function(g) g["summary",]))
labelCor$BH <- p.adjust(labelCor$pvalue, method = "BH")

########## Selecting Braak-related genes ##########

top10 <- floor(length(ahba.genes())*0.1) # number of genes in top 10%

# Get top 10% correlated genes
order <- rev(order(abs(labelCor$r))) # order absolute corr.

corrGenes <- rownames(labelCor)[order[1:top10]] # top 10% genes
corrGenes_neg <- corrGenes[labelCor[corrGenes, "r"] < 0]
corrGenes_pos <- corrGenes[labelCor[corrGenes, "r"] > 0]
max(abs(labelCor[corrGenes, "r"]))
min(abs(labelCor[corrGenes, "r"]))

# Top 10% mean Diff. expressed genes braak 1 vs. braak 6
order <- rev(order(abs(diffExpr$meanDiff))) # order absolute corr.
diffGenes1 <- rownames(diffExpr)[order[1:top10]] # top 10% genes
diffGenes1_pos <- diffGenes1[diffExpr[diffGenes1, "meanDiff"] > 0]
diffGenes1_neg <- diffGenes1[diffExpr[diffGenes1, "meanDiff"] < 0]
max(abs(diffExpr[diffGenes1, "meanDiff"]))
min(abs(diffExpr[diffGenes1, "meanDiff"]))

# Top 10% significant Diff. expressed genes braak 1 vs. braak 6
order <- order(diffExpr$BH) # order absolute corr.
diffGenes2 <- rownames(diffExpr)[order[1:top10]] # top 10% genes
diffGenes2_pos <- diffGenes2[diffExpr[diffGenes2, "meanDiff"] > 0]
diffGenes2_neg <- diffGenes2[diffExpr[diffGenes2, "meanDiff"] < 0]
max(diffExpr[diffGenes2, "BH"])
min(diffExpr[diffGenes2, "BH"])

# Venn diagram to visualize overlap of selections
criteria <- c("r", "fc", "pval_fc")
ll <- list(corrGenes, diffGenes1, diffGenes2)
ll1 <- list(corrGenes_neg, diffGenes1_neg, diffGenes2_neg)
ll2 <- list(corrGenes_pos, diffGenes1_pos, diffGenes2_pos)
names(ll) <- criteria
names(ll1) <- criteria
names(ll2) <- criteria
venn <- venn(ll)
venn1 <- venn(ll1)
venn2 <- venn(ll2)

# Selection based on 2 criteria and collect info of selected genes
braakGenes <- attr(venn, "intersections")[[paste0(criteria, collapse = ":")]]
braakGenes <- data.frame(
  entrez_id = braakGenes, 
  gene_symbol = entrezId2Name(braakGenes),
  r = labelCor[braakGenes, "r"], 
  r_BH = labelCor[braakGenes, "BH"], 
  fc = diffExpr[braakGenes, "meanDiff"], 
  fc_BH = diffExpr[braakGenes, "BH"]
)
braakGenes <- braakGenes[order(braakGenes$r),]
save(braakGenes, file = "resources/braakGenes.RData")
write.table(braakGenes, file = "braakGenes.txt", row.names = FALSE, quote = FALSE, sep = "\t")

########## Bar plot of selected genes ##########

# row labels
t1 <- paste("|r| >", floor(min(abs(braakGenes$r))*100)/100)
t2 <- paste("|FC| >", floor(min(abs(braakGenes$fc))*100)/100)
t3 <- paste("Pfc <", floor(as.numeric(max(braakGenes$fc_BH))*1e5)/1e5)
labels <- c(t1, t2, t3)
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
  theme_minimal() +
  theme(
    axis.text = element_text(size = 11),
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  )
p
pdf("braakgenes_barplot.pdf", 5, 2)
print(p)
dev.off()

########## Presence of PD-implicated genes ##########
pdGenes <- list(hiImpact = c("SNCA", "LRRK2", "GBA", "VPS35", "PARK2", "UCHL1", "PINK1", "PARK7", "ATP13A2", "PLA2G6", "FBXO7", "DNAJC6", "SYNJ1", 
                             "EIF4G1", "DNAJC13", "CHCHD2", "C20orf30", "RIC3", "LRP10"), #TMEM230 is C20orf30
                jansen2017 = c("INPP5F", "TMEM175", "ASH1L", "MAPT", "RIT1", "C14orf83", "STK39", "GPNMB", "BST1", 
                               "SIPA1L2", "DLG2", "NUCKS1", "GCH1", "MCCC1", "FAM47E", "BCKDK", "TMPRSS9", "UBOX5", 
                               "CCDC62", "SYNJ1", "EIF4G1", "FBXO7", "C20orf30", "POLG", "VPS13C", "PLA2G6"),
                hla = c("HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DQB1"),
                'Chang et al. 2017' = read.table("chang2017_riskgenes.txt", comment.char = "#", sep = "\n", row.names = NULL, stringsAsFactors = FALSE)[, 1], 
                'Nalls et al. 2014' = read.table("nalls2014_riskgenes.txt", comment.char = "#", sep = "\n", row.names = NULL, stringsAsFactors = FALSE)[, 1],
                parkingenes = read.table("parkin_genes_hgnc.txt", sep = "\t", header = TRUE)$Approved.Symbol
)
pdGenesID <- lapply(pdGenes, name2EntrezId)
tab <- lapply(names(pdGenesID), function(n){
  x <- pdGenesID[[n]]
  g <- intersect(braakGenes$entrez_id, x)
  cbind(study = rep(n, length(g)), braakGenes[braakGenes$entrez_id %in% g,])
})
tab <- Reduce(rbind, tab)
write.table(tab, file = "pdgenes_stats.txt", sep ="\t", quote = FALSE, row.names = FALSE)

########## Volcano plot for label correlation and differential expression ##########

braak_pos <- braakGenes$entrez_id[braakGenes$r>0]
braak_neg <- braakGenes$entrez_id[braakGenes$r<0]

# Plot theme
theme <- theme(legend.position = "none",
               panel.background = element_blank(),
               axis.line = element_line(colour = "black"),
               axis.title =  element_text(size = 12),
               plot.title = element_text(size = 12, face = "bold"))

prepare.tab <- function(tab){
  tab$'logp' <- -log10(tab$BH)
  tab$info <- sapply(rownames(tab), function(x){
    if (x %in% braak_pos) 1
    else if (x %in% braak_neg) 2
    else 0
  })
  tab$info <- as.factor(tab$info)
  order <- order(tab$info)# Plotting order of data points 
  tab <- tab[order, ]
  tab
}

# Braak correlation plot
tab <- prepare.tab(labelCor)
pdf(file = "volcanoplot_cor.pdf", 4, 3)
ggplot(tab, aes(r, logp, colour = info)) +
  geom_point(size = 0.25, alpha = 0.3) +
  scale_colour_manual(values = c("0"="grey", "1"="red", "2"="blue")) +
  labs(x = "r", y = "-log10 P-value") +
  ggtitle("Braak correlation") +
  theme
dev.off()

# Differential expression plot
ctab <- Reduce(rbind, summaryDiffExpr)
ctab$logp <- -log10(ctab$BH)
xmax <- max(ctab$meanDiff)
xmin <- min(ctab$meanDiff)
ymax <- ceiling(max(ctab$'logp'[is.finite(ctab$'logp')]))

plotll <- lapply(names(summaryDiffExpr), function(rp){
  tab <- summaryDiffExpr[[rp]]
  tab <- prepare.tab(tab)
  p <- ggplot(tab, aes(meanDiff, logp, colour = info)) +
    geom_point(size = 0.25, alpha = 0.3) +
    scale_colour_manual(values = c("0"="grey", "1"="red", "2"="blue")) +
    labs(x = "Fold-change", y = "-log10 P-value") +
    scale_x_continuous(limits = c(xmin, xmax), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, ymax), expand = c(0,0)) +
    ggtitle(paste("Braak stage ", gsub("braak", "", gsub("-", " vs ", rp)))) +
    theme
  p
  name <- paste0("DiffExpr_braak/volcanoplot_", rp, ".pdf")
  pdf(file = name, 4, 3)
  print(p)
  dev.off()
})

########## Heatmap expression of BRGs ##########

genes <- braakGenes$entrez_id
colsep <- which(braakGenes$r>0)[1]
colPal <- c("darkblue", "white", "darkred")
rampcols <- colorRampPalette(colors = colPal, space="Lab")(200)

pdf("heatmap_BRGs.pdf", 16, 12)
lapply(donorNames, function(d){
  # Subselect expression matrices
  samples <- unlist(braak_idx[[d]])
  df <- sampleInfo[[d]][samples,]
  labels <- braakLabels[[d]][samples]
  exprMat <- scale(t(brainExpr[[d]][genes, samples]))
  rowsep <- match(c(2:6), labels)# separate Braak regions
  
  heatmap.2(exprMat, col = rampcols, 
            labRow = df$acronym, labCol = entrezId2Name(colnames(exprMat)), 
            rowsep = rowsep, colsep = colsep, sepcolor = "black",
            Rowv=FALSE, Colv=FALSE, 
            cexCol = .5, cexRow = .5,
            scale = "none", trace = "none", dendrogram = "none", density.info = "none", 
            RowSideColors = df$color_hex_triplet, 
            main = paste0("Expression of BRGs in ", d),
            margins = c(5, 5))
})
dev.off()

# Functional enrichment of genes which expression correlate with Braak regions

braak <- lapply(c(positive = "pos", negative = "neg"), function(x){
  if (x == "pos") braakGenes$entrez_id[braakGenes$r > 0]
  else braakGenes$entrez_id[braakGenes$r < 0]
})

#Functional enrichment of  genes correlated greater or smaller than 0
david<-DAVIDWebService$new(email="D.L.Keo@tudelft.nl",
                           url="https://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
bg_list <- ahba.genes()
bg <- addList(david, bg_list, idType = "ENTREZ_GENE_ID", listName = "AHBA background", listType = "Background")
bg
t <- 0.05 # EASE p-value threshold
setTimeOut(david, 200000)

# Enrichment of positively and negatively correlated progression genes
lapply(names(braak), function(r){
  genes <- braak[[r]]
  result <- addList(david, genes, idType = "ENTREZ_GENE_ID", listName = r, listType = "Gene")
  print(result)
  setCurrentBackgroundPosition(david, 1)
  getFunctionalAnnotationChartFile(david, paste0("Functional_analyses/", r, "_goterms.txt"), threshold=t, count=2L)
  getClusterReportFile(david, paste0("Functional_analyses/", r, "_termclusters.txt"), type = c("Term"))
})
