# select significant genes based on significant summary estimate
library(metap)
library(gplots)
library(RDAVIDWebService)

########## Prepare data ##########

# Select results from differential expression and correlation analysis with summary statistics
diffExpr <- as.data.frame(summaryDiffExpr["R1-R6", "summary", ,])

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
order <- rev(order(abs(diffExpr$Estimate))) # order absolute corr.
diffGenes1 <- rownames(diffExpr)[order[1:top10]] # top 10% genes
diffGenes1_pos <- diffGenes1[diffExpr[diffGenes1, "Estimate"] > 0]
diffGenes1_neg <- diffGenes1[diffExpr[diffGenes1, "Estimate"] < 0]
max(abs(diffExpr[diffGenes1, "Estimate"]))
min(abs(diffExpr[diffGenes1, "Estimate"]))

# Top 10% significant Diff. expressed genes braak 1 vs. braak 6
order <- order(diffExpr$BH) # order absolute corr.
diffGenes2 <- rownames(diffExpr)[order[1:top10]] # top 10% genes
diffGenes2_pos <- diffGenes2[diffExpr[diffGenes2, "Estimate"] > 0]
diffGenes2_neg <- diffGenes2[diffExpr[diffGenes2, "Estimate"] < 0]
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
  gene_symbol = entrezId2Name(braakGenes),
  entrez_id = braakGenes, 
  r = labelCor[braakGenes, "r"], 
  r_BH = labelCor[braakGenes, "BH"], 
  fc = diffExpr[braakGenes, "Estimate"], 
  fc_BH = diffExpr[braakGenes, "BH"]
)
braakGenes <- braakGenes[order(braakGenes$r),]

table <- braakGenes
table[, c(3,5)] <- round(table[, c(3,5)], digits = 2)
table[, c(4,6)] <- format(table[, c(4,6)], digits = 3, scientific = TRUE)
write.table(table, file = "output/braakGenes.txt", row.names = FALSE, quote = FALSE, sep = "\t")

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
tab$y <- ifelse(tab$dir == "positive", tab$size+300, tab$size-300)

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
    panel.grid = element_blank(),
    legend.position = "top", 
    legend.title = element_blank()
  )
pdf("output/braakgenes_barplot.pdf", 4, 2)
p
dev.off()

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
tab$'logp' <- -log10(tab$BH)
pdf(file = "output/volcanoplot_cor.pdf", 4, 3)
ggplot(tab, aes(r, logp, colour = info)) +
  geom_point(size = 0.25, alpha = 0.3) +
  scale_colour_manual(values = c("0"="grey", "1"="red", "2"="blue")) +
  labs(x = bquote("Correlation with Braak stages ("*italic(r)*")"), y = bquote("-"*log[10]*" "*italic(P)*"-value")) +
  ggtitle("Braak correlation") +
  theme
dev.off()

# Functional enrichment of braak genes

braak <- lapply(c(positive = "pos", negative = "neg"), function(x){
  if (x == "pos") braakGenes$entrez_id[braakGenes$r > 0]
  else braakGenes$entrez_id[braakGenes$r < 0]
})

david<-DAVIDWebService$new(email="D.L.Keo@tudelft.nl",
                           url="https://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
setTimeOut(david, 200000)
bg_list <- ahba.genes()
bg <- addList(david, bg_list, idType = "ENTREZ_GENE_ID", listName = "AHBA background", listType = "Background")
bg
t <- 0.05 # EASE p-value threshold

# Enrichment of positively and negatively correlated progression genes
lapply(names(braak), function(r){
  genes <- braak[[r]]
  result <- addList(david, genes, idType = "ENTREZ_GENE_ID", listName = r, listType = "Gene")
  print(result)
  setCurrentBackgroundPosition(david, 1)
  getFunctionalAnnotationChartFile(david, paste0("output/Functional_analyses/", r, "_goterms.txt"), threshold=t, count=2L)
  getClusterReportFile(david, paste0("output/Functional_analyses/", r, "_termclusters.txt"), type = c("Term"))
})
