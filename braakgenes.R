# select significant genes based on significant summary estimate
library(metap)
library(gplots)
library(RDAVIDWebService)
library(ComplexHeatmap)
library(circlize)

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
pdf("output/brgs_venndiagram.pdf", 4.5, 4.5)
venn <- venn(ll)
venn1 <- venn(ll1)
venn2 <- venn(ll2)
dev.off()

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
table[, c(3)] <- format(table[, c(3)], digits = 2)
table[, c(5)] <- round(table[, c(5)], digits = 2)
table[, c(4,6)] <- format(table[, c(4,6)], digits = 3, scientific = TRUE)
colnames(table) <- c("Gene symbol", "Entrez ID", "Correlation with Braak r", 
"P-value of correlation with Braak r (BH-corrected)", "Fold-change between R1 and R6", 
"P-value of Fold-change between R1 and R6 (BH-corrected)")
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

# Braak volcano plot with correlations
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

bg <- list(
  'r < 0' = braakGenes$entrez_id[braakGenes$r < 0],
  'r > 0' = braakGenes$entrez_id[braakGenes$r > 0]
)

david<-DAVIDWebService$new(email="D.L.Keo@tudelft.nl",
                           url="https://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
setTimeOut(david, 200000)
background_list <- ahba.genes()
background <- addList(david, background_list, idType = "ENTREZ_GENE_ID", listName = "AHBA background", listType = "Background")
background
t <- 0.05 # EASE p-value threshold

# Enrichment of positively and negatively correlated progression genes
lapply(names(bg), function(r){
  genes <- bg[[r]]
  result <- addList(david, genes, idType = "ENTREZ_GENE_ID", listName = r, listType = "Gene")
  print(result)
  setCurrentBackgroundPosition(david, 1)
  if (r == "r < 0") r = "negative" else r = "positive"
  getFunctionalAnnotationChartFile(david, paste0("output/Functional_analyses/", r, "_goterms.txt"), threshold=t, count=2L)
  getClusterReportFile(david, paste0("output/Functional_analyses/", r, "_termclusters.txt"), type = c("Term"))
})

# Read files
lapply(names(bg), function(r){
  if (r == "r < 0") r = "negative" else r = "positive"
  go <- read.delim(paste0("output/Functional_analyses/", r, "_goterms.txt"))
  colnames(go) <- gsub("\\.", " ", colnames(go))
  colnames(go)[4] <- "%"
  print(paste0(r, ": ",  nrow(go)))
  # Select terms with BH-corrected P < 0.05
  go <- go[which(go$Benjamini < 0.05), ]
  print(paste0(r, ": ",  nrow(go)))
  # Remove terms with less than 20 counts
  go <- go[which(go$Count >= 20), ]
  print(paste0(r, ": ",  nrow(go)))
  write.table(go, file = paste0("output/Functional_analyses/", r, "_goterms_BH_count20.txt"), 
              quote = FALSE, sep = "\t", row.names = FALSE)
})

########## Heatmap expression of BRGs ##########

# Same heat colors for all BRG heatmaps
col_fun <- colorRamp2(c(-1.5, 0, 1.5), c("blue", "#EEEEEE", "red"))

expr <- lapply(donorNames, function(d){
  # Subselect expression matrices
  samples <- braak_idx[[d]]
  expr <- sapply(samples, function(s){ # Genes are sorted by correlation
    e <- brainExpr[[d]][unlist(bg), s]
    apply(e, 1, mean) # mean across samples in Braak region
   })
})
expr <- apply(simplify2array(expr), c(1,2), mean)
expr <- t(scale(t(expr))) # expr. is scaled across samples

pdf("output/heatmap_expr_BRGs_AHBA.pdf", 2.7, 10)
Heatmap(expr, name = 'Z-Score\nexpression',
        col = col_fun,
        row_split = rep(names(lengths(bg)), lengths(bg)),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        column_names_gp = gpar(fontsize = 10),
        column_names_rot = 0,
        column_names_centered = TRUE,
        row_title_rot = 0,
        width = unit(ncol(expr), "lines"),
        height = unit(nrow(expr)*.05, "lines")
)
dev.off()
