# Braak genes plots

# Volcano plot with fold-change (x) and correlation (y)
tab <- data.frame(r = labelCor$r, fc = diffExpr$Estimate)
rownames(tab) <- rownames(labelCor)
tab <- prepare.tab(tab)
pdf("BRGs_correlation_vs_foldchange.pdf", 4, 3)
ggplot(tab, aes(r, fc, colour = info)) +
  geom_point(size = 0.25, alpha = 0.3) +
  scale_colour_manual(values = c("0"="grey", "1"="red", "2"="blue")) +
  labs(x = "Braak correlation", y = "FC (R1-R6)") +
  ggtitle("Braak correlation") +
  theme
dev.off()

# Differential expression volcano plots for each pair or regions
xmax <- max(summaryDiffExpr[, , "Estimate"])
xmin <- min(summaryDiffExpr[, , "Estimate"])
ymax <- ceiling(-log10(min(summaryDiffExpr[, , "BH"])))

plotll <- lapply(dimnames(summaryDiffExpr)[[1]], function(rp){
  tab <- data.frame(summaryDiffExpr[rp,,])
  tab <- prepare.tab(tab)
  tab$'logp' <- -log10(tab$BH)
  p <- ggplot(tab, aes(Estimate, logp, colour = info)) +
    geom_point(size = 0.25, alpha = 0.3) +
    scale_colour_manual(values = c("0"="grey", "1"="red", "2"="blue")) +
    labs(x = "Fold-change", y = "-log10 P-value") +
    scale_x_continuous(limits = c(xmin, xmax), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, ymax), expand = c(0,0)) +
    ggtitle(gsub("-", " vs ", rp)) +
    theme
  p
  name <- paste0("DiffExpr_braak/volcanoplot_", rp, ".pdf")
  pdf(file = name, 4, 3)
  print(p)
  dev.off()
})

########## Heatmap expression of BRGs ##########

genes <- braakGenes$entrez_id # is sorted by correlation
colsep <- which(braakGenes$r>0)[1]
colPal <- c("darkblue", "white", "darkred")
rampcols <- colorRampPalette(colors = colPal, space="Lab")(200)

pdf("heatmap_BRGs.pdf", 16, 12)
lapply(donorNames, function(d){
  # Subselect expression matrices
  samples <- unlist(braak_idx[[d]])
  df <- sampleInfo[[d]][samples,]
  labels <- braakLabels[[d]][samples]
  # exprMat <- scale(t(brainExpr[[d]][genes, samples]))
  exprMat <- scale(t(brainExpr[[d]][name2EntrezId(c("SNCA", "GCH1", "TH", "SLC6A3", "SLC18A2")), samples]))
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