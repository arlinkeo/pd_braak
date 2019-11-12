# Differential expression analysis between Braak stage-related genes in the Allen Human Brain Atlas
library(metafor)
library(reshape2)
library(plyr)

# Pairwise combinations of Braak regions R1-R6
roiPairs <- t(combn(braakRoi, 2))
rownames(roiPairs) <- apply(roiPairs, 1, paste, collapse = "-")
colnames(roiPairs) <- c("region_A", "region_B")

# T-test to get p-values and CI's (only needed for forest plot of meta-analysis)
ttest <- lapply(donorNames, function(d){
  print(d)
  expr <- brainExpr[[d]]
  exprll <- lapply(braak_idx[[d]], function(b){ # expr. in Braak regions 1-6
    expr[, b]
  })
  tab <- alply(roiPairs, 1, function(r){
    print(unname(r))
    a <- exprll[[r[1]]] # matrix group 1
    b <- exprll[[r[2]]] # matrix group 2
    t.test.table(a,b) # Test for all genes
  }, .dims = TRUE) # keep names
  simplify2array(tab) # 3D array: genes x measures x region pairs
})
ttest <- simplify2array(ttest) # 4D array: genes x measures x region pairs x donors
saveRDS(ttest, file = "output/ttest.rds")

# Print number of diff. expr. genes
apply(ttest, c(3,4), function(x){
  meanDiff <- x[, "meanA"] - x[, "meanB"]
  sum(x[, "BH"] < 0.05 & abs(meanDiff) > 1)
})

# Meta-analysis across donors
summaryDiffExpr <- aaply(ttest, c(1,3), function(g){ # For each Braak region pair and gene
  gene <- as.data.frame(t(g))
  t <- escalc(measure = "MD",  # Get estimates, variance (needed for meta-analysis) and confidence intervals (region B vs. A)
              m1i = gene[, "meanB"], m2i = gene[, "meanA"], # estimate
              n1i = gene[, "sizeB"], n2i = gene[, "sizeA"], 
              sd1i = sqrt(gene[, "varB"]), sd2i = sqrt(gene[, "varA"]))
  t <- summary(t)[, -c(3,4)]
  rownames(t) <- rownames(gene)
  colnames(t) <- c("Estimate", "Var", "lower95", "upper95")
  summary <- rma(t$Estimate, t$Var, method = "DL", test = "t")
  gene$weight <- weights(summary)
  t <- cbind(t, pvalue = gene[, "pvalue"], weight = weights(summary))
  t <- rbind(t, 'summary' = list(summary$beta, summary$se^2 , summary$ci.lb, summary$ci.ub,
                            summary$pval, sum(weights(summary))))
  as.matrix(t)
}) # 4D-array: genes x regions x donors x measures
saveRDS(summaryDiffExpr, file = "output/summaryDiffExpr.rds")

# Filter summary estimates, and correct P-values
summaryDiffExpr <- aaply(summaryDiffExpr, c(2,3), function(t){ # P-value corrected for genes
  b <- p.adjust(t[, "pvalue"], method = "BH")
  cbind(t, BH = b)
}) # region x donors x genes x measures

diffGenes <- aaply(summaryDiffExpr, c(1,2), function(t){
  down <-  t[which(t[, "Estimate"] < -1 & t[,"BH"] < 0.05), ]
  up <- t[which(t[, "Estimate"] > 1 & t[,"BH"] < 0.05), ]
  n1 <- nrow(down)
  n2 <- nrow(up)
  n1 <- ifelse(is.null(n1), 0, n1)
  n2 <- ifelse(is.null(n2), 0, n2)
  c(downregulated = n1, upregulated = n2)
})

########## Bar plot of differentially expressed genes between all Braak regions ##########

df <- melt(diffGenes[, "summary", ])
colnames(df) <- c("region", "dir", "size")
df$size[df$dir == "downregulated"] <- df$size[df$dir == "downregulated"]*-1
df$region <- factor(df$region, levels = rev(unique(df$region)))
lab_offset <- max(df$size)*0.5
df$y <- ifelse(df$dir == "upregulated", df$size + lab_offset, df$size - lab_offset)

p <- ggplot(df) +
  geom_col(aes(x=region, y = size, fill=dir), size = 0.5, colour = "black") +
  geom_text(aes(x=region, y= y, label=format(abs(df$size), big.mark=","))) +
  scale_fill_manual(values = c("blue", "red")) +
  scale_y_continuous(expand = c(0.1,0.1)) +
  coord_flip() +
  labs(x = "", y = "Number of differentially expressed genes") +
  theme(
    axis.text = element_text(size = 11),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    legend.title = element_blank(),
    legend.position = "top"
  )
pdf("output/diff_expr_barplot.pdf", 5, 4)
p
dev.off()
