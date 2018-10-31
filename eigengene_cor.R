# Eigen gene differential expression
setwd("C:/Users/dkeo/surfdrive/pd_braak")
library(ggplot2)
library(ggrepel)
library(gplots)
source("PD/base_script.R")
load("../ABA_Rdata/BrainExprNorm.RData")
load("resources/modules.RData")
load("resources/braakInfo.RData")

#####  Functions #####

# eigen gene function for data matrix (samples x genes)
eigen.gene <- function(x){
  eg <- prcomp(x)$x[, 1]# 1st PC (eigen gene expr)
  mean <- apply(x, 1, mean)
  if (cor(eg, mean) > 0) eg else -eg # flip sign of eigen gene based on the data
}

# # Function for eigen gene expression for each module in the same region
# # x: data, l: list of modules with genes, s: logical vector of samples (columns)
# eigen.data <- function(x, l){
#   df <- as.data.frame(t(sapply(l, function(genes){ # For each module with genes (grouped gene rows)
#     eigen.gene(t(x[genes, s]))
#   })))
#   colnames(df) <- names(s)[s]
#   df
# }

# Braak-correlation for eigengenes
summary.braak.cor <- dget("PD/summary.braak.cor.R")

##### PCA eigen gene #####

# PCA first component of subselection expr. matrices
eigenExpr <- lapply(donorNames, function(d){
  s <- unlist(braak_idx[[d]])
  expr <- brainExprNorm[[d]][, s]
  df <- as.data.frame(t(sapply(modules, function(genes){ # For each module with genes (grouped gene rows)
    eigen.gene(t(expr[genes, ]))
  })))
  # colnames(df) <- names(s)
  df
  # eigen.data(expr, modules)
})
save(eigenExpr, file = "resources/eigenExpr.RData")

# Summary Braak correlation
# labels <- lapply(braakLabels, function(labels){labels[labels != "0"]})
labels <- lapply(donorNames, function(d){
  braakLabels[[d]][unlist(braak_idx[[d]])]
})
summaryLabelCorrEG <- summary.braak.cor(eigenExpr, labels)
save(summaryLabelCorrEG, file = "resources/summaryLabelCorrEG.RData")

#####  Volcano plot #####

theme <- theme(legend.position = "none",
               panel.background = element_blank(),
               axis.line = element_line(colour = "black"),
               axis.title =  element_text(size = 16),
               plot.title = element_text(size = 16),
               axis.text = element_text(size = 16),
               axis.title.y = element_text(face="italic")
)

labelCor <- do.call(rbind.data.frame, lapply(summaryLabelCorrEG, function(g) g["summary",]))
tab <- labelCor
tab$BH <- p.adjust(tab$pvalue, method = "BH")
tab$'logp' <- -log10(tab$BH)
mod_neg <- rownames(tab)[tab$BH < 0.001 & tab$r < 0] # significant correlated modules
mod_pos <- rownames(tab)[tab$BH < 0.001 & tab$r > 0] # significant correlated modules
eg <- rownames(tab)[tab$BH < 0.001 & tab$r < 0] # significant correlated modules
tab$info <- sapply(rownames(tab), function(m){
  if (m %in% mod_neg) 1
  else if (m %in% mod_pos) 2
  else 0
})
tab$label <- rownames(tab)
tab$label[!tab$label %in% c(mod_neg, mod_pos)] <- ""

# Plotting order of data points 
tab$info <- as.factor(tab$info)
order <- order(tab$info)
tab <- tab[order, ]

ymax <- max(tab$r)+.5
ymin <- min(tab$r)-.5
xmax <-  ceiling(max(tab$'logp'))

p <- ggplot(tab, aes(logp, r, colour = info)) +
  geom_point(size = 2, alpha = 0.5) +
  # geom_text(aes(label=label)) +
  geom_text_repel(aes(label=label), colour = "black", size = 4, nudge_x = 0) +
  scale_colour_manual(values = c("0"="grey", "1"="blue", "2"="red")) +
  labs(x = "-log10 p-value") +
  scale_y_continuous(limits = c(ymin, ymax), expand = c(0,0)) +
  scale_x_continuous(limits = c(0, xmax), expand = c(0,0)) +
  theme

pdf("eigengene_r_volcanoplot.pdf", 5, 4.5)
p
dev.off()

##### Heatmap Expression of eigen gene #####

# Order of modules
rowOrder <- order(-labelCor$r)

# Heatmap colors
colPal <- c("darkblue", "white", "darkred")
rampcols <- colorRampPalette(colors = colPal, space="Lab")(201)
rowColor <- rampcols[as.numeric(cut(labelCor$r, breaks = 201))][rowOrder]
rowsep <- tail(which(labelCor$r[rowOrder] > 0), 1)

pdf("heatmap_expr_eigengenes.pdf", 10, 5)
lapply(donorNames, function(d){
  samples <- unlist(braak_idx[[d]])
  df <- sampleInfo[[d]][samples,]
  labels <- braakLabels[[d]][samples]
  exprMat <- as.matrix(eigenExpr[[d]][rowOrder, ])
  colsep <- match(c(2:6), labels) -1# separate Braak regions
  colColor <- df$color_hex_triplet
  # rownames(exprMat) <- NULL
  heatmap.2(exprMat, col = rampcols, 
            labCol = df$name, 
            rowsep = rowsep, colsep = colsep, sepcolor = "black",
            Rowv=FALSE, Colv=FALSE, 
            cexCol = .1, cexRow = 1,
            scale = "none", trace = "none", dendrogram = "none", key = FALSE, 
            RowSideColors = rowColor, ColSideColors = colColor,
            main = d,
            margins = c(5, 20))
})
dev.off()

##### Print table with module info #####

module_info <- data.frame(
  Module = names(modules),
  Size = sapply(modules, length), 
  r = round(labelCor$r, digits = 2),
  "BH-corrected P" = format(labelCor$pvalue, digits = 3, scientific = TRUE),
  genes = sapply(modules, function(m) paste0(entrezId2Name(m), collapse = ","))
)
write.table(module_info, file = "module_info.txt", sep = "\t", quote = FALSE, row.names = FALSE)
