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

# Function for eigen gene expression for each module in the same region
# x: data, l: list of modules with genes, s: logical vector of samples (columns)
eigen.data <- function(x, l, s){ 
  df <- as.data.frame(t(sapply(l, function(genes){ # For each module with genes (grouped gene rows)
    eigen.gene(t(x[genes, s]))
  })))
  colnames(df) <- names(s)[s]
  df
}

# Braak-correlation for eigengenes
summary.braak.cor <- dget("PD/summary.braak.cor.R")

##### PCA eigen gene #####

# PCA first component of subselection expr. matrices
eigenExpr <- lapply(donorNames, function(d){
  s <- braakLabels[[d]] != 0 # logical vector
  expr <- brainExprNorm[[d]]
  eigen.data(expr, modules, s)
})
save(eigenExpr, file = "resources/eigenExpr.RData")

# Summary Braak correlation
labels <- lapply(braakLabels, function(labels){labels[labels != "0"]})
summaryLabelCorrEG <- summary.braak.cor(eigenExpr, labels)
save(summaryLabelCorrEG, file = "resources/summaryLabelCorrEG.RData")

#####  Volcano plot #####

theme <- theme(legend.position = "none",
               panel.background = element_blank(),
               axis.line = element_line(colour = "black"),
               axis.title =  element_text(size = 16),
               plot.title = element_text(size = 16),
               axis.text = element_text(size = 16),
               axis.title.x = element_text(face="italic")
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

xmax <- max(tab$r)+.2
xmin <- min(tab$r)-.2
ymax <-  ceiling(max(tab$'logp'))

p <- ggplot(tab, aes(r, logp, colour = info)) +
  geom_point(size = 2, alpha = 0.5) +
  # geom_text(aes(label=label)) +
  geom_text_repel(aes(label=label), colour = "black", size = 4, nudge_x = 0) +
  scale_colour_manual(values = c("0"="grey", "1"="blue", "2"="red")) +
  labs(y = "-log10 p-value") +
  scale_x_continuous(limits = c(xmin, xmax), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, ymax), expand = c(0,0)) +
  theme

pdf("eigengene_r_volcanoplot.pdf", 5, 4.5)
p
dev.off()

##### Heatmap Expression of eigen gene #####

ontology <- read.csv("../ABA_human_processed/Ontology_edited.csv")

# IDs and AHBA colors for each sample per donor
sampleInfo <- lapply(donorNames, function(d){
  sampleIds <- read.csv(paste("../ABA_human_processed/sample_info_normalized_microarray_", d, "_2014-11-11.csv", sep = ""))[ , 1]
  info <- ontology[match(sampleIds, ontology$id), ]
  info$color_hex_triplet <- sapply(info$color_hex_triplet, function(c){
    if (nchar(c) == 5) {paste("#0", c, sep = "")} else {paste("#", c, sep = "")}
  })
  info
})

# modules <- names(modules[["braak1-6"]])
colOrder <- order(labelCor$r)
colPal <- c("darkblue", "white", "darkred")
rampcols <- colorRampPalette(colors = colPal, space="Lab")(201)
colColor <- rampcols[as.numeric(cut(labelCor$r, breaks = 201))][colOrder]
colsep <- which(labelCor$r[colOrder] > 0)[1]

pdf("heatmap_expr_eigengenes.pdf", 8, 6)
lapply(donorNames, function(d){
  # Subselect expression matrices
  samples <- braakLabels[[d]] != 0
  df <- sampleInfo[[d]][samples,]
  exprMat <- as.matrix(eigenExpr[[d]])
  rowColor <- df$color_hex_triplet
  
  rowOrder <- order(labels[[d]], -df$graph_order)
  df <- df[rowOrder, ]
  exprMat <- t(exprMat[colOrder, rowOrder])
  colnames(exprMat) <- NULL
  rowColor <- rowColor[rowOrder]
  rowsep <- match(c(2:6), labels[[d]][rowOrder])# separate Braak regions
  
  heatmap.2(exprMat, col = rampcols, 
            labRow = df$name, 
            rowsep = rowsep, colsep = colsep, sepcolor = "red",
            Rowv=FALSE, Colv=FALSE, 
            cexCol = 1, cexRow = .1,
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