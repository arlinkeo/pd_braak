# Eigen gene differential expression
setwd("C:/Users/dkeo/surfdrive/pd_braak")
library(ggplot2)
library(ggrepel)
library(gplots)
source("PD/base_script.R")
load("resources/modules.RData")
load("resources/braakInfo.RData")

brainExpr <- readRDS("../AHBA_Arlin/gene_expr.RDS")

#####  Functions #####

# eigen gene function for data matrix (samples x genes)
eigen.gene <- function(x){
  eg <- prcomp(x, scale. = TRUE)$x[, 1]# 1st PC (eigen gene expr)
  mean <- apply(x, 1, mean)
  if (cor(eg, mean) > 0) eg else -eg # flip sign of eigen gene based on the data
}

# Braak-correlation for eigengenes
summary.braak.cor <- dget("PD/summary.braak.cor.R")

##### PCA eigen gene #####

# PCA first component of subselection expr. matrices
eigenExpr <- lapply(donorNames, function(d){
  s <- unlist(braak_idx[[d]])
  expr <- brainExpr[[d]][, s]
  as.data.frame(t(sapply(modules, function(genes){ # For each module with genes (grouped gene rows)
    eigen.gene(t(expr[genes, ]))
  })))
})
save(eigenExpr, file = "resources/eigenExpr.RData")

# Summary Braak correlation
labels <- lapply(donorNames, function(d){
  s <- unlist(braak_idx[[d]])
  braakLabels[[d]][s]
})
summaryLabelCorrEG <- summary.braak.cor(eigenExpr, labels)
save(summaryLabelCorrEG, file = "resources/summaryLabelCorrEG.RData")

labelCor <- do.call(rbind.data.frame, lapply(summaryLabelCorrEG, function(g) g["summary",]))
labelCor$BH <- p.adjust(labelCor$pvalue, method = "BH")
orderEG <- order(labelCor$r)
labelCor <- labelCor[orderEG, ]

mod_neg <- rownames(labelCor)[labelCor$BH < 0.0001 & labelCor$r < 0] # significant correlated modules
mod_pos <- rownames(labelCor)[labelCor$BH < 0.0001 & labelCor$r > 0] # significant correlated modules

braakModules <- list(down = mod_neg, up = mod_pos)
save(braakModules, file = "resources/braakModules.RData")

#####  Volcano plot #####

theme <- theme(legend.position = "none",
               panel.background = element_blank(),
               axis.line = element_line(colour = "black"),
               axis.title =  element_text(size = 16),
               plot.title = element_text(size = 16),
               axis.text = element_text(size = 16)
)

tab <- labelCor
tab$'logp' <- -log10(tab$BH)
# eg <- rownames(tab)[tab$BH < 0.001 & tab$r < 0] # significant correlated modules
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
ymax <-  ceiling(max(tab$'logp'))+.5

p <- ggplot(tab, aes(r, logp, colour = info)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_text_repel(aes(label=label), colour = "black", size = 4, nudge_x = 0) +
  scale_colour_manual(values = c("0"="grey", "1"="blue", "2"="red")) +
  labs(x = "Correlation with Braak r", y = "-log10 p-value") +
  scale_x_continuous(limits = c(xmin, xmax), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, ymax), expand = c(0,0)) +
  theme

pdf("eigengene_r_volcanoplot.pdf", 5, 5)
p
dev.off()

##### Heatmap Expression of eigen gene #####

# Order of modules
cols <- rownames(labelCor)

# Heatmap colors
colPal <- c("darkblue", "white", "darkred")
rampcols <- colorRampPalette(colors = colPal, space="Lab")(201)
colColor <- rampcols[as.numeric(cut(labelCor$r, breaks = 201))]
colsep <- tail(which(labelCor$r < 0), 1)

pdf("heatmap_expr_eigengenes.pdf", 5, 12)
lapply(donorNames[1], function(d){
  samples <- unlist(braak_idx[[d]])
  df <- sampleInfo[[d]][samples,]
  labels <- braakLabels[[d]][samples]
  exprMat <- t(as.matrix(eigenExpr[[d]][cols, ]))
  rowsep <- match(c(2:6), labels) -1# separate Braak regions
  rowColor <- df$color_hex_triplet
  # rownames(exprMat) <- NULL
  heatmap.2(exprMat, col = rampcols, 
            labRow = df$name, 
            rowsep = rowsep, colsep = colsep, sepcolor = "black",
            Rowv=FALSE, Colv=FALSE, 
            cexCol = 1, cexRow = .1,
            scale = "none", trace = "none", dendrogram = "none", key = FALSE, 
            RowSideColors = rowColor, ColSideColors = colColor,
            main = d,
            margins = c(20,5))
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
