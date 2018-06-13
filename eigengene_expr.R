# Eigen gene differential expression
setwd("C:/Users/dkeo/surfdrive/pd_braak")
library(ggplot2)
library(ggrepel)
source("PD/base_script.R")
load("../ABA_Rdata/BrainExprNorm.RData")
load("resources/modules.RData")
load("resources/braakLabels.RData")

##############################################################################################
# Functions

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

##############################################################################################

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

##############################################################################################

# Volcano plot
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
mod_pos <- rownames(tab)[tab$BH < 0.001 & tab$r > 0] # significant correlated modules
mod_neg <- rownames(tab)[tab$BH < 0.001 & tab$r < 0] # significant correlated modules
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
p
pdf("eigengene_r_volcanoplot.pdf", 6, 4.5)
p
dev.off()