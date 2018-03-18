# Volcano plot eigen gene braak label correlation
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")

library(ggplot2)
library(gridExtra)
library(ggrepel)

load("resources/summaryLabelCorrEG.RData")
labelCor <- do.call(rbind.data.frame, lapply(summaryLabelCorrEG, function(g) g["summary",]))

# Volcano plot
theme <- theme(legend.position = "none",
               panel.background = element_blank(),
               axis.line = element_line(colour = "black"),
               axis.title =  element_text(size = 12),
               plot.title = element_text(size = 12, face = "bold"))

tab <- labelCor
tab$benjamini_hochberg <- p.adjust(tab$pvalue, method = "BH")
tab$'logp' <- -log10(tab$benjamini_hochberg)
corrEG <- rownames(labelCor)[tab$benjamini_hochberg < 0.001] # significant correlated modules
tab$info <- as.numeric(rownames(tab) %in% corrEG)

braakMod <- read.table(file = "module_braakgene_overlap_braak1-6.txt", sep = "\t", header = TRUE)
braakMod <- braakMod$module[braakMod$benjamini_hochberg < 0.05] # significant overlap with braak genes
tab$labels <- rownames(tab)
tab$labels[!tab$labels %in% braakMod] <- ""
tab[braakMod, "info"] <- 2

# Plotting order of data points 
tab$info <- as.factor(tab$info)
order <- order(tab$info)
tab <- tab[order, ]

xmax <- max(tab$r)
xmin <- min(tab$r)
ymax <-  max(tab$'logp'[is.finite(tab$'logp')])
p <- ggplot(tab, aes(r, logp, colour = info)) +
  geom_point(size = 0.25) +
  scale_colour_manual(values = c("0"="#F8766D", "1"="#00BFC4", "2"="black")) +
  geom_text_repel(label = tab$labels, colour = "black", size = 3, nudge_x = 0.2) +
  labs(x = "Pearson's r", y = "-log10 p-value") +
  scale_x_continuous(limits = c(xmin, xmax), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, ymax), expand = c(0,0)) +
  ggtitle("Eigen gene-Braak label correlation") +
  theme

order <- order(tab$r)
label <- rownames(tab)[order]
color <- ifelse(label %in% corrEG, "#F8766D", "#00BFC4")
barplot(tab$r[order], horiz = TRUE, col = color)
hist(tab$r, breaks = 30)

name <- paste0("eigengene_r_volcanoplot.pdf")
pdf(file = name, 4, 3)
print(p)
dev.off()