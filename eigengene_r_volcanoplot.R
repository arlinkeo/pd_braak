# Volcano plot eigen gene braak label correlation
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")

library(ggplot2)
library(gridExtra)
library(ggrepel)
library(ggpubr)
load("resources/module_enrichment.RData")
load("resources/summaryLabelCorrEG.RData")
labelCor <- do.call(rbind.data.frame, lapply(summaryLabelCorrEG, function(g) g["summary",]))

# Volcano plot
theme <- theme(legend.position = "none",
               panel.background = element_blank(),
               axis.line = element_line(colour = "black"),
               axis.title =  element_text(size = 16),
               plot.title = element_text(size = 16),
               axis.text = element_text(size = 16),
               axis.title.x = element_text(face="italic")
               )

tab <- labelCor
tab$benjamini_hochberg <- p.adjust(tab$pvalue, method = "BH")
tab$'logp' <- -log10(tab$benjamini_hochberg)
corrEG <- rownames(labelCor)[tab$benjamini_hochberg < 0.001] # significant correlated modules
tab$info <- as.numeric(rownames(tab) %in% corrEG)

braakMod <- module_enrichment$module[module_enrichment$pvalue_progression < 0.001]
tab$labels <- module_enrichment$module
tab$labels[!tab$labels %in% braakMod] <- ""
tab[tab$labels %in% braakMod, "info"] <- 2
tab$module <- paste0("M", rownames(tab))

# Plotting order of data points 
tab$info <- as.factor(tab$info)
order <- order(tab$info)
tab <- tab[order, ]

xmax <- max(tab$r)+.2
xmin <- min(tab$r)-.2
ymax <-  ceiling(max(tab$'logp'))
p1 <- ggplot(tab, aes(r, logp, colour = info)) +
  geom_point(size = 4) +
  scale_colour_manual(values = c("0"="#F8766D", "1"="#00BFC4", "2"="black")) +
  geom_text_repel(label = tab$labels, colour = "black", size = 6, nudge_x = 0) +
  labs(y = "-log10 p-value") +
  scale_x_continuous(limits = c(xmin, xmax), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, ymax), expand = c(0,0)) +
  annotate("label", label ="Downregulated", x=0.7, y=0.5) +
  annotate("label", label ="Upregulated", x=-0.7, y=0.5) +
  theme
p1
####################################################################################

# Braak label correlation plot
orderEG <- order(tab$r)
tab <- tab[orderEG, ]

# significance stars
star <- function(v){
  sapply(v, function(x) 
    if (x<0.001) "*" 
    # else if (x<0.01) "**" 
    # else if (x<0.05) "*" 
    else ""
  )
}

tab$star <- star(tab$benjamini_hochberg)
tab$module <- factor(tab$module, levels = unique(tab$module))
offset <- 0.1
tab$y <- sapply(tab$r, function(x) if (x>0) x+offset else x-offset)

colPal <- c("darkblue", "white", "darkred")
rampcols <- colorRampPalette(colors = colPal, space="Lab")(201)
tab$color <- rampcols[as.numeric(cut(tab$r, breaks = 201))]
tab$color<- factor(tab$color, levels = unique(tab$color))

p2 <- ggplot(tab) + 
  geom_col(aes(x=module, y = r, fill=r)) +
  geom_text(aes(x=module, y=y, label=star), size = 4, vjust = 0.75) +
  scale_y_continuous(limits = c(xmin, xmax), expand = c(0,0)) +
  labs(y ="Pearson's r") +
  scale_fill_gradientn(colours = rampcols) +
  theme + theme(axis.text.y=element_blank(), axis.ticks.y = element_blank()) +
  coord_flip()
p2

p <- ggarrange(p1, p2+coord_flip(), nrow = 2, heights = c(4, 3),  labels = c("A", "B"), align = "v")
p
pdf("eigengene_r_volcanoplot.pdf", 8, 6)
print(p)
p1
dev.off()
