# Basic ggplot theme# Forest plot of differential expression in each brain for a single gene
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(ggplot2)
library(gridExtra)
library(gtable)
load("resources/summaryDiffExpr.RData")
load("resources/summaryLabelCorr.RData")

gene <- "SNCA"
geneId <- name2EntrezId(gene)
theme.grid <- theme(panel.background = element_blank(), panel.grid = element_blank(), axis.ticks.y = element_blank(), axis.line.x = element_line(colour = "black"), 
                    legend.position = "none", axis.title.y = element_blank())

####################################
# Forest plot expression change, multiple region pairs

tab <- summaryDiffExpr$`braak1-braak6`[[geneId]]

tab$donors <- factor(tab$donors, levels = unique(rev(tab$donors)))
tab$rmd <- paste0(round(tab$meanDiff, digits = 2), " (", round(tab$lower95, digits = 2), ", ", round(tab$upper95, digits = 2), ")")
tab$pvalue <- signif(tab$pvalue, digits = 2) # Uncorrected
tab$is.summary <- tab$donors == "Summary"

# Forest plot
table <- tableGrob(tab[, c("rmd", "weight")], rows = NULL, theme = ttheme_minimal())
table$heights <- unit(rep(1/(nrow(table)), nrow(table)), "npc")

intercepts <- if (min(tab$meanDiff) > 0) 1.53 else if (max(tab$meanDiff) < 0) -1.53 else c(-1.53,1.53)

plot <- ggplot(data = tab, aes(meanDiff, donors))  +
  geom_point(aes(size = weight, shape = is.summary, color = is.summary)) +
  geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = intercepts, linetype = "dotted") +
  scale_shape_manual(values = c(15,18)) +
  scale_color_manual(values = c('#00CCCC', '#FF8000')) +
  labs(title = bquote(italic(.(gene))), x = "Fold-change") +
  theme.grid
g=ggplotGrob(plot)
g$heights[7] = grobHeight(1)
grid.arrange(g, table, ncol = 2)

pdf(file = paste0("forestplot_", gene, ".pdf"), 12, 6)
leftPanel
dev.off()

####################################
# Forest plot expression change, one region pair

tab <- summaryDiffExpr$`braak1-braak6`[[geneId]]

tab$donors <- factor(tab$donors, levels = unique(rev(tab$donors)))
# Side table info
tab$rmd <- paste0(round(tab$meanDiff, digits = 2), " (", round(tab$lower95, digits = 2), ", ", round(tab$upper95, digits = 2), ")")
tab$pvalue <- signif(tab$pvalue, digits = 2)
tab$is.summary <- tab$donors == "Summary"

# Forest plot
fp <- ggplot(data = tab, aes(meanDiff, donors)) 

leftPanel <- fp +
  geom_point(aes(size = weight, shape = is.summary, color = is.summary)) +
  geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  # geom_vline(xintercept = 1, linetype = "dotted") +
  geom_vline(xintercept = -1, linetype = "dotted") +
  scale_shape_manual(values = c(15,18)) +
  scale_color_manual(values = c('#00CCCC', '#FF8000')) +
  labs(title = bquote(italic(.(gene))), x = "Mean difference") +
  scale_y_discrete(labels = rev(tab$donor)) +
  # scale_x_continuous(limits = c(x.min, x.max)) +
  geom_text(aes(x =  1, label = rmd), size = 3) +
  geom_text(aes(x = 1.9, label = weight), size = 3) +
  # geom_text(data = colLabels, aes(x, y, label = label)) +
  theme.grid + theme(plot.margin = unit(c(0,4,0,4), "lines"))
pdf(file = paste0("forestplot_b1b6_", gene, ".pdf"), 8,2)
leftPanel
dev.off()

##########################################

# Forrest plot Braak correlations



#Get values for a gene to plot
tab <- summaryLabelCorr[[geneId]]
tab$donors <- factor(tab$donors, levels = unique(rev(tab$donors)))
# Side table info
tab$rci <- paste0(round(tab$r, digits = 2), " (", round(tab$lower95, digits = 2), ", ", round(tab$upper95, digits = 2), ")")
tab$is.summary <- tab$donors == "Summary"

# Forest plot
fp <- ggplot(data = tab, aes(r, donors)) 

# x.max <- round(max(tab$upper95), digits = 2)
# x.min <- round(min(tab$lower95), digits = 2)
# x.positions <- c(0.25, 0.55, 0.65) + x.max#c(1:3)
# colLabels <- data.frame(x = x.positions, y = x.positions, label = c("Raw mean difference", "N", "Weight"))

leftPanel <- fp +
  geom_point(aes(size = braakSize, shape = is.summary, color = is.summary)) +
  geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dotted") +
  geom_vline(xintercept = -1, linetype = "dotted") +
  # geom_vline(xintercept = -1, linetype = "dotted") +
  scale_shape_manual(values = c(15,18)) +
  scale_color_manual(values = c('#00CCCC', '#FF8000')) +
  labs(title = bquote(italic(.(gene))), x = "Correlation") +
  scale_y_discrete(labels = rev(tab$donor)) +
  # scale_x_continuous(limits = c(x.min, x.max)) +
  geom_text(aes(x = 1.6, label = rci), size = 3) +
  geom_text(aes(x = 2.2, label = weight), size = 3) +
  # geom_text(data = colLabels, aes(x, y, label = label)) +
  theme.grid + theme(plot.margin = unit(c(0,4,0,4), "lines"))
pdf(file = paste0("forestplot_corr_", gene, ".pdf"), 8, 2)
leftPanel
dev.off()