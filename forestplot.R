# Basic ggplot theme# Forest plot of differential expression in each brain for a single gene
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(ggplot2)
library(gridExtra)
library(gtable)
load("resources/summaryDiffExpr.RData")
load("resources/summaryLabelCor.RData")

gene <- "SNCA"
geneId <- name2EntrezId(gene)
theme.grid <- theme(panel.background = element_blank(), panel.grid = element_blank(), axis.ticks.y = element_blank(), axis.line.x = element_line(colour = "black"), 
                    legend.position = "none", axis.title.y = element_blank())

####################################
# # Forest plot expression change, multiple region pairs
# 
# tab <- summaryDiffExpr
# tab$donors <- factor(tab$donors, levels = unique(rev(tab$donors)))
# tab$rmd <- paste0(round(tab$meanDiff, digits = 2), " (", round(tab$lower95, digits = 2), ", ", round(tab$upper95, digits = 2), ")")
# tab$pvalue <- signif(tab$pvalue, digits = 2) # Uncorrected
# tab$is.summary <- tab$donors == "Summary"
# 
# # Forest plot
# table <- tableGrob(tab[, c("rmd", "weight")], rows = NULL, theme = ttheme_minimal())
# table$heights <- unit(rep(1/(nrow(table)), nrow(table)), "npc")
# 
# intercepts <- if (min(tab$meanDiff) > 0) 1.53 else if (max(tab$meanDiff) < 0) -1.53 else c(-1.53,1.53)
# 
# plot <- ggplot(data = tab, aes(meanDiff, donors))  +
#   geom_point(aes(size = weight, shape = is.summary, color = is.summary)) +
#   geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0) +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   geom_vline(xintercept = intercepts, linetype = "dotted") +
#   scale_shape_manual(values = c(15,18)) +
#   scale_color_manual(values = c('#00CCCC', '#FF8000')) +
#   labs(title = bquote(italic(.(gene))), x = "Fold-change") +
#   theme.grid
# g=ggplotGrob(plot)
# g$heights[7] = grobHeight(1)
# grid.arrange(g, table, ncol = 2)
# 
# pdf(file = paste0("forestplot_", gene, ".pdf"), 12, 6)
# leftPanel
# dev.off()

####################################
# Forest plot expression change, one region pair

tab <- data.frame(summaryDiffExpr[geneId, "R1-R6", , ])
tab$donors <- gsub("donor", "Donor ", rownames(tab))
tab$donors <- factor(tab$donors, levels = unique(rev(tab$donors)))
# Side table info
tab$rmd <- paste0(round(tab$Estimate, digits = 2), " (", round(tab$lower95, digits = 2), ", ", round(tab$upper95, digits = 2), ")")
tab$pvalue <- signif(tab$pvalue, digits = 2)
tab$is.summary <- tab$donors == "Summary"
tab$weight <- round(tab$weight, digits = 2)

# Forest plot
fp1 <- ggplot(data = tab, aes(Estimate, donors)) +
  geom_point(aes(size = weight, shape = is.summary, color = is.summary)) +
  geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  # geom_vline(xintercept = -1, linetype = "dotted") +
  scale_shape_manual(values = c(15,18)) +
  scale_color_manual(values = c('#00CCCC', '#FF8000')) +
  labs(title = bquote(italic(.(gene))), x = "FC") +
  scale_y_discrete(labels = rev(tab$donor)) +
  # scale_x_continuous(limits = c(x.min, x.max)) +
  geom_text(aes(x =  3, label = rmd), size = 3) +
  geom_text(aes(x = 4, label = weight), size = 3) +
  # geom_text(data = colLabels, aes(x, y, label = label)) +
  theme.grid + theme(plot.margin = unit(c(0,4,0,4), "lines"))

##########################################

# Forrest plot Braak correlations
tab <- summaryLabelCor[[geneId]]
tab$donors <- factor(tab$donors, levels = unique(rev(tab$donors)))
# Side table info
tab$rci <- paste0(round(tab$r, digits = 2), " (", round(tab$lower95, digits = 2), ", ", round(tab$upper95, digits = 2), ")")
tab$is.summary <- tab$donors == "Summary"

# Forest plot
fp2 <- ggplot(data = tab, aes(r, donors)) +
  geom_point(aes(size = weight, shape = is.summary, color = is.summary)) +
  geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
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
pdf(file = paste0("forestplot_", gene, "_2.pdf"), 8,2)
fp1
fp2
dev.off()
