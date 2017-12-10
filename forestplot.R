# Forest plot of differential expression in each brain for a single gene

setwd("C:/Users/dkeo/surfdrive/Parkinson")
library(ggplot2)
library(gridExtra)

source("PD/base_script.R")
load("resources/summaryDiffExprNonBraak.RData")

gene <- "GCH1"
geneId <- name2EntrezId(gene)

# Basic ggplot theme
theme.grid <- theme(panel.background = element_blank(), panel.grid = element_blank(), axis.ticks.y = element_blank(), axis.line.x = element_line(colour = "black"), 
                    legend.position = "none", axis.title.y = element_blank())

####################################
# Forest plot expression change

nonBraakRef <- "nonBraakA2"

#Get values for a gene to plot
braakList <- summaryMeanDiff[[nonBraakRef]][[geneId]]
braakList <- braakList[braakNamesMerged1]
braakList <- lapply(names(braakList), function(b){
  df <- braakList[[b]]
  df$braak <- gsub("braak", "Braak ", b)
  df
})

# Combine df's for each braak stage into one df
tab <- do.call(rbind, braakList)
# tab$name <- rownames(tab)
tab$donors <- factor(tab$donors, levels = unique(rev(tab$donors)))
tab$braak <- factor(tab$braak, levels = unique(tab$braak))
# Side table info
tab$rmd <- paste0(round(tab$meanDiff, digits = 2), " (", round(tab$lower95, digits = 2), ", ", round(tab$upper95, digits = 2), ")")
tab$benjamini_hochberg <- signif(tab$benjamini_hochberg, digits = 2)
tab$is.summary <- tab$donors == "Summary"

# x.max <- round(max(tab$upper95), digits = 2)
# x.min <- round(min(tab$lower95), digits = 2)

# Basic ggplot grid for braak stages
facet.braak <- facet_grid(braak~., scales = 'free', space = 'free', switch = "y")

# Forest plot
fp <- ggplot(data = tab, aes(meanDiff, donors)) 

# x.positions <- c(0.25, 0.55, 0.65) + x.max#c(1:3)
# colLabels <- data.frame(x = x.positions, y = x.positions, label = c("Raw mean difference", "N", "Weight"))

leftPanel <- fp +
  geom_point(aes(size = size, shape = is.summary, color = is.summary)) +
  geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dotted") +
  geom_vline(xintercept = -1, linetype = "dotted") +
  scale_shape_manual(values = c(15,18)) +
  scale_color_manual(values = c('#00CCCC', '#FF8000')) +
  labs(title = bquote(italic(.(gene))), x = "Raw mean difference") +
  scale_y_discrete(labels = rev(tab$donor)) +
  # scale_x_continuous(limits = c(x.min, x.max)) +
  geom_text(aes(x = 3.1, label = rmd), size = 3) +
  geom_text(aes(x = 3.9, label = size), size = 3) +
  geom_text(aes(x = 4.4, label = weight), size = 3) +
  geom_text(aes(x = 5, label = benjamini_hochberg), size = 3) +
  # geom_text(data = colLabels, aes(x, y, label = label)) +
  theme.grid + theme(plot.margin = unit(c(0,4,0,4), "lines")) +
  facet.braak
pdf(file = paste0("forestplot_", gene, ".pdf"), 12, 3.5)
leftPanel
dev.off()
##########################################

# Forrest plot Braak correlations

load("resources/summaryLabelCorr.RData")

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
  geom_text(aes(x = 2.2, label = braakSize), size = 3) +
  geom_text(aes(x = 2.6, label = weight), size = 3) +
  # geom_text(data = colLabels, aes(x, y, label = label)) +
  theme.grid + theme(plot.margin = unit(c(0,4,0,4), "lines"))
pdf(file = paste0("forestplot_corr_", gene, ".pdf"), 8, 2)
leftPanel
dev.off()