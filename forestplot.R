# Forest plot of differential expression in each brain for a single gene

setwd("C:/Users/dkeo/surfdrive/Parkinson")
library(ggplot2)
library(gridExtra)

source("PD/base_script.R")
# load("resources/diffGenesBraak.RData")
load("resources/sumEffectSize.RData")

# load("resources/braakStages.RData")
# braakNames <- names(braakStages)
# names(braakNames) <- braakNames
# sizes <- sapply(braakStages, function(bs){
#   sapply(bs, sum)
# })
# sizesNB <- sapply(nonBraak, sum)
# 
# 
# genes <- rownames(diffGenesList$braak1$donor9861)
# names(genes) <- genes
# 
# braakNames <- braakNames[1:8]
# sumEffectSize <- lapply(sumEffectSize, function(b)b[1:8])

gene <- "GBA"
geneId <- name2EntrezId(gene)

#Get values for a gene to plot
braakList <- sumEffectSize[[geneId]]

tab <- do.call(rbind, braakList)
# tab = braakList[[1]]
tab$name <- rownames(tab)
tab$name <- factor(tab$name, levels = unique(rev(tab$name)))
tab$braak <- factor(tab$braak, levels = unique(tab$braak))
# Side table info
tab$rmd <- paste0(round(tab$meanDiff, digits = 2), " (", round(tab$lower95, digits = 2), ", ", round(tab$upper95, digits = 2), ")")

x.max <- round(max(tab$upper95), digits = 2)
x.min <- round(min(tab$lower95), digits = 2)

# Basic ggplot theme
theme.grid <- theme(panel.background = element_blank(), axis.ticks.y = element_blank(), axis.line.x = element_line(colour = "black"), 
                    legend.position = "none", axis.title.y = element_blank())
# Basic ggplot grid for braak stages
facet.braak <- facet_grid(braak~., scales = 'free', space = 'free', switch = "y")

# Forest plot
fp <- ggplot(data = tab, aes(meanDiff, name)) 

x.positions <- c(0.25, 0.55, 0.65) + x.max#c(1:3)
colLabels <- data.frame(x = x.positions, y = x.positions, label = c("Raw mean difference", "N", "Weight"))

leftPanel <- fp +
  geom_point(aes(size = size, shape = isSum, color = isSum)) +
  geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dotted") +
  geom_vline(xintercept = -1, linetype = "dotted") +
  scale_shape_manual(values = c(15,18)) +
  scale_color_manual(values = c('#00CCCC', '#FF8000')) +
  labs(title = paste0("Summary effect size of ", gene), x = "Raw mean difference") +
  scale_y_discrete(labels = rev(tab$donor)) +
  # scale_x_continuous(limits = c(x.min, x.max)) +
  geom_text(aes(x = x.max + 0.3, label = rmd), size = 3) +
  geom_text(aes(x = x.max + 0.7, label = size), size = 3) +
  geom_text(aes(x = x.max + 1.0, label = weight), size = 3) +
  # geom_text(data = colLabels, aes(x, y, label = label)) +
  theme.grid + theme(plot.margin = unit(c(0,4,0,4), "lines")) +
  facet.braak
pdf(file = paste0("forestplot_", gene, ".pdf"), 12, 8)
leftPanel
dev.off()

rightPanel <- fp +
  geom_text(aes(x = 2, label = rmd), size = 3) +
  geom_text(aes(x = 3, label = size), size = 3) +
  geom_text(aes(x = 4, label = weight), size = 3) +
  # geom_text(data = colLabels, aes(x, y, label = label), vjust = 0.5) +
  theme.grid +
  theme(axis.text.x = element_text(colour = "white"), axis.text.y = element_blank(),
        axis.line.x = element_line(colour = "white"), axis.ticks.x = element_line(colour = "white"),
        axis.title.x = element_text(colour = "white"), strip.text = element_blank()) +
  facet.braak

# rightPanel <- fp + 
#   annotate("text", x = rep(1, nrow(tab)), y = c(1:nrow(tab)), label = tab$rmd) +
#   annotate("text", x = c(1:3), y = c(nrow(tab)+1),label = c("Raw mean difference", "N", "Weight")) + # Headers
#   theme.grid
# rightPanel

p <- grid.arrange(leftPanel, rightPanel, widths = c(3,1), ncol = 2)
p