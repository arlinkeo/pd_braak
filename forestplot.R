# Forest plot of differential expression in each brain for a single gene

setwd("C:/Users/dkeo/surfdrive/Parkinson")
library(ggplot2)
library(metafor)

source("PD/base_script.R")
load("resources/diffGenesBraak.RData")

load("resources/braakStages.RData")
braakNames <- names(braakStages)
names(braakNames) <- braakNames
sizes <- sapply(braakStages, function(bs){
  sapply(bs, sum)
})
sizesNB <- sapply(nonBraak, sum)

gene <- "PINK1"
geneId <- name2EntrezId(gene)

#Get values for a gene
braakList <- lapply(braakNames, function(bs){
  donorList <- diffGenesList[[bs]]
  tab <- as.data.frame(t(sapply(donorList, function(t) unlist(t[geneId, ]))))
  tab$meanDiff <- unlist(tab$meanB) - unlist(tab$meanNB)# mean difference
  donors <- sapply(rownames(tab), function(n){ paste0("Donor ", unlist(strsplit(n, split = "donor"))[2])})
  tab$donor <- donors
  tab$size <- sizes[, bs]
  tab$varDiff <- (tab$varB / tab$size) + (tab$varNB / sizesNB) # variance of mean difference
  tab <- tab[, c("lower95", "upper95", "meanDiff", "varDiff", "donor", "size")]
  sumEffect <- rma(tab$meanDiff, tab$varDiff) # Summary effect size
  tab$weight <- round(weights(sumEffect), digits = 2)
  tab <- rbind(tab, 'SummaryEff' = list(sumEffect$ci.lb, sumEffect$ci.ub, sumEffect$beta, sumEffect$se^2, "Summary", sum(tab$size), sum(tab$weight)))
  tab$isSum <- tab$donor == "Summary"
  tab$braak <- bs
  tab
})

tab <- do.call(rbind, braakList)
# tab = braakList[[1]]
tab$name <- rownames(tab)
tab$name <- factor(tab$name, levels = unique(rev(tab$name)))
tab$braak <- factor(tab$braak, levels = unique(tab$braak))
# Side table plot
tab$rmd <- paste0(round(tab$meanDiff, digits = 2), " (", round(tab$lower95, digits = 2), ", ", round(tab$upper95, digits = 2), ")")

x.max <- max(tab$upper95)

# Forest plot
fp <- ggplot(data = tab, aes(meanDiff, name)) +
  geom_point(aes(size = size, shape = isSum, color = isSum)) +
  geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_shape_manual(values = c(15,18)) +
  scale_color_manual(values = c('#00CCCC', '#FF8000')) +
  scale_x_continuous(name = "Raw mean difference") +
  scale_y_discrete(labels = rev(tab$donor)) +
  geom_text(aes(x = x.max + 0.25, label = rmd), size = 3) +
  geom_text(aes(x = x.max + 0.55, label = size), size = 3) +
  geom_text(aes(x = x.max + 0.7, label = weight), size = 3) +
  theme(panel.background = element_blank(), axis.ticks.y = element_blank(), axis.line.x = element_line(colour = "black"), 
        legend.position = "none", axis.title.y = element_blank()) +
  facet_grid(braak~., scales = 'free', space = 'free', switch = "y")
fp


