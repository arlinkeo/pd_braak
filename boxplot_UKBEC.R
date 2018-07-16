# Validate differential expression in UK Brain bank data

setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(ggplot2)

# Load datasets
load("../UKBEC/regionExpr.RData")
load("../UKBEC/probe2gene.map.RData")

##############################################################################################
# All genes
probeNames <- probe2gene.map$exprID

# Expression in Braak regions 1, 3-6
braakRegions <- c("MEDU", "SNIG", "OCTX", "TCTX", "FCTX")
braakExpr <- regionExpr[braakRegions]

# Entrez ID in AHBA of Braak genes to exprID's in UKBEC
load("resources/braakGenes.RData")
braakGenes <- lapply(braakGenes, function(x){
  rows <- match(entrezId2Name(x), probe2gene.map$AHBA)
  rows <- rows[!is.na(rows)]
  probe2gene.map[rows, "exprID"]
})
bgAll <- unlist(braakGenes)
# overlap <- intersect(allBG, degSymbols)
bgPos <- braakGenes$positive_r
bgNeg <- braakGenes$negative_r

# Mean expression braak genes across donors
meanExprll <- lapply(braakRegions, function(b) {
  expr <- braakExpr[[b]]
  geneMean <- apply(expr, 1, function(x) mean(x, na.rm = TRUE)) # mean of genes
  meanExpr <- c(geneMean[bgPos], geneMean[bgNeg])
  corr <- c(rep("r>0", length(bgPos)), rep("r<0", length(bgNeg)))
  region <- rep(b, length(meanExpr))
  data.frame(meanExpr, corr, region)
})
df <- Reduce(rbind, meanExprll)
df$region <- factor(df$region, levels = braakRegions)

# Boxplot braak genes expression in UK brainBank
theme <- theme(panel.background = element_blank(), panel.grid = element_blank(), 
               axis.line = element_line(colour = "black"),
               legend.title = element_blank())
p1 <- ggplot(df) + 
  geom_boxplot(aes(x = factor(corr), y = meanExpr, fill = factor(region, levels = braakRegions))) +
  labs(x = "Correlation", y = "Expression") +
  ggtitle("Expression of Braak-related genes in UKBEC") +
  theme

p2 <- ggplot(df) + 
  geom_boxplot(aes(y = meanExpr, x = region, fill = region)) +
  labs(x = "Braak stage region", y = "Expression") +
  ggtitle("Expression of Braak-related genes in UKBEC") +
  scale_x_discrete(expand=c(0.2,0), labels = c(1, 3:6)) +
  facet_grid(.~corr, scales = 'free', space = 'free', switch = "y") +
  theme

pdf("boxplot_UKBEC.pdf", 8,6)
print(p1)
print(p2)
dev.off()
