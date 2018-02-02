# Validate differential expression in UK Brain bank data

setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(ggplot2)

# Load datasets
load("../UKBEC_RData/regionExpr.RData")
load("../UKBEC/expr.maps2.RData")

##############################################################################################
# All genes
probeNames <- rownames(regionExpr$CRBL)

# Expression in Braak regions 1, 3-6
braakRegions <- c("MEDU", "SNIG", "OCTX", "TCTX", "FCTX")
braakExpr <- regionExpr[braakRegions]

# Gene probe mapping functions
geneProbeMap <- expr.map2[probeNames, ]
tID2expxrID <- function(x){geneProbeMap$exprID[match(x, geneProbeMap$tID)]}
AHBA2exprID<- function(x){geneProbeMap$exprID[match(x, geneProbeMap$AHBA)]}
entrezId2exprID <- function(x) AHBA2exprID(entrezId2Name(x))
# exprID2AHBA <- function(x){geneProbeMap$AHBA[match(x, geneProbeMap$exprID)]}

# Entrez ID in AHBA of Braak genes to exprID's in UKBEC
load("resources/braakGenes.RData")
braakGenes <- lapply(braakGenes, entrezId2exprID)
braakGenes <- lapply(braakGenes, na.exclude)
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