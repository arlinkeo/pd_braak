# Validate differential expression in UK Brain bank data

setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
# library("reshape2")
library(WGCNA)
library(ggplot2)
load("../UKBEC/expr.maps.rda",verbose=T)

# Load dataets
regions <- c("CRBL", "FCTX", "HIPP","MEDU", "OCTX", "PUTM", "SNIG", "TCTX", "THAL", "WHMT")
regionExpr <- sapply(regions, function(r){
  fName <- paste("../UKBEC/expr_", r, ".txt", sep = "")
  read.csv(fName, header = TRUE, sep = " ", row.names = 1)
}, simplify = FALSE)

#Concatenate data matrices
exprConcat <- Reduce(cbind, regionExpr)

# Probes to gene based on concatenated data
probeSelection <- collapseRows(exprConcat, expr.map$tID, rownames(exprConcat), method = "maxRowVariance", connectivityBasedCollapsing = TRUE)

# Select probes for each data matrix
probeNames <- names(probeSelection$selectedRow)[probeSelection$selectedRow] # exprID
regionExpr <- lapply(regionExpr, function(x) x[probeNames, ])
# regionExpr <- lapply(regionExpr, function(x) t(scale(t(x), center = TRUE)) ) # Normallize data

save(regionExpr, file = "../UKBEC_RData/regionExpr.RData")
# load("../UKBEC_RData/regionExpr.RData")

##############################################################################################

# Function to split character string of genes by comma
genes.split <- function(x) {unlist(strsplit(x, split = ","))}

# Add names of AHBA genes in t.map and expr.map
AHBAgenes <- probeInfo$gene_symbol # all AHBA genes
mappedGenes<- apply(t.map, 1, function(t){ # For each gene, check overlap in name symbols
  genes <- unique(c(genes.split(t[[6]]), genes.split(t[[8]])))
  overlap <- intersect(genes, AHBAgenes)
  if (length(overlap) != 1) NA else overlap
})
t.map2 <- cbind(t.map, AHBA = mappedGenes)
expr.map2 <- cbind(expr.map, AHBA = t.map2$AHBA[match(expr.map$tID, t.map2$tID)])

##############################################################################################

# Gene probe mapping functions
geneProbeMap <- expr.map2[probeNames, ]
tID2expxrID <- function(x){geneProbeMap$exprID[match(x, geneProbeMap$tID)]}
AHBA2exprID<- function(x){geneProbeMap$exprID[match(x, geneProbeMap$AHBA)]}
entrezId2exprID <- function(x) AHBA2exprID(entrezId2Name(x))

# Entrez ID in AHBA of Braak genes to exprID's in UKBEC
load("resources/braakGenes.RData")
braakGenes <- lapply(braakGenes, entrezId2exprID)
braakGenes <- lapply(braakGenes, na.exclude)
bgAll <- unlist(braakGenes)
# overlap <- intersect(allBG, degSymbols)
bgDown <- braakGenes$greater0
bgUp <- braakGenes$smaller0

# Mean expression across donors
braakRegions <- c("MEDU", "SNIG", "OCTX", "TCTX", "FCTX")
braakExpr <- regionExpr[braakRegions]
meanExprll <- lapply(braakRegions, function(b) {
  expr <- braakExpr[[b]]
  geneMean <- apply(expr, 1, function(x) mean(x, na.rm = TRUE)) # mean of genes
  meanExpr <- c(geneMean[bgDown], geneMean[bgUp])
  corr <- c(rep("r>0", length(bgUp)), rep("r<0", length(bgDown)))
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