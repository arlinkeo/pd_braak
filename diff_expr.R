# Differential expression in braak regions of PD
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(metafor)
library(reshape2)
library(ggplot2)
library(plyr)
load("resources/braakStages.RData")

# Pairwise combinations of Braak regions 1-6
braakPairs <- t(combn(braakNames, 2))
rownames(braakPairs) <- apply(braakPairs, 1, paste, collapse = "-")
colnames(braakPairs) <- c("region_A", "region_B")

# Function T-test for each gene
ttestGene <- function(a, b) {
  test2tail <- t.test(a, b) # two-sided
  estimate <- test2tail$estimate
  names(estimate) <- NULL
  confidence95 <- test2tail$conf.int
  c('pvalue' = test2tail$p.value, 
    'meanA' = estimate[1], 'varA' = var(unlist(a)), 
    'meanB' = estimate[2], 'varB' = var(unlist(b)), 
    'sizeA' = length(a), 'sizeB' = length(b),
    'lower95' = confidence95[1], 'upper95' = confidence95[2])
}

# T-test to get p-values and CI's
# Donors -> Braak region pairs -> genes (table)
brainExpr <- readRDS("../ABA_Rdata/BrainExprNorm.rds")
genes <- ahba.genes()
ttest <- lapply(donorNames, function(d){
  print(d)
  labels <- braakStages[[d]]
  expr <- brainExpr[[d]]
  exprll <- apply(labels, 2, function(v){ # expr. in Braak regions 1-6
    expr[, v]
  })
  alply(braakPairs, 1, function(p){
    print(p)
    region_a <- exprll[[p[1]]]
    region_b <- exprll[[p[2]]]
    genesTab <- as.data.frame(t(sapply(genes, function(g){
      ttestGene(region_a[g, ], region_b[g, ])
    })))
    genesTab$benjamini_hochberg <- p.adjust(genesTab$'pvalue' , method = "BH")
    genesTab
  }, .dims = TRUE) # keep names
})
save(ttest, file = "resources/ttest.RData")
# load("resources/ttest.RData")

# Number of diff. expr. genes
sapply(ttest, function(d){
  sapply(d, function(rp){
    meanDiff <- rp$meanA - rp$meanB
    sum(rp$benjamini_hochberg < 0.05 & abs(meanDiff) > 1)
  })
})

# Braak region pairs -> Genes -> Donors (table)
diffExpr <- sapply(rownames(braakPairs), function(p){
  sapply(genes, function(g){
    as.data.frame(t(sapply(donorNames, function(d){
      res <- unlist(ttest[[d]][[p]][g, ])
    })))
  }, simplify = FALSE)
}, simplify = FALSE)
save(diffExpr, file = "resources/diffExpr.RData")

# Meta-analysis of mean expression difference across donors
summaryDiffExpr <- sapply(names(diffExpr), function(rp){ # For each Braak region pair
  print(rp)
  rp <- diffExpr[[rp]]
  lapply(rp, function(gene){
    
    # Get variance and confidence intervals
    m1i <- gene$meanA
    m2i <- gene$meanB
    n1i <- gene$sizeA
    n2i <- gene$sizeB
    sd1i <- sqrt(gene$varA)
    sd2i <- sqrt(gene$varB)
    # Get mean difference and its variance
    t <- escalc(measure = "MD", m1i = m1i, m2i = m2i, n1i = n1i, n2i = n2i, sd1i = sd1i, sd2i = sd2i)
    t <- summary(t)
    
    # Get summary estimate
    summary <- rma(t$yi, t$vi, method = "DL", test = "t") # Summary effect size
    
    donors <- sapply(rownames(gene), function(n){ gsub("donor", "Donor ", n)})
    meanDiff <- as.numeric(t$yi)
    varDiff <- t$vi
    lower95 <- t$ci.lb
    upper95 <- t$ci.ub
    weight <- round(weights(summary), digits = 2)
    pvalue <- gene$pvalue
    t <- data.frame(donors, meanDiff, varDiff, lower95, upper95, pvalue, weight)
    
    # Combine into table
    rbind(t, 'summary' = list("Summary", summary$beta, summary$se^2 , summary$ci.lb, summary$ci.ub, 
                                   summary$pval, sum(weight)))
  })
}, simplify = FALSE)
save(summaryDiffExpr, file = "resources/summaryDiffExpr.RData")

# Bar plot
load("resources/summaryDiffExpr.RData")
summTables <- lapply(summaryDiffExpr, function(l){
  t <- do.call(rbind.data.frame, lapply(l, function(g) g["summary",]))
  t$BH <- p.adjust(t$pvalue, method = "BH")
  t
})
degLists <- t(sapply(summTables, function(t){
  positive_r <- sum(t$meanDiff < -1 & t$BH < 0.05)
  negative_r <- -sum(t$meanDiff > 1 & t$BH < 0.05)
  c('higher' = positive_r, 'lower' = negative_r)
}))
rownames(degLists) <- sapply(rownames(degLists), function(x){
  x <- gsub("braak", "Braak ", x)
  gsub("-", " - ", x)
})
df <- melt(degLists)
colnames(df) <- c("region_pair", "dir", "deg")
df$region_pair <- factor(df$region_pair, levels = rev(unique(df$region_pair)))
df$y <- ifelse(df$dir == "higher", df$deg+700, df$deg-800)

p1 <- ggplot(df) + 
  geom_col(aes(x=region_pair, y = deg, fill=dir)) + 
  geom_text(aes(x=region_pair, y= y, label=format(abs(df$deg), big.mark=","))) + 
  scale_fill_manual(values = c("red", "blue")) +
  scale_y_continuous(expand = c(0.1,0.1)) +
  coord_flip() +
  labs(x = "", y = "Number of differentially expressed genes") +
  theme(
      axis.text = element_text(size = 11),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_blank(),
      legend.title = element_blank()
    )
pdf("diff_expr_barplot.pdf", 6, 4)
print(p1)
dev.off()