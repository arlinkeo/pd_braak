# Differential expression of Left & right hemisphere
setwd("C:/Users/dkeo/surfdrive/pd_braak")
# library(metafor)
library(ggplot2)
library(ggpubr)
source("PD/base_script.R")
source("PD/sample.ids.R")
source("PD/t.test.table.R")
load("../ABA_Rdata/BrainExpr.RData")

roi <- c("myelencephalon", "pontine tegmentum", "substantia nigra", "CA2 field", 
         "occipito-temporal gyrus", "cingulate gyrus", "temporal lobe", 
         "frontal lobe", "parietal lobe")
regionIDs <- sample.ids(roi, hemisphere = TRUE)

# Select columns/samples per donor based on AHBA probe IDs
roiSamples <- lapply(regionIDs, function(r){
  lapply(donorNames[1:2], function(d){ # Only for 2 donors with samples from both hemispheres
    expr <- brainExpr[[d]]
    colnames <- colnames(expr)
    samples <- sapply(r, function(ids){
      ids <- intersect(ids, colnames)
      cols <- which(colnames %in% ids) # columnn indices
      names(cols) <- colnames[cols]
      cols # column indices
    })
  })
})
lapply(roiSamples, function(r){
  sapply(r, function(d){
    sapply(d, length)
  })
})

# T-test for two donors
ttest <- lapply(roiSamples, function(roi){
  tab <- lapply(donorNames[1:2], function(d){
    print(d)
    samples <- roi[[d]] # indices L & R
    expr <- brainExpr[[d]]
    l <- expr[, samples$left]
    r <- expr[, samples$right]
    t.test.table(l, r)
  })
  tab <-simplify2array(tab) # 3D array: genes x measures x donors
})
ttest <- simplify2array(ttest) # 4D array: genes x measures x donors x regions

# Number of differentially expresssed genes
t(apply(ttest, c(3,4), function(x){
  meanDiff <- x[, "meanA"] - x[, "meanB"]
  sum(x[, "BH"] < 0.05)
}))

#Volcano plots
theme <- theme(legend.position = "none",
               panel.background = element_blank(),
               axis.line = element_line(colour = "black"),
               axis.title =  element_text(size = 12),
               plot.title = element_text(size = 12, face = "bold"))

# braakGenes <- unlist(braakGenes)

# volcano plot for every region and donor
pdf("diff_expr_hemispheres.pdf", 8, 4)
plotll <- lapply(dimnames(ttest)[[4]], function(r){
  p <- lapply(donorNames[1:2], function(d){
    tab <- as.data.frame(ttest[, , d, r])
    tab$meanDiff <- tab$meanB-tab$meanA
    tab$logp <- -log10(tab$BH)
    
    ggplot(tab, aes(meanDiff, logp)) +
      geom_point(size = 0.25, alpha = 0.3) +
      labs(x = "fold-change", y = "-log10 p-value") +
      ggtitle(paste0(r, ", ", d)) +
      theme
  })
  p <- ggarrange(p[[1]], p[[2]])
  print(p)
})
dev.off()
