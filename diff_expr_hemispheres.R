# Differential expression of Left & right hemisphere
setwd("C:/Users/dkeo/surfdrive/pd_braak")
library(metafor)
source("PD/base_script.R")
load("../ABA_Rdata/BrainExprNorm.RData")
load("resources/roiSamples.RData")

ontology <- read.csv("../ABA_human_processed/Ontology_edited.csv")

# roi <- c("myelencephalon", "substantia nigra") 
roi <- c("myelencephalon", "pontine tegmentum", "substantia nigra", "CA2 field", 
         "occipito-temporal gyrus", "cingulate gyrus", "temporal lobe", 
         "frontal lobe", "parietal lobe") 

# select probe IDs for left and right hemisphere of roi
selectIds <- function(r){ # for a single structure
  row <- match(r, ontology$name)
  id <- ontology$id[row]
  rows <- grep(id, ontology$structure_id_path)
  hemisphere <- ontology$hemisphere[rows]
  rows_left <- rows[hemisphere == "L"]
  rows_right <- rows[hemisphere == "R"]
  id_left <- ontology$id[rows_left]
  id_right <- ontology$id[rows_right]
  list(left = id_left, right = id_right)
}

regionIDs <- sapply(roi, selectIds, simplify = FALSE)

# Select columns/samples per donor based on probe IDs
roiSamples <- lapply(regionIDs, function(r){
  lapply(donorNames, function(d){
    expr <- brainExprNorm[[d]]
    colnames <- colnames(expr)
    samples <- sapply(r, function(ids){
      ids <- intersect(ids, colnames)
      cols <- colnames %in% ids
      cols #as.integer(cols)
    })
    rownames(samples) <- colnames
    samples
  })
})
sizes <- lapply(roiSamples, function(r){
  sapply(r, function(d){
    apply(d, 2, sum)
  })
})

# T-test function
ttestGene <- function(a, b) {
  test2tail <- t.test(a, b) # two-sided
  estimate <- test2tail$estimate
  names(estimate) <- NULL
  confidence95 <- test2tail$conf.int
  c('pvalue' = test2tail$p.value, 
    'meanA' = estimate[1], 'varA' = var(unlist(a)), 
    'meanB' = estimate[2], 'varB' = var(unlist(b)), 
    'lower95' = confidence95[1], 'upper95' = confidence95[2],
    'sizeA' = length(a), 'sizeB' = length(b))
}

genes <- ahba.genes()
load("resources/braakGenes.RData")

# T-test for two donors
ttest <- lapply(roiSamples, function(roi){
  lapply(donorNames[1:2], function(d){
    print(d)
    samples <- roi[[d]]
    expr <- brainExprNorm[[d]]
    l <- expr[, samples[, "left"]]
    r <- expr[, samples[, "right"]]
    # lapply(braakGenes, function(cor){ # For progression genes
      t <- as.data.frame(t(sapply(genes, function(g){
        ttestGene(unlist(l[g, ]), unlist(r[g,  ]))
      })))
      t$BH <- p.adjust(t$pvalue)
      t$meanDiff <- t$meanA - t$meanB
      t$logp <- -log10(t$BH)
      t
    # })
  })
})

# # Prepare data for meta-analysis: Tables with effect sizes per donor for each gene
# geneTab <- lapply(ttest, function(r){
#   sapply(genes, function(g){
#     as.data.frame(t(sapply(r, function(d){
#       unlist(d[g, ])
#     })))
#   }, simplify = FALSE)
# })
# 
# summary <- lapply(geneTab, function(r){
#   lapply(r, function(gene){
#     
#     # Get variance and confidence intervals
#     m1i <- gene$meanA
#     m2i <- gene$meanB
#     n1i <- gene$sizeA
#     n2i <- gene$sizeB
#     sd1i <- sqrt(gene$varA)
#     sd2i <- sqrt(gene$varB)
#     # Get mean difference and its variance
#     t <- escalc(measure = "MD", m1i = m1i, m2i = m2i, n1i = n1i, n2i = n2i, sd1i = sd1i, sd2i = sd2i)
#     t <- summary(t)
#     # Get summary estimate
#     summary <- rma(t$yi, t$vi, method = "DL", test = "t") # Summary effect size
#     
#     donors <- sapply(rownames(gene), function(n){ gsub("donor", "Donor ", n)})
#     meanDiff <- as.numeric(t$yi)
#     varDiff <- t$vi
#     lower95 <- t$ci.lb
#     upper95 <- t$ci.ub
#     weight <- round(weights(summary), digits = 2)
#     pvalue <- gene$pvalue
#     t <- data.frame(donors, meanDiff, varDiff, lower95, upper95, pvalue, weight)
#     
#     # Combine into table
#     t <- rbind(t, 'summary' = list("Summary", summary$beta, summary$se^2 , summary$ci.lb, summary$ci.ub, 
#                                    summary$pval, sum(weight)))
#     t
#   })
# })
# 
# diffExpr <- lapply(summary, function(roi){
#   t <- do.call(rbind.data.frame, lapply(roi, function(g) g["summary",]))
#   t$benjamini_hochberg <- p.adjust(t$pvalue)
#   print(sum(t$pvalue < 0.05))
#   t
# })

#Volcano plots

library(ggplot2)
library(ggpubr)
# library(gridExtra)

theme <- theme(legend.position = "none",
               panel.background = element_blank(),
               axis.line = element_line(colour = "black"),
               axis.title =  element_text(size = 12),
               plot.title = element_text(size = 12, face = "bold"))

# # axis limits for all plots
# ctab <- Reduce(rbind, diffExpr)
# xmax <- max(ctab$meanDiff)
# xmin <- min(ctab$meanDiff)
# ymax <- ceiling(max(ctab$'logp'[is.finite(tab$'logp')]))
braakGenes <- unlist(braakGenes)

# volcano plot for every region and donor
pdf("diff_expr_hemispheres.pdf", 8, 4)
plotll <- sapply(names(ttest), function(r){
  
  p <- lapply(donorNames[1:2], function(d){
    tab <- ttest[[r]][[d]]
    tab$info <- as.numeric(rownames(tab) %in% braakGenes)
    tab$info <- as.factor(tab$info) # Plotting order of data points 
    order <- order(tab$info)
    tab <- tab[order, ]
    
    ggplot(tab, aes(meanDiff, logp, colour = info)) +
      geom_point(size = 0.25, alpha = 0.3) +
      scale_colour_manual(values = c("0"="#F8766D", "1"="#00BFC4", "2"="black")) +
      labs(x = "fold-change", y = "-log10 p-value") +
      ggtitle(paste0(r, ", ", d)) +
      theme
  })
  p <- ggarrange(p[[1]], p[[2]])
  print(p)
}, simplify = FALSE)
dev.off()
