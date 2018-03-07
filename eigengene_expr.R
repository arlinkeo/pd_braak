# Eigen gene differential expression
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")

load("../ABA_Rdata/BrainExprNorm.RData")
load("resources/modules_absCor.RData")
load("resources/braakLabels.RData")

regions <- c("braak1", "braak6", "braak1-6")

samples <- lapply(braakLabels, function(labels){
  braak1 <- labels == "1"
  braak6 <- labels == "6"
  braak1to6 <- labels != "0"
  list(braak1 = braak1, braak6 = braak6, 'braak1-6' = braak1to6)
})

# PCA first component of subselection expr. matrices
eigenExpr <- lapply(donorNames, function(d){
  expr <- brainExprNorm[[d]]
  
  sapply(regions, function(r){ # For list of modules found in different brain regions
    m <-modules[[r]]
    s <- samples[[d]][[r]]
    df <- as.data.frame(t(sapply(m, function(genes){ # For each module with genes (grouped gene rows)
      modExpr <- t(expr[genes, s]) # Expr of genes in modules across samples in one braak region
      prcomp(modExpr)$x[, 1]# 1st PC (eigen gene expr)
    })))
    colnames(df) <- names(s)[s]
    df
  }, simplify = FALSE)
})
save(eigenExpr, file = "resources/eigenExpr.RData")

# plot 

ontology <- read.csv("../ABA_human_processed/Ontology_edited.csv")
load("resources/modules_braak_absCor.RData")

load("resources/hierclust_tree_absCor.RData")

pdf("eigengene_expr.pdf", 8, 4)
lapply(regions, function(b){
  braakmods <- names(modules_braak[[b]])
  tree <- hierclust_tree[[b]][["average"]]
  color <- unique(cbind(tree$module, tree$color))
  rownames(color) <- color[,1]
  color <- color[braakmods, 2]
  
  lapply(donorNames, function(d){
    mat <- eigenExpr[[d]][[b]]
    
    labels <- colnames(mat)
    graph_order <- ontology$graph_order[match(labels, ontology$id)]
    braakorder <- braakLabels[[d]][labels]
    order <- order(braakorder, -graph_order)
    mat <- mat[braakmods, order]
  
    ahbacolor <- paste0("#", ontology$color_hex_triplet[match(labels, ontology$id)])[order]
    
    matplot(t(mat), type = "l", 
            col = color, xlab = "Braak regions", ylab = "Expression",
            xaxt = "n")
    title(paste0(b, ", ", d))
    lapply(1:length(labels), function(x){
      axis(1, at = x, col = ahbacolor[x], labels = c(""), lwd = 10, lwd.ticks = 0)
    })
    
  })
})
dev.off()


# # Function T-test for each gene
# ttestGene <- function(a, b) {
#   test2tail <- t.test(a, b) # two-sided
#   estimate <- test2tail$estimate
#   names(estimate) <- NULL
#   confidence95 <- test2tail$conf.int
#   c('pvalue' = test2tail$p.value, 
#     'meanA' = estimate[1], 'varA' = var(unlist(a)), 
#     'meanB' = estimate[2], 'varB' = var(unlist(b)), 
#     'lower95' = confidence95[1], 'upper95' = confidence95[2])
# }
# 
# # T-test eigen expression in Braak 1 vs. 6
# ttest <- lapply(donorNames, function(d){
#   labels <- braakStages[[d]]
#   modulesets <- eigenExpr[[d]]
#   lapply(modulesets, function(expr){ # For list of modules found in Braak 1 and 6 resp.
#     genesTab <- as.data.frame(t(apply(expr, 1, function(eg){ # For each eigen gene row
#       region_a <- eg[as.logical(labels[, "braak1"])]
#       region_b <- eg[as.logical(labels[, "braak6"])]
#       ttestGene(region_a, region_b)
#     })))
#     genesTab$'benjamini_hochberg' <- p.adjust(genesTab$'pvalue' , method = "BH", n = nrow(genesTab)) #corrected p
#     genesTab[order(genesTab$benjamini_hochberg),]
#   })
# })
# 
# # Saved datastructure for meta-analysis: module set -> Eigen genes -> table of Donors
# diffExpr_eigengene <- lapply(braakNames, function(b){
#   mods <- names(modules[[b]])
#   list <- sapply(mods, function(eg){
#     as.data.frame(t(sapply(donorNames, function(d){
#       unlist(ttest[[d]][[b]][eg, ])
#     })))
#   }, simplify = FALSE)
# })
# save(diffExpr_eigengene, file = "resources/diffExpr_eigengene.RData")