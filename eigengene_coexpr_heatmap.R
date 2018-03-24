#Eigen gene co-expression
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(gplots)

load("resources/eigenExpr.RData")
eigenExpr <- eigenExpr$`braak1-6`
load("resources/summaryLabelCorrEG.RData")
labelCor <- do.call(rbind.data.frame, lapply(summaryLabelCorrEG, function(g) g["summary",]))
# load("resources/hierclust_tree.RData")
# tree <- hierclust_tree$`braak1-6`$average

# Co-expression for each donor
coexpr <- lapply(eigenExpr, function(x){
  cor(t(x))
})
coexpr <- apply(simplify2array(coexpr), 1:2, mean)
rownames(coexpr) <- paste0("M", rownames(coexpr))
colnames(coexpr) <- rownames(coexpr)

colPal <- c("darkblue", "white", "darkred")
rampcols <- colorRampPalette(colors = colPal, space="Lab")(201)

colColor <- rampcols[as.numeric(cut(labelCor$r, breaks = 201))]
# 
# modColor <- unique(data.frame(as.numeric(tree$module), tree$color))
# modColor <- modColor[order(modColor$as.numeric.tree.module.), ]
# rownames(modColor) <- modColor$as.numeric.tree.module.
# modColor$as.numeric.tree.module. <- NULL
# rowColor <- modColor[rownames(modColor)!=0,]

pdf("eigengene_coexpr_heatmap.pdf", 8, 8)
heatmap.2(coexpr, col = rampcols, #RowSideColors = rowColor, 
          ColSideColors = colColor,
          cexCol = 0.4, cexRow = 0.4, revC = TRUE,
          scale = "none", trace = "none", key = TRUE, 
          main = paste0("Co-expression of module eigen genes in Braak 1-6")
)
dev.off()

