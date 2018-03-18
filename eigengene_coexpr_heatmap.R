#Eigen gene co-expression
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
load("resources/eigenExpr.RData")
eigenExpr <- eigenExpr$`braak1-6`
load("resources/summaryLabelCorrEG.RData")
labelCor <- do.call(rbind.data.frame, lapply(summaryLabelCorrEG, function(g) g["summary",]))

# Co-expression for each donor
coexpr <- lapply(eigenExpr, function(x){
  cor(t(x))
})
coexpr <- apply(simplify2array(coexpr), 1:2, mean)

colPal <- c("darkblue", "white", "darkred")
rampcols <- colorRampPalette(colors = colPal, space="Lab")(200)

rowColor <- rampcols[as.numeric(cut(labelCor$r, breaks = 201))]

pdf("eigengene_coexpr_heatmap.pdf", 8, 8)
heatmap.2(coexpr, col = rampcols, RowSideColors = rowColor,
          cexCol = 0.5, cexRow = 0.5,# revC = TRUE,
          scale = "none", trace = "none", key = TRUE, 
          main = paste0("Co-expression of module eigen genes in Braak 1-6")
)
dev.off()
