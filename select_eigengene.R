# Select modules

setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
load("resources/eigenExpr.RData")
load("resources/summaryLabelCorrEG.RData")
load("resources/braakLabels.RData")
load("resources/hierclust_tree.RData")


# Get top 10% correlated genes
labelCor <- do.call(rbind.data.frame, lapply(summaryLabelCorrEG, function(g) g["summary",]))
order <- rev(order(labelCor$r)) # order absolute corr.
corrEG <- rownames(labelCor)[order[1:floor(nrow(labelCor)*0.1)]] # top 10% genes
plot(labelCor[order,]$r, -log10(labelCor[order,]$pvalue))
hist(labelCor[order,]$r, breaks = 30)

# plot 
ontology <- read.csv("../ABA_human_processed/Ontology_edited.csv")

pdf("eigengene_expr.pdf", 8, 4)
# lapply(regions, function(b){
b="braak1-6"
  # braakmods <- names(modules_braak[[b]])
  tree <- hierclust_tree[[b]][["average"]]
  color <- unique(cbind(tree$module, tree$color))
  rownames(color) <- color[,1]
  color <- color[corrEG, 2]
  
  lapply(donorNames, function(d){
    mat <- eigenExpr[[d]][[b]][corrEG, ]
    
    labels <- colnames(mat)
    onto_rows <- match(labels, ontology$id)
    graph_order <- ontology$graph_order[onto_rows]
    braak_order <- braakLabels[[d]][labels]
    order <- order(braak_order, -graph_order)
    mat <- mat[corrEG, order] # module selection by id?
    
    ahbacolor <- paste0("#", ontology$color_hex_triplet[match(labels, ontology$id)])[order]
    matplot(t(mat), type = "l",
            col = color, xlab = "Braak regions", ylab = "Expression",
            xaxt = "n")
    title(paste0(b, ", ", d))
    lapply(1:length(labels), function(x){
      axis(1, at = x, col = ahbacolor[x], labels = c(""), lwd = 10, lwd.ticks = 0);
    })
    
  })
# })
dev.off()
