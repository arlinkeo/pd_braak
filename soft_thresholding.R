# Soft thresholding before obtaining clusters
library("WGCNA")

avgCoexpr <- readRDS("output/avgCoexpr_wholeBraak.rds")

# Soft thresholding of co-expression network
# sft <- pickSoftThreshold.fromSimilarity(avgCoexpr) # on server
# saveRDS(sft, file = "sft.rds")
sft <- readRDS("output/sft_similarity.rds")

# Plot the results
pdf("output/picksoftthreshold.pdf", 9, 5)
powers = c(seq(1, 10, by = 1), seq(12, 20, by = 2))
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# Clustering
hierclust_tree_power8 <- avgCoexpr^8
saveRDS(hierclust_tree, file = "output/avgCoexpr_power8.rds")
hierclust_tree_power8 <- hier.clust(hierclust_tree_power8)
saveRDS(hierclust_tree_power8, file = "output/hierclust_tree_power8.rds")
hierclust_tree_power8 <- readRDS("output/hierclust_tree_power8.rds")
  
# count modules and missed genes
t(sapply(hierclust_tree_power8, function(t){
  n <- max(t$module)
  missed <- table(t$module)["0"]
  c(modules= n, missed_genes = missed)
}))

# Modules with gene names
t <- hierclust_tree_power8$average
modNames <- sort(unique(t$module)) # Unique module names, remove module 0 (gray)
modNames <- modNames[modNames!="0"]
names(modNames) <- paste0("M", modNames)
modules_power8 <- lapply(modNames, function(m){
  t$labels[t$module == m]
})

# intersect compare
t <- sapply(modules, function(m1){
  sapply(modules_power8, function(m2){
    length(intersect(m1, m2))
  })
})
colOrder <- hclust(dist(t(t)))[["order"]]
t <- t[, colOrder]
t[t == 0] <- NA

pdf("output/module_overlap_softthreshold.pdf", 18.2, 7.6)
Heatmap(t, name = "module overlap",
        col = colorRamp2(c(0, quantile(t, 0.9, na.rm = T)), c("#EEEEEE", "red")), #c(1, max(t)),
        na_col = "gray",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8),
        row_names_side = "left",
        column_names_side = "top",
        width = unit(ncol(t)*.5, "lines"),
        height = unit(nrow(t)*.5, "lines")
)
dev.off()
