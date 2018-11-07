#  Differential expression based on linear regression
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(metafor)
library(reshape2)
library(ggplot2)
library(abind)
library(plyr)
library("RDAVIDWebService")
load("resources/braakInfo.RData")
load("../ABA_Rdata/BrainExpr.RData")

# Load function
source("PD/diff.expr.lm.R")
source("PD/plot.deg.numbers.R")

# Cell-type genes
celltypes <- sapply(c("Neurons", "Astrocytes", "Oligodendrocytes", "Microglia", "Endothelial_cells"), function(type){
  file = paste0("brainscope_celltypes/", type, ".txt")
  as.character(read.csv(file, header = TRUE)$entrez_id)
}, simplify = FALSE)

# Cell-type eigengene expression across all samples (whole brain)
ct_ref <- lapply(donorNames, function(d){
  t(sapply(celltypes, function(ct){
    x <- t(brainExpr[[d]][ct, ])
    eg <- prcomp(x)$x[, 1]# 1st PC (eigen gene expr)
    mean <- apply(x, 1, mean)
    if (cor(eg, mean) > 0) eg else -eg # flip sign of eigen gene based on the data
  }))
})

# Fit linear model & get coefficients
diffExpr <- lapply(donorNames, function(d){
  idx <- unlist(braak_idx[[d]]) # idx for samples
  braak <- paste0("braak", braakLabels[[d]][idx]) # braak labels
  ct <- t(ct_ref[[d]][, idx]) # samples x cell-types
  expr <- t(brainExpr[[d]][, idx]) # samples x genes
  diff.expr.lm(expr, ct, braak) # returns 3D array: braak x measures x genes
})
diffExpr <- simplify2array(diffExpr) # 4D-array: braak x measures x genes x donors

apply(diffExpr, c(1,4), function(b){
  sum(b["BH", ] < 0.05 & abs(b["Estimate",]) > 1)
})

# Meta-analysis of mean expression difference across donors
summaryCoef <- aaply(diffExpr, c(1,3), function(g){ # For each Braak region and gene
  gene <- as.data.frame(t(g))
  summary <- rma(yi = gene$Estimate, vi = gene$`Std. Error`^2, method = "DL", test = "t") # Summary effect size
  gene$weight <- weights(summary)
  t <- rbind(gene, 'summary' = list(summary$beta, summary$se, summary$pval, NA, sum(gene$weight)))
  as.matrix(t)
}) # 4D-array: braak region x genes x donors x measures
saveRDS(summaryCoef, file = "resources/summaryCoef.rds")

# Filter summary estimates, and correct P-values
summaryCoef2 <- summaryCoef[,,"summary", -which(dimnames(summaryCoef)[[4]] %in% c("BH", "weight"))] # braak region x genes x measures
summaryCoef2 <- aaply(summaryCoef2, 1, function(t){ # P-value corrected for genes
  b <- p.adjust(t[, "Pr(>|t|)"], method = "BH")
  cbind(t, BH = b)
}) # groups x genes x measures
# df <- melt(summaryCoef2[,,"Pr(>|t|)"])
# df$value <- p.adjust(df$value, method = "BH") # P-value corrected for genes AND braak regions
# bh <- dcast(df, X1 ~ X2)[, -1] # Braak regions x genes
# summaryCoef2 <- abind(summaryCoef2, 'BH' = bh, along = 3) # Returns 3D-array: groups x measures x rows (genes)

diffGenes <- apply(summaryCoef2, 1, function(t){
  # rownames(t)[which(abs(t[, "Estimate"])>1 & t[,"BH"] < 0.05)]
  down <-  rownames(t)[which(t[, "Estimate"] < -1 & t[,"BH"] < 0.05)]
  up <- rownames(t)[which(t[, "Estimate"] > 1 & t[,"BH"] < 0.05)]
  list(down = down, up = up)
})
sapply(diffGenes, function(x){
  deg <- sapply(x, length)
  c(deg, sum = sum(deg))
})

degR6 <- diffGenes$braak6
braakgenes2 <- degR6 
save(braakgenes2, file = "resources/braakGenes2.RData")

#Functional enrichment of genes down- and upregulated between Braak region 1 and 6
david<-DAVIDWebService$new(email="D.L.Keo@tudelft.nl",
                           url="https://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
bg_list <- ahba.genes()
bg <- addList(david, bg_list, idType = "ENTREZ_GENE_ID", listName = "AHBAbackground", listType = "Background")
bg
t <- 0.05 # EASE p-value threshold
setTimeOut(david, 200000)

lapply(names(degR6), function(dir){
  genes <- degR6[[dir]]
  result <- addList(david, genes, idType = "ENTREZ_GENE_ID", listName = dir, listType = "Gene")
  print(result)
  setCurrentBackgroundPosition(david, 1)
  getFunctionalAnnotationChartFile(david, paste0("Functional_analyses/deg_lm_", dir, "_goterms.txt"), threshold=t, count=2L)
  getClusterReportFile(david, paste0("Functional_analyses/deg_lm_", dir, "_termclusters.txt"), type = c("Term"))
})

# Correlation across fold-changes per braak region
# For genes differentially expressed between R1 and R6
braakCor <- sapply(degR6, function(g){
  t <- summaryCoef2[,g, "Estimate"]
  apply(t, 2, function(x) cor(c(0,x), c(1:6)))# cor(c(0,x), c(1:6))
})
thresholds <- seq(0.1, 0.9, 0.1)
names(thresholds) <- thresholds
sapply(braakCor, function(r){
  sapply(thresholds, function(x) sum(abs(r)> x))
})
braakgenes <- sapply(braakCor, function(r){
  names(r)[abs(r)>0.8]
})

# line plot of fc's across regions
df <- melt(t[,unlist(braakgenes)])
colnames(df) <- c("region", "gene", "fc")
ggplot(df, aes(x=region, y=fc, group=gene, color=gene)) + geom_line() + geom_point()

# Presence of PD-implicated genes
pdg <- lapply(pdGenesID, function(l){
  g <- unlist(diffGenes$braak6)
  entrezId2Name(g[which(g %in% l)])
})
fcTab <- t[, name2EntrezId(unique(unlist(pdg)))]
colnames(fcTab) <- entrezId2Name(colnames(fcTab))
fcTab

##############################################################
# Bar plot of DEGs

numbers <- t(sapply(diffGenes, function(l){
  sapply(l, length)
}))
p <- plot.deg.numbers(numbers)
pdf("diff_expr_lm_eg_barplot.pdf", 6, 3)
print(p)
dev.off()