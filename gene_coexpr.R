# Gene co-expression in whole Braak region

wholeBraak <- lapply(braakLabels, function(d){ d != 0 })# Merge labels 1-6
sapply(wholeBraak, sum)

# Average co-expression for whole region
gene_coexpr <- lapply(donorNames, function(d){
  print(d)
  e <- brainExpr[[d]][, wholeBraak[[d]]]
  cor(t(e), method = "pearson") #r-matrix
})
gene_coexpr <- simplify2array(gene_coexpr) # genes x genes x donors
saveRDS(gene_coexpr, file = "output/gene_coexpr.rds")
avgCoexpr <- apply(gene_coexpr, 1:2, mean)
saveRDS(avgCoexpr, file = "output/avgCoexpr_wholeBraak.rds")
