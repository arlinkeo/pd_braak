# Gene co-expression in whole Braak region

# Average co-expression for whole region
gene_coexpr <- lapply(donorNames, function(d){
  print(d)
  e <- brainExpr[[d]][, unlist(braak_idx[[d]])] # whole braak region
  cor(t(e), method = "pearson") #r-matrix
})
gene_coexpr <- simplify2array(gene_coexpr) # genes x genes x donors
saveRDS(gene_coexpr, file = "output/gene_coexpr.rds")
avgCoexpr <- apply(gene_coexpr, 1:2, mean)
saveRDS(avgCoexpr, file = "output/avgCoexpr_wholeBraak.rds")
