# Gene co-expression in whole Braak region
setwd("/tudelft.net/staff-bulk/ewi/insy/DBL/Arlin/pd_braak")
load("../ABA_Rdata/BrainExpr.RData")

# Brain donor names
donorNames <- c("donor9861", "donor10021", "donor12876", "donor14380", "donor15496", "donor15697")
names(donorNames) <- donorNames

# Labels for whole Braak region
load("resources/braakInfo.RData") # Braak stage label vectors
wholeBraak <- lapply(braakLabels, function(d){ d != 0 })# Merge labels 1-6
sapply(wholeBraak, sum)

# Average co-expression for whole region
gene_coexpr <- lapply(donorNames, function(d){
  print(d)
  e <- brainExpr[[d]][, wholeBraak[[d]]]
  cor(t(e), method = "pearson") #r-matrix
})
gene_coexpr <- simplify2array(gene_coexpr) # genes x genes x donors
avgCoexpr <- apply(gene_coexpr, 1:2, mean)
saveRDS(avgCoexpr, file = "resources/avgCoexpr_wholeBraak.rds")