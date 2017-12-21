# Average co-expression
setwd("/tudelft.net/staff-bulk/ewi/insy/DBL/Arlin/pd_braak")

source("PD/base_script.R")
load("resources/gene_coexpr.Rdata") # list with corr matries for each donor

gene_coexpr <- lapply(gene_coexpr, function(x)x$r)
avgCor <- apply(simplify2array(gene_coexpr), 1:2, mean)
save(avgCor, file = "resources/avgCor.RData")