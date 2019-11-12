# Correlate expression  with Braak stage labels

labels <- lapply(donorNames, function(d){
  idx <- unlist(braak_idx[[d]])
  braakLabels[[d]][idx]
})

expr_braak <- lapply(donorNames, function(d){
  idx <- unlist(braak_idx[[d]])
  brainExpr[[d]][, idx] # Expression in Braak regions
})
summaryLabelCor <- summary.braak.cor(expr_braak, labels)
saveRDS(summaryLabelCor, file = "output/summaryLabelCor.rds")
load("output/summaryLabelCor.RData")
