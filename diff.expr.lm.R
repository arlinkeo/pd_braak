# Differential expression based on linear model

diff.expr.lm <- function(data, ref, group){
  fit <- lm(data ~ ref + group)
  summary <- summary(fit)
  names(summary) <- colnames(fit$coefficients)
  
  # List of tables with coefficients for each group
  groups <- paste0("group", fit$xlevels$group[-1])
  tab_list <- lapply(summary, function(g){
    t <- g$coefficients[groups, -3]
    rownames(t) <- fit$xlevels$group[-1]
    t
  })
  arr <- simplify2array(tab_list) # 3D-array: groups x measures x genes
  arr <- aaply(arr, 1, function(t){ # Correct P for number of genes (per group)
    bh <- p.adjust(t["Pr(>|t|)", ], method = "BH")
    rbind(t, 'BH' = bh)
  })
  # df <- melt(arr[, "Pr(>|t|)", ])
  # df$value <- p.adjust(df$value, method = "BH") # P-value corrected for genes and braak regions
  # bh <- dcast(df, Var1 ~ Var2)[, -1]
  # arr <- abind(arr, 'BH' = bh, along = 2) # Returns 3D-array: groups x measures x rows (genes)
  arr
}
