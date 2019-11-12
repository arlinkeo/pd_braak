# T-test between two regions
# Prepare table for meta-analysis

t.test.table <- function(a, b) { # Data matrices a and b
  identical_rows <- identical(rownames(a), rownames(b))
  if (identical_rows){
    rows <- rownames(a)
    table <- t(sapply(rows, function(x){
      a1 <- unlist(a[x, ])
      b1 <- unlist(b[x, ])
      test2tail <- t.test(a1, b1) # two-sided
      estimate <- unname(test2tail$estimate)
      confidence95 <- test2tail$conf.int
      if (var(a1, na.rm = TRUE) != 0 | var(b1, na.rm = TRUE) != 0) {
        c('meanDiff' = estimate[2] - estimate[1], 
          'FC' = log2(estimate[2] / estimate[1]),
          'meanA' = estimate[1], 'varA' = var(a1),
          'meanB' = estimate[2], 'varB' = var(b1),
          'sizeA' = length(a1), 'sizeB' = length(b1),
          'lower95' = confidence95[1], 'upper95' = confidence95[2],
          'pvalue' = test2tail$p.value)
      } else {
        stop(paste("Variance of", x, "is 0"))
      }
    }))
    cbind(table, "BH" = p.adjust(table[, "pvalue"], method = "BH"))
  } else {
    stop("non-identical rows")
  }
}