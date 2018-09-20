setwd("C:/Users/dkeo/surfdrive/pd_braak")
library(ggplot2)

tab <- read.table("vumc results/ttest_regions.txt", sep = "\t")
cols <- seq(1, ncol(tab), 3)
tab <- tab[, cols]

apply(tab, 2, function(x) x < )


tab <- tab < 0.05
as.data.frame(apply(tab2, 2, sum))
