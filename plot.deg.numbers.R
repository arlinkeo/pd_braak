# Bar plot with  number of differentially expressed genes

plot.deg.numbers <- function(l){
  
  numbers <- t(sapply(l, function(b){
    t <- do.call(rbind.data.frame, lapply(b, function(g) g["summary",]))
    t$BH <- p.adjust(t$pvalue, method = "BH")
    t
    negative_r <- -sum(t$estimate < -1 & t$BH < 0.05)
    positive_r <- sum(t$estimate > 1 & t$BH < 0.05)
    c('negative' = negative_r, 'positive' = positive_r)
    
  }))
  rownames(numbers) <- gsub("-", " - ", rownames(numbers))
  
  df <- melt(numbers)
  colnames(df) <- c("region", "dir", "n")
  df$region <- factor(df$region, levels = rev(unique(df$region)))
  df$y <- ifelse(df$dir == "positive", df$n+500, df$n-500)
  
  ggplot(df) + 
    geom_col(aes(x=region, y = n, fill=dir), size = 0.5, colour = "black") + 
    geom_text(aes(x=region, y= y, label=format(abs(df$n), big.mark=","))) + 
    scale_fill_manual(values = c("blue", "red")) +
    scale_y_continuous(expand = c(0.1,0.1)) +
    coord_flip() +
    labs(x = "", y = "Number of differentially expressed genes") +
    theme(
      axis.text = element_text(size = 11),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_blank(),
      legend.title = element_blank()
    )
}
