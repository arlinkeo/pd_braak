# Bar plot with  number of differentially expressed genes

plot.deg.numbers <- function(t){
  df <- melt(t)
  colnames(df) <- c("region", "dir", "size")
  df$size[df$dir == "downregulated"] <- df$size[df$dir == "downregulated"]*-1
  df$region <- factor(df$region, levels = rev(unique(df$region)))
  lab_offset <- max(df$size)*0.5
  df$y <- ifelse(df$dir == "upregulated", df$size + lab_offset, df$size - lab_offset)
  
  ggplot(df) +
    geom_col(aes(x=region, y = size, fill=dir), size = 0.5, colour = "black") +
    geom_text(aes(x=region, y= y, label=format(abs(df$size), big.mark=","))) +
    scale_fill_manual(values = c("blue", "red")) +
    scale_y_continuous(expand = c(0.1,0.1)) +
    coord_flip() +
    labs(x = "", y = "Number of differentially expressed genes") +
    theme(
      axis.text = element_text(size = 11),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_blank(),
      legend.title = element_blank(),
      legend.position = "top"
    )
}
