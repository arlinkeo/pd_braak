# Plot reference signal

plot.ref.signal <- function(data){
  colnames(data) <- make.names(colnames(data), unique = TRUE)
  
  df <- melt(as.matrix(expr))
  colnames(df) <- c("celltype", "sample", "expr")
  df$sample <- factor(df$sample, levels = unique(df$sample))
  
  intercepts <- match(c(1:6), braakLabels[[d]][idx])[-1]
  
  ggplot(df, aes(x=sample, y=expr, color=celltype)) + 
    geom_point() +
    geom_smooth(aes(group=celltype)) +
    geom_vline(xintercept = intercepts) +
    ggtitle(d) +
    theme(panel.background = element_blank(), 
          panel.grid = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
}
