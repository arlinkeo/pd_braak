# Plot reference signal

plot.ref.signal <- function(data, title, vline = NULL){
  colnames(data) <- make.names(colnames(data), unique = TRUE)
  
  df <- melt(data)
  colnames(df) <- c("celltype", "sample", "expr")
  df$sample <- factor(df$sample, levels = unique(df$sample))
  
  ggplot(df, aes(x=sample, y=expr, color=celltype)) + 
    geom_point() +
    geom_smooth(aes(group=celltype)) +
    geom_vline(xintercept = vline) +
    ggtitle(title) + theme_classic() + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
}
