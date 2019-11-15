# Boxplot of expression in Braak regions for a single gene

########## Functions for boxplots ##########

# Default theme for boxplot
theme <- theme(panel.background = element_blank(), panel.grid = element_blank(), 
               axis.line = element_line(colour = "black"),
               legend.title = element_blank(),
               legend.key = element_blank())

# Boxplot function
box.plot <- function(df, title){
  p <- ggplot(df) + 
    geom_boxplot(aes(x = label, y = expr, alpha = donor, fill = label), outlier.size = 1) +
    labs(x = "Brain region", y = bquote("Expression ("*log[2]*"-transformed)")) +
    guides(alpha=guide_legend(override.aes=list(fill=hcl(c(15,195),100,0), colour=NA))) +
    scale_alpha_discrete(labels = gsub("donor", "Donor ", donorNames)) +
    scale_fill_manual(values = unname(braakColors), guide = FALSE) +
    ggtitle(title) +
    theme
  p
}

# Prepare data (mean across genes)
prepare.data <- function(g){
  df <- lapply(donorNames, function(d) {
    expr <- brainExpr[[d]][g, ]
    expr <- apply(expr, 2, mean) # Mean across genes
    label <- braakLabels[[d]]
    # donor <- rep(d, length(expr))
    data.frame(expr, label)#, donor)
  })
  df <- melt(df, value.name = "expr")
  colnames(df)[4] <- "donor"
  df <- df[df$label != "0", ]#Remove Braak 0
  df$label <- paste0("R", df$label)
  df$label <- factor(df$label, levels = sort(unique(df$label)))
  df$donor <- factor(df$donor, levels = unique(df$donor))
  df
}

# Boxplot for a gene
boxplot.gene <- function(g, title){
  df <- prepare.data(g)
  box.plot(df, title)
}

plot.pdf <- function(name, genes){
  pdf(name, 6, 4)
  lapply(genes, function(g){
    print(entrezId2Name(g))
    r <- format(summaryLabelCor[[g]]["summary", c("r", "pvalue")], digits = 2)
    title <- bquote(italic(.(entrezId2Name(g)))*", "*italic(r)*"="*.(r$r))
    df <- prepare.data(g)
    p <- box.plot(df, title)
    print(p)
  })
  dev.off()
}

###################################################

# Boxplot for mean expression of -ve and +ve Braak genes

meanExpr <- lapply(bg, prepare.data) 
df <- melt(meanExpr)
colnames(df) <- c("label", "variable", "donor", "expr", "dir")
y_max <- max(sapply(meanExpr, function(x) max(x$expr)))
y_min <- min(sapply(meanExpr, function(x) min(x$expr)))

pdf("output/boxplot_AHBA.pdf", 6, 4)
box.plot(df, "Mean BRGs") + facet_grid(.~dir, space = "free", scales = "free") +
  scale_y_continuous(limits = c(y_min, y_max))
dev.off()

# # Mean expression in one donor
# meanExprD1 <- lapply(meanExpr, function(t){
#   rows <- which(t$donor %in% "donor9861")
#   t[rows,]
# })
# df <- melt(meanExprD1)
# colnames(df) <- c("label", "variable", "donor", "expr", "dir")
# pdf("boxplot_AHBA_donor9861.pdf", 4.8, 4)
#   ggplot(df) + 
#   geom_boxplot(aes(x = label, y = expr,fill = label), outlier.size = 1) +
#   labs(x = "Brain region", y = "Expression (log2-transformed)") +
#   scale_fill_manual(values = unname(braakColors), guide = FALSE) +
#   ggtitle("Mean BRGs in donor 9861") +
#   theme + 
#   facet_grid(.~dir, space = "free", scales = "free") +
#   scale_y_continuous(limits = c(y_min, y_max))
# dev.off()