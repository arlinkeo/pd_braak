# PSEA
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(PSEA)
library(plyr)
library(ggplot2)
load("resources/braakInfo.RData")
brainExpr <- readRDS("../AHBA_Arlin/gene_expr.RDS")
load("resources/braakGenes.RData")

expr_concat <- brainExpr$donor9861#Reduce(cbind, brainExpr)
labels_concat <- braakLabels$donor9861#unlist2(braakLabels)

# Cell-type genes
celltypes <- sapply(c("Neurons", "Astrocytes", "Oligodendrocytes", "Microglia", "Endothelial_cells"), function(type){
  file = paste0("brainscope_celltypes/", type, ".txt")
  as.character(read.csv(file, header = TRUE)$entrez_id)
}, simplify = FALSE)

# Cell-type mean expression across all samples (whole brain)
ct_ref <- t(sapply(celltypes, function(ct){
  x <- t(expr_concat[ct, ])
  mean <- apply(x, 1, mean)
}))

##### PSEA #####
# R1-R6  
idx <- which(labels_concat %in% c("1", "6"))

expr <- expr_concat[, idx]
groups <- labels_concat[idx]
groups <- as.numeric(groups == "6")
ct <- ct_ref[, idx]

# Reference signals
neurons <- ct["Neurons", ]
astrocytes <- ct["Astrocytes", ]
oligodendrocytes <- ct["Oligodendrocytes", ]
microglia <- ct["Microglia", ]
endothelial_cells <- ct["Endothelial_cells", ]

# Interaction regressors
neurons_diff <- groups * neurons
astrocytes_diff <- groups * astrocytes
oligodendrocytes_diff <- groups * oligodendrocytes
microglia_diff <- groups * microglia
endothelial_cells_diff <- groups * endothelial_cells

psea <- sapply(braakGenes$entrez_id, function(g){
  gene <- t(expr[g, ])
  
  fit <- lm(gene ~ neurons + astrocytes + oligodendrocytes + microglia + endothelial_cells, subset = which(groups==0))
  par(mfrow=c(2,3), mex=0.8)
  # crplot(fit, "neurons", newplot = FALSE)
  # crplot(fit, "astrocytes", newplot = FALSE)
  # crplot(fit, "oligodendrocytes", newplot = FALSE)
  # crplot(fit, "microglia", newplot = FALSE)
  # crplot(fit, "endothelial_cells", newplot = FALSE)
  summary <- summary(fit)
  coefficients <- summary$coefficients
  celltype_fit <- coefficients[-1, c(1,4)]
  colnames(celltype_fit) <- c("celltype_beta", "celltype_pval")
  
  fit_neurons <- lm(gene ~ neurons + neurons_diff)
  # crplot(fit_neurons, "neurons", g = "neurons_diff")
  summary_neurons <- summary(fit_neurons)
  pval_neurons <- summary_neurons$coefficients["neurons_diff", "Pr(>|t|)"]
  foldchange_neurons <- (fit_neurons$coefficients[2] + fit_neurons$coefficients[3]) / fit_neurons$coefficients[2]
  
  fit_astrocytes <- lm(gene ~ astrocytes + astrocytes_diff)
  # crplot(fit_astrocytes, "astrocytes", g = "astrocytes_diff")
  summary_astrocytes <- summary(fit_astrocytes)
  pval_astrocytes <- summary_astrocytes$coefficients["astrocytes_diff", "Pr(>|t|)"]
  foldchange_astrocytes <- (fit_astrocytes$coefficients[2] + fit_astrocytes$coefficients[3]) / fit_astrocytes$coefficients[2]
  
  fit_oligodendrocytes <- lm(gene ~ oligodendrocytes + oligodendrocytes_diff)
  # crplot(fit_oligodendrocytes, "oligodendrocytes", g = "oligodendrocytes_diff")
  summary_oligodendrocytes <- summary(fit_oligodendrocytes)
  pval_oligodendrocytes <- summary_oligodendrocytes$coefficients["oligodendrocytes_diff", "Pr(>|t|)"]
  foldchange_oligodendrocytes <- (fit_oligodendrocytes$coefficients[2] + fit_oligodendrocytes$coefficients[3]) / fit_oligodendrocytes$coefficients[2]
  
  fit_microglia <- lm(gene ~ microglia + microglia_diff)
  # crplot(fit_microglia, "microglia", g = "microglia_diff")
  summary_microglia <- summary(fit_microglia)
  pval_microglia <- summary_microglia$coefficients["microglia_diff", "Pr(>|t|)"]
  foldchange_microglia <- (fit_microglia$coefficients[2] + fit_microglia$coefficients[3]) / fit_microglia$coefficients[2]
  
  fit_endothelial_cells <- lm(gene ~ endothelial_cells + endothelial_cells_diff)
  # crplot(fit_endothelial_cells, "endothelial_cells", g = "endothelial_cells_diff")
  summary_endothelial_cells <- summary(fit_endothelial_cells)
  pval_endothelial_cells <- summary_endothelial_cells$coefficients["endothelial_cells_diff", "Pr(>|t|)"]
  foldchange_endothelial_cells <- (fit_endothelial_cells$coefficients[2] + fit_endothelial_cells$coefficients[3]) / fit_endothelial_cells$coefficients[2]
  
  fc <- c(foldchange_neurons, foldchange_astrocytes, foldchange_oligodendrocytes, foldchange_microglia, foldchange_endothelial_cells)
  pval <- c(pval_neurons, pval_astrocytes, pval_oligodendrocytes, pval_microglia, pval_endothelial_cells)
  t=cbind(celltype_fit, group_fc = fc, group_pval = pval)
}, simplify = FALSE)
psea <- simplify2array(psea)

# Celltype-specific expression
m1 <- t(psea[, "celltype_pval", ])
m1 <- apply(m1, 2, function(x)p.adjust(x, method = "BH"))
apply(m1 < 0.05, 2, sum)

# Group-dependent expression
m2 <- t(psea[, "group_pval", ])
m2 <- apply(m2, 2, function(x) p.adjust(x, method = "BH"))
apply(m2 < 0.05, 2, sum)

# Plot number of significant Betas
sum1 <- apply(m1 < 0.05, 2, sum)
sum2 <- apply(m1 < 0.05 & m2 < 0.05, 2, sum)
celltype <- paste0(toupper(substr(names(sum),1,1)), substring(names(sum), 2))
celltype <- gsub("_", " ", celltype)
df <- data.frame(celltype = celltype, no_group_diff = sum1-sum2, group_diff = sum2)
df <- melt(df)
df$celltype <- factor(df$celltype, levels = unique(df$celltype)[order(sum1)])
# df$y <- ifelse(df$variable == "group_diff", df$value-50, df$value+50)
ggplot(df) + 
  geom_col(aes(x=celltype, y=value, fill=celltype, alpha = variable)) +# stat_count() +
  geom_text(aes(x=celltype, y= value-50, label=value)) +
  coord_flip() +
  labs(x = "", y = "Number of BRGs with significant Betas") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 11),
    # axis.text.x = element_blank(),
    panel.grid = element_blank()
  ) +
  guides(fill = guide_legend(reverse = TRUE))
df

