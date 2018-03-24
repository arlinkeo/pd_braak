# Braak gene enrichment of modules
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(ggplot2)
load("resources/modules.RData")
load("resources/braakGenes.RData")
load("resources/summaryLabelCorrEG.RData")

braakGenes <- unlist(braakGenes, use.names = FALSE)
labelCor <- do.call(rbind.data.frame, lapply(summaryLabelCorrEG, function(g) g["summary",]))
labelCor$pvalue <- p.adjust(labelCor$pvalue, method = "BH")
total <- length(ahba.genes())

m <- modules[["braak1-6"]]
module_size <- sapply(m, length)

# Initialize table
tab <- data.frame(module = paste0("M", names(m)), module_size)

# Add eigengene-braak correlations
tab$eigengene_r <- labelCor$r
tab$eigengene_pvalue <- labelCor$pvalue

# Add numbers of significant GO terms
correctedTerms <- sapply(names(m), function(l){
  file <- paste0("Functional_analyses/braak1-6_modules/", l, "_goterms.txt")
  terms <- read.csv(file, header = TRUE, sep = "\t", colClasses = "character")
  rows <- which(as.numeric(terms$Benjamini) < 0.05)
  terms <- terms[rows, c("Category", "Term", "Bonferroni", "Benjamini")]
  terms
}, simplify = FALSE)
tab$number_goterms <- sapply(correctedTerms, nrow)

# Overlap with Braak genes and cell types
hyper.test <- function(a, b){
  genes <- intersect(a, b)
  overlap <- length(genes)
  ns1 <- length(a)
  ns2 <- length(b)
  p <- phyper(overlap - 1, ns1, total - ns1, ns2, lower.tail = FALSE)
  c(overlap = overlap, pvalue = p)
}


# Cell-type enrichment numbers
geneLists <- append(list(progression = braakGenes), celltype_genes)
overlap <- lapply(names(geneLists), function(n){
  type_genes <- geneLists[[n]]
  df <- as.data.frame(t(sapply(m, function(mod_genes){
    hyper.test(mod_genes, type_genes)
  })))
  df$pvalue <- p.adjust(df$pvalue, method = "BH") # corrected
  colnames(df) <- sapply(c("Overlap_", "pvalue_"), function(x) paste0(x, n))
  df
})
overlap <- Reduce(cbind, overlap)
overlap <- overlap[, order(!c(1:ncol(overlap))%%2)]
tab <- cbind(tab, overlap)

# Add lists of PD-implicated genes (and cell-types)
geneLists <- Reduce(append, list(pdGenesID, list(progression = braakGenes), celltype_genes))
pdOverlap <-  sapply(geneLists, function(pd){
  sapply(m, function(genes){
    res <- intersect(genes, pd)
    paste0(entrezId2Name(res), collapse = ",")
  })
})
tab <- cbind(tab, pdOverlap)

# Order genes based on summary correlation with Braak labels
orderEG <- rev(order(labelCor$r))
tab <- tab[orderEG, ]

# Print results in text-file
write.table(tab, file = paste0("module_braakgene_overlap_", b, ".txt"), sep ="\t", quote = FALSE, row.names = FALSE)

# # Plot barplots 
# pdf("celltypes_enrichment.pdf", 8, 9)
# layout(matrix(c(1:5)),  5, 1)
# par(oma = c(0,0,0,0), mai = c(0.2,0.5,0.2,0))
# lapply(names(celltype_genes), function(ct){
#   name <- paste0("Overlap_", ct)
#   x <- tab[, name]
#   names(x) <- tab$module
#   barplot(x, axisnames = TRUE)
#   title(ct)
# })
# dev.off()

###############################################################
# Plot barplots 

star <- function(v){
  sapply(v, function(x) if (x<0.001) "***" else if (x<0.01) "**" else if (x<0.05) "*" else "")
}

theme <- theme(
  axis.text.x = element_text(angle = 90, vjust = 0.4, size = 7), 
  axis.ticks.x = (element_blank()),
  axis.line.y = element_line(),
  panel.background = element_blank()
)

colPal <- c("darkblue", "white", "darkred")
rampcols <- colorRampPalette(colors = colPal, space="Lab")(201)

# Cell-type plot
df <- lapply(names(celltype_genes), function(ct){
  t <- cbind(tab[, c("module", paste0("Overlap_", ct), paste0("pvalue_", ct))], rep(ct, nrow(tab)))
  colnames(t) <- c("module", "overlap", "pvalue", "celltype")
  t
})
df <- Reduce(rbind, df)
df$logp <- -log10(df$pvalue)
df$star <- star(df$pvalue)
df$module <- factor(df$module, levels = unique(df$module))
df$celltype <- factor(df$celltype, levels = unique(df$celltype))
ymax <- max(df$logp)+30

p1 <- ggplot(df) + 
  geom_col(aes(x=module, y = logp), fill = "steelblue") +
  geom_text(aes(x=module, y=logp, label=star), size = 4, angle = 90, vjust = 0.75,  hjust = -0.2) +
  scale_y_continuous(limits = c(0, ymax)) +
  theme +
  facet_grid(celltype~.)

pdf("celltypes_enrichment.pdf", 16, 6)
print(p1)
dev.off()

# Braak label correlation plot
df2 <- cbind(tab[, c("module", "eigengene_r", "eigengene_pvalue")], type = rep("Braak label correlation", nrow(tab)))
df2$star <- ifelse(df2$eigengene_pvalue<0.001, "*", "")
df2$star <- star(df2$eigengene_pvalue)
df2$module <- factor(df2$module, levels = unique(df2$module))
offset <- 0.25
df2$y <- sapply(df2$eigengene_r, function(x) if (x>0) x+offset else x-offset)
# df2$celltype <- factor(df2$celltype, levels = unique(df$celltype))
df2$color <- rampcols[as.numeric(cut(df2$eigengene_r, breaks = 201))]
df2$color<- factor(df2$color, levels = unique(df2$color))

p2 <- ggplot(df2) + 
  geom_col(aes(x=module, y = eigengene_r, fill=eigengene_r)) +
  geom_text(aes(x=module, y=y, label=star), size = 4, angle = 90, vjust = 0.75) +
  scale_y_continuous(limits = c(min(df2$eigengene_r)-0.3, max(df2$eigengene_r)+0.3)) +
  scale_fill_gradientn(colours = rampcols) +
  theme + theme(legend.position="none") +
  facet_grid(type~.)

# Module size plot
df3 <- cbind(tab[, c(1:2)], type = rep("Module size", nrow(tab)))
df3$module <- factor(df3$module, levels = unique(df3$module))

p3 <- ggplot(df3) + 
  geom_col(aes(x=module, y = module_size), fill = "steelblue") +
  theme +
  facet_grid(type~.) 

# Progression plot
df4 <- cbind(tab[, c("module", "Overlap_progression", "pvalue_progression")], type = rep("Progression genes", nrow(tab)))
df4$module <- factor(df4$module, levels = unique(df4$module))
df4$logp <- -log10(df4$pvalue_progression)
df4$star <- star(df4$pvalue_progression)
ymax <- max(df4$logp)+5

p4 <- ggplot(df4) + 
  geom_col(aes(x=module, y = logp), fill = "steelblue") +
  geom_text(aes(x=module, y=logp, label=star), size = 4, angle = 90, vjust = 0.75,  hjust = -0.2) +
  scale_y_continuous(limits = c(0, ymax)) +
  theme +
  facet_grid(type~.) 
p4

pdf("eigengene_r.pdf", 16, 2)
print(p2)
print(p3)
print(p4)
dev.off()