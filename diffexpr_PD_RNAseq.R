setwd("M:/doorgeefluik/Arlin Keo doorgeefuik/RNAseq")
options(stringsAsFactors = FALSE)
library(ggplot2)
library(ggsignif)
library(reshape2)
library(RColorBrewer)
library(PSEA)

##############################################################################
# Prepare and load data

# Braak region colors
braakColors <- brewer.pal(6, "Set2")
names(braakColors) <- paste0("R", c(1:6))

# function to select colors for each region in this dataset
braak.color <- function(x) {
  unname(unlist(lapply(x, function(r){
    if (r == "SN") braakColors[3]
    else braakColors[5]
  })))
}

# Sample info
samples <- list(control = 
                  list(
                    SN = c("A53.3", "A34.3", "Anke3", "A36.3", "A37.2", "A77.2", "A45.3", "Anke4", "Anke5", "Anke7", "Anke6"),
                    GTM = c("A54.3", "A11.3",  "A10.3", "A38.3", "A03.3",  "A15.2",  "A50.2", "A71.3", "A74.2")
                  ), 
                PD = 
                  list(
                    SN = c("A30.3", "A24.3", "A43.2", "Anke2", "Anke1", "A60.2", "A57.1", "A55.2", "A51.2", "A47.2"),
                    GTM = c("A21.3", "A05.3", "A08.3", "A01.3", "A23.3")
                  ))
sapply(samples, function(x)sapply(x, length))
sample_info <- melt(samples)
colnames(sample_info) <- c("code", "region", "disease")
rownames(sample_info) <- sample_info$code

#Load PD and control data
data_PD_region <- read.csv2("region comparison/compare_condition_to_SN Park PD Late_VS_GTM Park PDD Late.csv", skip = 3)
rownames(data_PD_region) <- data_PD_region$ID
data_CTRL_region <- read.csv2("region comparison/compare_condition_to_SN nonDem CTL_VS_GTM nonDem CTL.csv", skip = 3)
rownames(data_CTRL_region) <- data_CTRL_region$ID

# Select expression data and combine into 1 data matrix
expr <- cbind(data_PD_region[, unlist(samples$PD)], data_CTRL_region[rownames(data_PD_region), unlist(samples$control)])

# Load region data
data_SN_disease <- read.csv2("disease comparison/regroup_compare_condition_to_SN Park PD Late_VS_SN nonDem CTL.csv")
rownames(data_SN_disease) <- data_SN_disease$ID
data_GTM_disease <- read.csv2("disease comparison/compare_condition_to_GTM Park PDD Late_VS_GTM nonDem CTL.csv")
rownames(data_GTM_disease) <- data_GTM_disease$ID

# tables from excel
data_list <- list('R3_PD_vs_CTRL' = data_SN_disease, 'R4/R5_PD_vs_CTRL' = data_GTM_disease, 
                  'PD_SN_vs_GTM' = data_PD_region, 'CTRL_SN_vs_GTM' = data_CTRL_region)
# Info
# genes <- rownames(expr)
regions <- c("SN", "GTM")
names(regions) <- regions
diseases <- c("control", "PD")
names(diseases) <- diseases

##############################################################################

# Number of DEGs between diseases and between regions
degs <- lapply(data_list, function(t){
  down <- rownames(t)[which(t$log2FoldChange < -1 & t$padj<0.05)]
  up <-  rownames(t)[which(t$log2FoldChange > 1 & t$padj<0.05)]
  list(down = down, up = up)
})
no_degs <- sapply(degs, function(x) sapply(x, length))

# Heatmap number of differentially expressed genes (DESeq2)
m1 <- data.frame('R3.vs.R4/R5' = apply(no_degs[, c(3,4)], 2, sum), row.names = rev(diseases))
m2 <- data.frame('PD.vs.controls' = apply(no_degs[, c(1,2)], 2, sum), row.names = c("R3", "R4/R5"))
max_degs <- max(m1, m2)

pdf("heatmap_diffgenes.pdf", 2.3, 2)
t1 <- melt(as.matrix(m1))
t1$Var1 <- factor(t1$Var1, levels = rev(unique(t1$Var1)))
ggplot(t1) +
  geom_tile(aes(Var2, Var1, fill=value), color = "black") +
  geom_text(aes(Var2, Var1, label=value)) +
  scale_x_discrete(position = "top") +
  scale_fill_gradient(low="white", high = "chocolate", limits = c(0,max_degs)) +
  theme(axis.title=element_blank(), axis.ticks=element_blank(), panel.background=element_blank()) +
  ggtitle("PD vs. controls")
t2 <- melt(as.matrix(m2))
t2$Var1 <- factor(t2$Var1, levels = rev(unique(t2$Var1)))
ggplot(t2) +
  geom_tile(aes(Var2, Var1, fill=value), color = "black") +
  geom_text(aes(Var2, Var1, label=value)) +
  scale_x_discrete(position = "top") +
  scale_fill_gradient(low="white", high = "chocolate", limits = c(0,max_degs)) +
  theme(axis.title=element_blank(), axis.ticks=element_blank(), panel.background=element_blank()) +
  ggtitle("R3 vs. R4/R5")
dev.off()

#####################################################################################
#Volcano plot of fold-change and p-values

volcano.theme <- theme(legend.position = "none",
                       panel.background = element_blank(),
                       axis.line = element_line(colour = "black"),
                       axis.title = element_text(size = 12),
                       plot.title = element_text(size = 12, face = "bold")
)

volcano.plot <- function(t, title){
  diffGenes <- rownames(t)[t$padj < 0.05 & abs(t$log2FoldChange) > 1]
  t$info <- ifelse(rownames(t) %in% diffGenes, 1, 0)
  t$logp <- -log10(t$padj)
  
  p <- ggplot(t, aes(log2FoldChange, logp, colour = info)) +
    geom_point(size = 0.25, alpha = 0.3) +
    labs(x = "Log2 FC", y = "-log10 P") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    ggtitle(title) +
    volcano.theme
  p
}

pdf("volcanoplots.pdf", 3, 2)
lapply(names(data_list), function(n){
  t <- data_list[[n]]
  volcano.plot(t, n)
})
dev.off()

#####################################################################################
# Boxplots of PD genes

# Conversion table gene IDs
conversion <- read.csv("../ahba_entrez2ensembl.txt", sep = "\t")

# PD genes
pd_genes <- c(SCARB2 = "950", ELOVL7 = "79993", SH3GL2 = "6456", 
              SNCA = "6622", BAP1 = "8314", ZNF184 = "7738", SLC18A2 = "6571")
pd_names <- names(pd_genes)
pd_genes <- conversion$ensembl_gene_id[match(pd_genes, conversion$entrezgene)]
names(pd_genes) <- pd_names

box.theme <- theme(panel.background = element_blank(),
                   axis.line = element_line(colour = "black"),
                   axis.text = element_text(size = 10),
                   axis.title = element_text(size = 10),
                   plot.title = element_text(size = 12, face = "bold"),
                   strip.text.x = element_text(size = 12),
                   legend.position = "bottom")

pdf("boxplot_RNAseq_PDgenes.pdf", 2.6, 4)
plots <- lapply(pd_genes, function(g){
  print(names(pd_genes)[match(g, pd_genes)])
  t <- lapply(diseases, function(d){
    t <- lapply(regions, function(r){
      s <- samples[[d]][[r]]
      t <- data.frame(t(expr[g, s]))
      colnames(t) <- "expr"
      t$region <- r
      t$disease <- d
      t
    })
    Reduce(rbind, t)
  })
  t <- Reduce(rbind, t)
  
  t$region <- factor(t$region, levels = c("SN", "GTM"))
  t$disease <- factor(t$disease, levels = rev(unique(t$disease)))
  
  max_y <- max(t$expr)
  min_y <- min(t$expr) 
  range_y <- max_y-min_y
  
  pval <- format(sapply(data_list, function(t){ # DESeq2 results
    t[g, 'padj']
  }), digits = 2, scientific = TRUE)
  pval_region <- pval[3:4]
  pval_disease <- pval[1:2]
  
  p1 <- ggplot(t, aes(x=region, y=expr, fill=region)) + 
    geom_boxplot() +
    geom_signif(comparisons = list(c("SN", "GTM")),
                step_increase = 0.1) +
    scale_y_continuous(limits = c(min_y, max_y+range_y/6)) +
    scale_fill_manual(values = unname(braak.color(levels(t$region)))) +
    labs(x = "", y = "Expression (CPM)") +
    ggtitle(paste(names(pd_genes)[match(g, pd_genes)], "per group")) +
    box.theme +
    facet_grid(~disease)
  pg <- ggplot_build(p1)
  pg$data[[2]]$annotation <- pval_region
  pg$data[[2]]$textsize <- 2.5
  pg$data[[2]]$colour <- ifelse(as.numeric(pg$data[[2]]$annotation) < 0.05, "red", "black")
  q1 <- ggplot_gtable((pg))
  p1 <- plot(q1)
  
  p2 <- ggplot(t, aes(x=disease, y=expr, fill = disease)) + 
    geom_boxplot() +
    geom_signif(comparisons = list(c("PD", "control")),
                step_increase = 0.1) +
    scale_y_continuous(limits = c(min_y, max_y+range_y/6)) +
    scale_fill_manual(values = rep("#FFFFFF",3)) +
    labs(x = "", y = "Expression (CPM)") +
    ggtitle(paste(names(pd_genes)[match(g, pd_genes)], "per region")) +
    box.theme +
    facet_grid(~region)
  pg <- ggplot_build(p2)
  pg$data[[1]]$fill <- unlist(lapply(unname(braak.color(levels(t$region))), function(x)rep(x, 2)))
  pg$data[[2]]$annotation <- pval_disease
  pg$data[[2]]$textsize <- 2.5
  pg$data[[2]]$colour <- ifelse(as.numeric(pg$data[[2]]$annotation) < 0.05, "red", "black")
  q2 <- ggplot_gtable((pg))
  p2 <- plot(q2)
  
})
dev.off()

#####################################################################################
# Box plots for Braak genes
load("../braakGenes.RData")

# Convert entrez IDs to ensembl IDs (and intersect genes in PD data)
bg <- list( # Split Braak genes
    down = braakGenes$entrez_id[braakGenes$r < 0],
    up = braakGenes$entrez_id[braakGenes$r > 0]
)
bg <- lapply(bg, function(g){
    rows <- which(conversion$entrezgene %in% g)
    id <- conversion[rows, "ensembl_gene_id"]
  })

# Check diff. expressed Braak genes
lapply(bg, function(gl){
  sapply(degs, function(x){
    sapply(x, function(y){
      length(intersect(y, gl))
    })
  })
})

pdf("boxplot_RNAseq_braakgenes.pdf", 2.6, 4)
lapply(names(bg), function(dir){
    gl <- bg[[dir]]
    
    t <- lapply(diseases, function(d){
      t <- lapply(regions, function(r){
        s <- samples[[d]][[r]]
        t <- expr[gl, s]
        t <- data.frame('expr' = apply(t, 2, function(x) mean(x, na.rm = TRUE))) # Average across genes
        t$region <- r
        t$disease <- d
        t
      })
      # t <- melt(t)
      Reduce(rbind, t)
    })
    # t <- melt(t)
    t <- Reduce(rbind, t)
    t$region <- factor(t$region, levels = c("SN", "GTM"))
    t$disease <- factor(t$disease, levels = rev(unique(t$disease)))
    
    p1 <- ggplot(t, aes(x=region, y=expr, fill=region)) + 
      geom_boxplot() +
      scale_fill_manual(values = unname(braak.color(levels(t$region)))) +
      labs(x = "", y = "Expression (CPM)") +
      ggtitle(paste0(dir, "regulated \n Braak genes per group")) +
      box.theme +
      facet_grid(~disease)
    plot(p1)
    
    p2 <- ggplot(t, aes(x=disease, y=expr, fill = disease)) + 
      geom_boxplot() +
      scale_fill_manual(values = rep("#FFFFFF",3)) +
      labs(x = "", y = "Expression (CPM)") +
      ggtitle(paste0(dir, "regulated \n Braak genes per region")) +
      box.theme +
      facet_grid(~region)
    pg <- ggplot_build(p2)
    pg$data[[1]]$fill <- unlist(lapply(unname(braak.color(levels(t$region))), function(x)rep(x, 2)))
    q2 <- ggplot_gtable((pg))
    p2 <- plot(q2)
})
dev.off()

# #####################################################################################
# # check presence of genes in module M47
# m47 <- unlist(read.table("../m47_genes.txt", colClasses = 'character'))
# m47 <- conversion$ensembl_gene_id[match(m47, conversion$entrezgene)]
# m47 <- m47[!is.na(m47)]

########## PSEA: cell-type correction ##########

# Cell-type genes
celltypes <- sapply(c("Neurons", "Astrocytes", "Oligodendrocytes", "Microglia", "Endothelial_cells"), function(type){
  file = paste0("../brainscope_celltypes/", type, ".txt")
  entrez_ids <- as.character(read.csv(file, header = TRUE)$entrez_id)
  ensembl_ids <- conversion$ensembl_gene_id[match(entrez_ids, conversion$entrezgene)]
  ensembl_ids[!is.na(ensembl_ids)]
}, simplify = FALSE)

# Cell-type mean expression 
ct_mean <- t(sapply(celltypes, function(ct){
  x <- expr[ct, ]
  colMeans(x, na.rm = TRUE)
}))

psea <- function(gene, ct, groups){
  
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
  
  # fit <- lm(gene ~ neurons + astrocytes + oligodendrocytes + microglia + endothelial_cells, subset = which(groups==0))
  # # par(mfrow=c(2,3), mex=0.8)
  # # crplot(fit, "neurons", newplot = FALSE)
  # # crplot(fit, "astrocytes", newplot = FALSE)
  # # crplot(fit, "oligodendrocytes", newplot = FALSE)
  # # crplot(fit, "microglia", newplot = FALSE)
  # # crplot(fit, "endothelial_cells", newplot = FALSE)
  # summary <- summary(fit)
  # coefficients <- summary$coefficients
  # celltype_fit <- coefficients[-1, c(1,4)]
  # colnames(celltype_fit) <- c("celltype_beta", "celltype_pval")
  
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
  # if (length(pval) !=5) break
  # cbind(celltype_fit, group_fc = fc, group_pval = pval)
  cbind(group_pval = pval, group_fc = fc)
}

geneset <- conversion$ensembl_gene_id[conversion$entrezgene %in% braakGenes$entrez_id]
geneset <- intersect(genes, geneset)
  
# PSEA per region between disease groups
brgs_regions <- sapply(regions, function(r){
  s <- sample_info$code[sample_info$region == r]
  groups <- sample_info[s, "disease"]
  groups <- as.numeric(groups == "PD")
  
  # BRGs
  psea_brgs <- sapply(geneset, function(g){
    gene <- unlist(expr[g, s])
    ct <- ct_mean[, s]
    psea(gene, ct, groups)
  }, simplify = FALSE)
  psea_brgs <- simplify2array(psea_brgs)
  
  # # Celltype-specific expression
  # m1 <- t(psea_brgs[, "celltype_pval", ])
  # m1 <- apply(m1, 2, function(x)p.adjust(x, method = "BH"))
  # colSums(m1 < 0.05, na.rm = T)
  
  # Group-dependent expression
  m2 <- t(psea_brgs[, "group_pval", ])
  m2 <- apply(m2, 2, function(x) p.adjust(x, method = "BH"))
  colSums(m2 < 0.05, na.rm = T)
})

# PSEA per patient group between regions
brgs_diseases <- sapply(diseases, function(d){
  s <- sample_info$code[sample_info$disease == d]
  groups <- sample_info[s, "region"]
  groups <- as.numeric(groups == "GTM")
  
  # BRGs
  psea_brgs <- sapply(geneset, function(g){
    gene <- unlist(expr[g, s])
    ct <- ct_mean[, s]
    psea(gene, ct, groups)
  }, simplify = FALSE)
  psea_brgs <- simplify2array(psea_brgs)
  
  # # Celltype-specific expression
  # m1 <- t(psea_brgs[, "celltype_pval", ])
  # m1 <- apply(m1, 2, function(x)p.adjust(x, method = "BH"))
  # colSums(m1 < 0.05, na.rm = T)
  
  # Group-dependent expression
  m2 <- t(psea_brgs[, "group_pval", ])
  m2 <- apply(m2, 2, function(x) p.adjust(x, method = "BH"))
  colSums(m2 < 0.05, na.rm = T)
})
max_brgs <- max(brgs_diseases, brgs_regions)

pdf("heatmap_diffgenes_BRGs.pdf", 3.3, 2)
m <- brgs_regions
rownames(m) <- gsub("_", " ", rownames(m))
colnames(m) <- c("R3", "R4/R5")
t <- melt(m)
t$Var2 <- factor(t$Var2, levels = rev(unique(t$Var2)))
ggplot(t) +
  geom_tile(aes(Var2, Var1, fill=value), color = "black") +
  geom_text(aes(Var2, Var1, label=value)) +
  scale_x_discrete(position = "top") +
  scale_fill_gradient(low="white", high = "chocolate", limits = c(0,max_brgs)) +
  theme(axis.title=element_blank(), axis.ticks=element_blank(), panel.background=element_blank()) +
  ggtitle("PD vs. controls")
m <- brgs_diseases
rownames(m) <- gsub("_", " ", rownames(m))
t <- melt(m)
t$Var2 <- factor(t$Var2, levels = rev(unique(t$Var2)))
ggplot(t) +
  geom_tile(aes(Var2, Var1, fill=value), color = "black") +
  geom_text(aes(Var2, Var1, label=value)) +
  scale_x_discrete(position = "top") +
  scale_fill_gradient(low="white", high = "chocolate", limits = c(0,max_brgs)) +
  theme(axis.title=element_blank(), axis.ticks=element_blank(), panel.background=element_blank()) +
  ggtitle("R1 vs. R3")
dev.off()
