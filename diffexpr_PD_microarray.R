setwd("M:/doorgeefluik/Arlin Keo doorgeefuik/Arrays data_coded for Arlin")
options(stringsAsFactors = FALSE)
library(ggplot2)
library(ggsignif)
library(reshape2)
library(RColorBrewer)
library(PSEA)

########## Prepare data ##########

# Braak region colors
braakColors <- brewer.pal(6, "Set2")
names(braakColors) <- paste0("R", c(1:6))

braak.color <- function(x) {
  unname(unlist(lapply(x, function(r){
    if (r == "MO") braakColors[1]
    else if (r == "LC") braakColors[2]
    else if (r == "SN") braakColors[3]
    else "#FFFFFF"
  })))
}

# Load expression data genes (entrez IDs) x samples
expr <- read.csv("rma_genes.csv")
rownames(expr) <- unlist(expr[,1])
expr <- expr[, -1]

#Load sample info
info <- read.csv("clinicalTraits_CODED_dis_num_batch_no commas.csv")
info$Region <- sapply(info$Region, function(r){
  if (r == 1) "OB"
  else if (r == 2) "SN"
  else if (r == 3) "LC"
  else "MO"
})
info$Disease <- sapply(info$Disease, function(r){
  if (r == 0) "control"
  else if (r == 1) "iLBD"
  else "PD"
})
rownames(info) <- info$Array.code

genes <- rownames(expr)
regions <- unique(info$Region)
names(regions) <- regions
diseases <- sort(unique(info$Disease))
names(diseases) <- diseases

# Sample IDs per region and disease group
samples <- sapply(regions, function(r){
  sapply(diseases, function(d){
    info$Array.code[info$Disease == d & info$Region == r]
  }, simplify = FALSE)
}, simplify = FALSE)
t(sapply(samples, function(x)sapply(x, length))) # print sample size table

########## T-test between disease groups and regions ##########

# For each region do t-test between disease groups
ttest_disease <- sapply(regions, function(r){
  print(r)
  # T-test controls vs. iLBD and controls vs. PD
  s <- lapply(diseases, function(d){
    samples[[r]][[d]]
  })
  tab <- as.data.frame(t(sapply(genes, function(g){
    a <- expr[g, s[["control"]]]
    b <- expr[g, s[["iLBD"]]]
    c <- expr[g, s[["PD"]]]
    ttest1 <- t.test(a,b)
    ttest2 <- t.test(a,c)
    c(ilbd.meanDiff = unname(ttest1$estimate[1] - ttest1$estimate[2]), ilbd.pvalue = unname(ttest1$p.value), 
      pd.meanDiff = unname(ttest2$estimate[1] - ttest2$estimate[2]), pd.pvalue = unname(ttest2$p.value))
  })))
  # Corrected P-values
  tab$ilbd.pvalue.BH <- p.adjust(tab$ilbd.pvalue, method = "BH")
  tab$pd.pvalue.BH <- p.adjust(tab$pd.pvalue, method = "BH")
  
  colnames(tab) <- paste0(r, ".", colnames(tab)) # Add region to column name
  tab[, c(1,2,5,3,4,6)]
}, simplify = FALSE)

# For each subject type, T-test between regions
ttest_region <- sapply(diseases, function(d){
  print(d)
  # T-test between regions
  s <- lapply(regions, function(r){
    samples[[r]][[d]]
  })
  tab <- as.data.frame(t(sapply(genes, function(g){
    a <- expr[g, s[["MO"]]]
    b <- expr[g, s[["LC"]]]
    c <- expr[g, s[["SN"]]]
    ttest1 <- t.test(a,b)
    ttest2 <- t.test(a,c)
    ttest3 <- t.test(b,c)
    c(mo.lc.meanDiff = unname(ttest1$estimate[1] - ttest1$estimate[2]), mo.lc.pvalue = unname(ttest1$p.value), 
      mo.sn.meanDiff = unname(ttest2$estimate[1] - ttest2$estimate[2]), mo.sn.pvalue = unname(ttest2$p.value),
      lc.sn.meanDiff = unname(ttest3$estimate[1] - ttest3$estimate[2]), lc.sn.pvalue = unname(ttest3$p.value))
  })))
  # Corrected P-values
  tab$mo.lc.pvalue.BH <- p.adjust(tab$mo.lc.pvalue, method = "BH")
  tab$mo.sn.pvalue.BH <- p.adjust(tab$mo.sn.pvalue, method = "BH")
  tab$lc.sn.pvalue.BH <- p.adjust(tab$lc.sn.pvalue, method = "BH")
  
  colnames(tab) <- paste0(d, ".", colnames(tab)) # Add disease type to column name
  tab[, c(1,2,7,3,4,8,5,6,9)]
}, simplify = FALSE)

# Combine tables for all regions and save as csv
tab_all_disease <- Reduce(cbind, ttest_disease)
write.table(tab_all_disease, file = "ttest_disease.txt", quote = FALSE, sep = "\t")

# Combine tables for all regions and save as csv
tab_all_regions <- Reduce(cbind, ttest_regions)
write.table(tab_all_regions, file = "ttest_regions.txt", quote = FALSE, sep = "\t")

# Number of DEGs between diseases
sapply(ttest_disease, function(t){
  ilbd_up <- sum(t[, 1]>1 & t[, 3]<0.05)
  ilbd_down <- sum(t[, 1]< -1 & t[, 3]<0.05)
  pd_up <- sum(t[, 4]>1 & t[, 6]<0.05)
  pd_down <- sum(t[, 4]< -1 & t[, 6]<0.05)
  c(ilbd_down=ilbd_down, ilbd_up=ilbd_up, pd_down=pd_down, pd_up=pd_up)
})

# Number of DEGs between regions
sapply(ttest_region, function(t){
  mo_lc_up <- sum(t[, 1]>1 & t[, 3]<0.05)
  mo_lc_down <- sum(t[, 1]< -1 & t[, 3]<0.05)
  mo_sn_up <- sum(t[, 4]>1 & t[, 6]<0.05)
  mo_sn_down <- sum(t[, 4]< -1 & t[, 6]<0.05)
  lc_sn_up <- sum(t[, 7]>1 & t[, 9]<0.05)
  lc_sn_down <- sum(t[, 7]< -1 & t[, 9]<0.05)
  c(mo_lc_up=mo_lc_up, mo_lc_down=mo_lc_down, mo_sn_up=mo_sn_up, mo_sn_down=mo_sn_down, lc_sn_up=lc_sn_up, lc_sn_down=lc_sn_down)
})

##############################################################################
# Heatmap number of differentially expressed genes (t-test)

# Number of DEGs between diseases
degs_disease <- sapply(ttest_disease, function(t){
  ilbd <- sum(abs(t[, 1])>1 & t[, 3]<0.05)
  pd <- sum(abs(t[, 4])>1 & t[, 6]<0.05)
  c(ilbd=ilbd, pd=pd)
})

# Number of DEGs between regions
degs_regions <- sapply(ttest_region, function(t){
  mo_lc <- sum(abs(t[, 1])>1 & t[, 3]<0.05)
  mo_sn <- sum(abs(t[, 4])>1 & t[, 6]<0.05)
  lc_sn <- sum(abs(t[, 7])>1 & t[, 9]<0.05)
  c(mo_lc=mo_lc, mo_sn=mo_sn, lc_sn=lc_sn)
})
max_degs <- max(degs_regions, degs_disease)

pdf("heatmap_diffgenes.pdf", 3, 2)
lapply(diseases, function(d){
  m <- matrix(NA, 3, 3)
  rownames(m) <- paste0("R", c(1:3))
  colnames(m) <- paste0("R", c(1:3))
  m[1,2] <- degs_regions[1, d] 
  m[1,3] <- degs_regions[2, d]
  m[2,3] <- degs_regions[3, d]
  m[lower.tri(m)] <- m[upper.tri(m)]
  diag(m) <- 0
  t <- melt(m)
  t$Var2 <- factor(t$Var2, levels = rev(unique(t$Var2)))
  ggplot(t) +
    geom_tile(aes(Var1, Var2, fill=value), color = "black") +
    geom_text(aes(Var1, Var2, label=value)) +
    scale_x_discrete(position = "top") +
    scale_fill_gradient(low="white", high = "chocolate", limits = c(0,max_degs)) +
    theme(axis.title=element_blank(), axis.ticks=element_blank(), panel.background=element_blank()) +
    ggtitle(d)
})
m <- degs_disease[,c(3,4,2)]
rownames(m) <- c("iLBD", "PD")
colnames(m) <- paste0("R", c(1:3))
t <- melt(m)
t$Var2 <- factor(t$Var2, levels = rev(unique(t$Var2)))
ggplot(t) +
  geom_tile(aes(Var1, Var2, fill=value), color = "black") +
  geom_text(aes(Var1, Var2, label=value)) +
  scale_x_discrete(position = "top") +
  scale_fill_gradient(low="white", high = "chocolate", limits = c(0,max_degs)) +
  theme(axis.title=element_blank(), axis.ticks=element_blank(), panel.background=element_blank())
dev.off()

########## Volcano plot of fold-change and p-values ########## 

volcano.theme <- theme(legend.position = "none",
               panel.background = element_blank(),
               axis.line = element_line(colour = "black"),
               axis.title = element_text(size = 12),
               plot.title = element_text(size = 12, face = "bold")
)

volcano.plot <- function(t, title){
  diffGenes <- rownames(t)[t$BH < 0.05 & abs(t$meanDiff) > 1]
  t$info <- ifelse(rownames(t) %in% diffGenes, 1, 0)
  t$logp <- -log10(t$pvalue)
  
  p <- ggplot(t, aes(meanDiff, logp, colour = info)) +
    geom_point(size = 0.25, alpha = 0.3) +
    labs(x = "Fold-change", y = "-log10 p-value") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    ggtitle(title) +
    volcano.theme
  p
}

pdf("volcanoplots_disease.pdf", 3, 2)
lapply(regions, function(r){
  tab <- ttest_disease[[r]]
  tabll <- list(iLBD = tab[,c(1:3)], PD = tab[,c(4:6)])
  
  lapply(names(tabll), function(d){
    t <- tabll[[d]]
    colnames(t) <- c("meanDiff", "pvalue", "BH")
    title <- paste0("controls vs. ", d, " in ", r)
    volcano.plot(t, title)
  })
})
dev.off()

pdf("volcanoplots_regions.pdf", 3, 2)
lapply(diseases, function(d){
  tab <- ttest_region[[d]]
  tabll <- list(mo.lc = tab[,c(1:3)], mo.sn = tab[,c(4:6)], lc.sn = tab[,c(7:9)])
  
  lapply(names(tabll), function(r){
    t <- tabll[[r]]
    colnames(t) <- c("meanDiff", "pvalue", "BH")
    title <- paste(paste(unlist(strsplit(r, split = "\\.")), collapse = " vs. "), "in", d)
    volcano.plot(t, title)
  })
})
dev.off()

#####################################################################################
# Boxplots of PD genes

# Box plot of PD Braak genes
pd_genes <- c(SCARB2 = "950", ELOVL7 = "79993", SH3GL2 = "6456", 
              SNCA = "6622", BAP1 = "8314", ZNF184 = "7738", SLC18A2 = "6571")

# sizes <- sapply(samples, function(x)sapply(x, length))[c(3:1), c("OB", "MO", "LC", "SN")]

box.theme <- theme(panel.background = element_blank(),
                   axis.line = element_line(colour = "black"),
                   axis.text = element_text(size = 10),
                   axis.title = element_text(size = 10),
                   plot.title = element_text(size = 12, face = "bold"),
                   strip.text.x = element_text(size = 12),
                   legend.position = "bottom")

pdf("boxplot_microarray_PDgenes.pdf", 5, 4)
plots <- lapply(pd_genes[-6], function(g){
  print(names(pd_genes)[match(g, pd_genes)])
  t <- lapply(diseases, function(d){
    s <- unlist(lapply(samples, function(x) x[[d]]))
    data.frame(expr = unlist(expr[g, s]), region = info[s, "Region"])
  })
  t <- melt(t)
  colnames(t) <- c("region", "variable", "expr", "disease")
  t <- t[t$region != "OB", ]
  t$region <- factor(t$region, levels = c("MO", "LC", "SN"))
  t$disease <- factor(t$disease, levels = rev(unique(t$disease)))
  
  max_y <- max(t$expr)
  min_y <- min(t$expr) 
  range_y <- max_y-min_y
  y_pos <- c(max_y+range_y/12, max_y+range_y/5, max_y+range_y/3)
  pval_region <- unlist(sapply(rev(diseases), function(d){
    format(unlist(ttest_region[[d]][g, c(3,6,9)]), digits = 2, scientific = TRUE)
  }, simplify = FALSE))
  pval_disease <- unlist(sapply(regions[c(3,4,2)], function(r){
    format(unlist(ttest_disease[[r]][g, c(6,3)]), digits = 2, scientific = TRUE)
  }, simplify = FALSE))
  
  p1 <- ggplot(t, aes(x=region, y=expr, fill=region)) + 
    geom_boxplot() +
    geom_signif(comparisons = list(c("MO", "LC"), c("LC", "SN"), c("MO", "SN")),
                step_increase = 0.1) +
    scale_y_continuous(limits = c(min_y, max_y+range_y/3)) +
    scale_fill_manual(values = unname(braak.color(levels(t$region)))) +
    labs(x = "", y = "Expression \n (log2-transformed)") +
    ggtitle(paste(names(pd_genes)[match(g, pd_genes)], "per group")) +
    box.theme +
    facet_grid(~disease)
  pg <- ggplot_build(p1)
  pg$data[[2]]$annotation <- unlist(lapply(pval_region, function(x)rep(x, 3)))
  pg$data[[2]]$textsize <- 2.5
  pg$data[[2]]$colour <- ifelse(as.numeric(pg$data[[2]]$annotation) < 0.05, "red", "black")
  q1 <- ggplot_gtable((pg))
  p1 <- plot(q1)
  
  p2 <- ggplot(t, aes(x=disease, y=expr, fill = disease)) + 
    geom_boxplot() +
    geom_signif(comparisons = list(c("iLBD", "control"), c("PD", "control")),
                step_increase = 0.1) +
    scale_y_continuous(limits = c(min_y, max_y+range_y/3)) +
    scale_fill_manual(values = rep("#FFFFFF",3)) +
    labs(x = "", y = "Expression \n (log2-transformed)") +
    ggtitle(paste(names(pd_genes)[match(g, pd_genes)], "per region")) +
    box.theme +
    facet_grid(~region)
  pg <- ggplot_build(p2)
  pg$data[[1]]$fill <- unlist(lapply(unname(braak.color(levels(t$region))), function(x)rep(x, 3)))
  pg$data[[2]]$annotation <- unlist(lapply(pval_disease, function(x)rep(x, 3)))
  pg$data[[2]]$textsize <- 2.5
  pg$data[[2]]$colour <- ifelse(as.numeric(pg$data[[2]]$annotation) < 0.05, "red", "black")
  q2 <- ggplot_gtable((pg))
  p2 <- plot(q2)
  
})
dev.off()

#####################################################################################
# Box plots for Braak genes
load("../braakGenes.RData")
# load("../braakGenes2.RData")
# load("../braakGenes3.RData")

bg <- list(
  bg1 = list( # Braak genes selected WITHOUT cell-type correction
    down = braakGenes$entrez_id[braakGenes$r < 0],
    up = braakGenes$entrez_id[braakGenes$r > 0]
  )#,
  # bg2 = list( # Braak genes selected WITH cell-type correction
  #   down = braakGenes2$entrez_id[braakGenes2$braak6 < 0],
  #   up = braakGenes2$entrez_id[braakGenes2$braak6 > 0]
  # ),
  # bg3 = braakGenes3 # Intersection of corrected and uncorrected results
)

# Convert entrez IDs to ensembl IDs (and intersect genes in PD data)
bg <- lapply(bg, function(s){
  lapply(s, function(g){
    intersect(g, genes)
  })
})

# Check diff. expressed Braak genes
sapply(bg$bg1, function(l){
  print(sapply(ttest_disease, function(t){
    ilbd_up <- rownames(t)[t[, 1]>1 & t[, 3]<0.05]
    ilbd_down <- rownames(t)[t[, 1]< -1 & t[, 3]<0.05]
    pd_up <- rownames(t)[t[, 4]>1 & t[, 6]<0.05]
    pd_down <- rownames(t)[t[, 4]< -1 & t[, 6]<0.05]
    c(ilbd_down=length(intersect(ilbd_down, l)), 
      ilbd_up=length(intersect(ilbd_up, l)), 
      pd_down=length(intersect(pd_down, l)), 
      pd_up=length(intersect(pd_up,l)))
  }))
  sapply(ttest_region, function(t){
    mo_lc_up <- rownames(t)[t[, 1]>1 & t[, 3]<0.05]
    mo_lc_down <- rownames(t)[t[, 1]< -1 & t[, 3]<0.05]
    mo_sn_up <- rownames(t)[t[, 4]>1 & t[, 6]<0.05]
    mo_sn_down <- rownames(t)[t[, 4]< -1 & t[, 6]<0.05]
    lc_sn_up <- rownames(t)[t[, 7]>1 & t[, 9]<0.05]
    lc_sn_down <- rownames(t)[t[, 7]< -1 & t[, 9]<0.05]
    c(mo_lc_up=length(intersect(mo_lc_up,l)), 
      mo_lc_down=length(intersect(mo_lc_down,l)), 
      mo_sn_up=length(intersect(mo_sn_up,l)), 
      mo_sn_down=length(intersect(mo_sn_down,l)), 
      lc_sn_up=length(intersect(lc_sn_up,l)),
      lc_sn_down=length(intersect(lc_sn_down,l)))
  })
}, simplify = FALSE)

pdf("boxplot_microarray_braakgenes.pdf", 5, 4)
lapply(names(bg), function(n){
  bgs <- bg[[n]]
  
  lapply(names(bgs), function(dir){
    gl <- bgs[[dir]]
    
    t <- lapply(diseases, function(d){
      s <- unlist(lapply(samples, function(x) x[[d]]))
      t <- expr[gl, s]
      mean <- apply(t, 2, mean) # Average across genes
      data.frame(expr = mean, region = info[s, "Region"])
    })
    t <- melt(t)
    colnames(t) <- c("region", "variable", "expr", "disease")
    t <- t[t$region != "OB", ]
    t$region <- factor(t$region, levels = c("MO", "LC", "SN"))
    t$disease <- factor(t$disease, levels = rev(unique(t$disease)))
    
    p1 <- ggplot(t, aes(x=region, y=expr, fill=region)) + 
      geom_boxplot() +
      scale_fill_manual(values = unname(braak.color(levels(t$region)))) +
      labs(x = "", y = "Expression \n (log2-transformed)") +
      ggtitle(paste0(n, "\n", dir, "regulated Braak genes per group")) +
      box.theme +
      facet_grid(~disease)
    plot(p1)
    
    p2 <- ggplot(t, aes(x=disease, y=expr, fill = disease)) + 
      geom_boxplot() +
      scale_fill_manual(values = rep("#FFFFFF",3)) +
      labs(x = "", y = "Expression \n (log2-transformed)") +
      ggtitle(paste0(n, "\n", dir, "regulated Braak genes per region")) +
      box.theme +
      facet_grid(~region)
    pg <- ggplot_build(p2)
    pg$data[[1]]$fill <- unlist(lapply(unname(braak.color(levels(t$region))), function(x)rep(x, 3)))
    q2 <- ggplot_gtable((pg))
    p2 <- plot(q2)
  })
})
dev.off()

########## PSEA: cell-type correction ##########

# Cell-type genes
celltypes <- sapply(c("Neurons", "Astrocytes", "Oligodendrocytes", "Microglia", "Endothelial_cells"), function(type){
  file = paste0("../brainscope_celltypes/", type, ".txt")
  as.character(read.csv(file, header = TRUE)$entrez_id)
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
  
  fit <- lm(gene ~ neurons + astrocytes + oligodendrocytes + microglia + endothelial_cells, subset = which(groups==0))
  # par(mfrow=c(2,3), mex=0.8)
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
  # if (length(pval) !=5) break
  # cbind(celltype_fit, group_fc = fc, group_pval = pval)
  cbind(group_pval = pval, group_fc = fc)
}

# PSEA per region between disease groups
brgs_regions <- sapply(regions, function(r){
  samples <- info$Array.code[info$Region == r & info$Disease %in% c("control", "PD")]
  groups <- info[samples, "Disease"]
  groups <- as.numeric(groups == "PD")

  # BRGs
  psea_brgs <- sapply(intersect(braakGenes$entrez_id, genes), function(g){
    gene <- unlist(expr[g, samples])
    ct <- ct_mean[, samples]
    psea(gene, ct, groups)
  }, simplify = FALSE)
  psea_brgs <- simplify2array(psea_brgs)
  
  # # Celltype-specific expression
  # m1 <- t(psea_brgs[, "celltype_pval", ])
  # m1 <- apply(m1, 2, function(x)p.adjust(x, method = "BH"))
  # colSums(m1 < 0.05)
  
  # Group-dependent expression
  m2 <- t(psea_brgs[, "group_pval", ])
  m2 <- apply(m2, 2, function(x) p.adjust(x, method = "BH"))
  colSums(m2 < 0.05)
})

# PSEA per patient group between regions
brgs_diseases <- sapply(diseases, function(d){
  samples <- info$Array.code[info$Disease == d & info$Region %in% c("MO", "SN")]
  groups <- info[samples, "Region"]
  groups <- as.numeric(groups == "SN")
  
  # BRGs
  psea_brgs <- sapply(intersect(braakGenes$entrez_id, genes), function(g){
    gene <- unlist(expr[g, samples])
    ct <- ct_mean[, samples]
    psea(gene, ct, groups)
  }, simplify = FALSE)
  psea_brgs <- simplify2array(psea_brgs)
  
  # Celltype-specific expression
  # m1 <- t(psea_brgs[, "celltype_pval", ])
  # m1 <- apply(m1, 2, function(x)p.adjust(x, method = "BH"))
  # colSums(m1 < 0.05)
  
  # Group-dependent expression
  m2 <- t(psea_brgs[, "group_pval", ])
  m2 <- apply(m2, 2, function(x) p.adjust(x, method = "BH"))
  colSums(m2 < 0.05)
})
max_brgs <- max(brgs_diseases, brgs_regions)

# Heatmap
pdf("heatmap_diffgenes_BRGs.pdf", 4, 2)
m <- brgs_regions[, c(3,4,2)]
rownames(m) <- gsub("_", " ", rownames(m))
colnames(m) <- paste0("R", c(1:3))
t <- melt(m)
t$Var1 <- factor(t$Var1, levels = rev(unique(t$Var1)))
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
t$Var1 <- factor(t$Var1, levels = rev(unique(t$Var1)))
t$Var2 <- factor(t$Var2, levels = rev(unique(t$Var2)))
ggplot(t) +
  geom_tile(aes(Var2, Var1, fill=value), color = "black") +
  geom_text(aes(Var2, Var1, label=value)) +
  scale_x_discrete(position = "top") +
  scale_fill_gradient(low="white", high = "chocolate", limits = c(0,max_brgs)) +
  theme(axis.title=element_blank(), axis.ticks=element_blank(), panel.background=element_blank()) +
  ggtitle("R1 vs. R3")
dev.off()