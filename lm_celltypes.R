# Reference signal for cell-type correction
setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")
library(reshape2)
library(ggplot2)
library(gplots)
library(plyr)
load("resources/braakInfo.RData")
load("../ABA_Rdata/BrainExpr.RData")

# Cell-type markers from different studies
celltypes <- list(
  Brainscope = sapply(c("Neurons", "Astrocytes", "Oligodendrocytes", "Microglia", "Endothelial_cells"), function(type){
    file = paste0("brainscope_celltypes/", type, ".txt")
    as.character(read.csv(file, header = TRUE)$entrez_id)
  }, simplify = FALSE),
  Darmanis2015 = sapply(c("Neurons", "Astrocytes", "Oligodendrocytes", "Microglia", "Endothelial_cells", "OPC"), function(type){
    conn <- file(paste0("Darmanis2015_celltypes/", type, ".txt"), "r")
    lines <- readLines(conn)
    close(conn)
    name2EntrezId(lines)
  }, simplify = FALSE),
  PSEA = lapply(list(Neurons = c("NEFL", "ENO2", "SLC12A5", "KCNQ2", "SCN3A"),
                     Astrocytes = c("GFAP", "AQP4", "GJA1"),
                     Oligodendrocytes = c("MOG", "MAG", "MOBP", "MBP"),
                     Microglia = c("CD37", "CD53"),
                     Endothelial_cells = c("PECAM1")), name2EntrezId)
)

# Cell-type colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
ggcolor <- gg_color_hue(6)
names(ggcolor) <- unique(unlist(sapply(celltypes, names)))

# Correlation heatmap of marker genes
pdf("celltype_markers_correlation_heatmap.pdf", 12, 8)
lapply(names(celltypes), function(n){
  print(n)
  s <- celltypes[[n]]
  markers <- unlist(s)
  colors <- unlist(lapply(names(s), function(c){
    rep(ggcolor[c], length(s[[c]]))
  }))
  
  lapply(donorNames[1], function(d){
    x <- t(brainExpr[[d]][markers, ])
    r <- cor(x)
    df <- melt(r)
    df$Var1 <- factor(df$Var1, levels = unique(df$Var1))
    df$Var2 <- factor(df$Var2, levels = unique(df$Var2))
    ggplot(df, aes(x=Var1, y=Var2, fill = value)) + 
      geom_tile() +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                           midpoint = 0, limit = c(-1, 1), space = "Lab", 
                           name="r") +
      theme(axis.text = element_text(color = colors)) +
      ggtitle(n)
  })
})
dev.off()

# Variance explained by PCA and SVD
ref_var <- lapply(celltypes, function(s){
  simplify2array(lapply(donorNames, function(d){
    t <- sapply(s, function(ct){
      x <- t(brainExpr[[d]][ct, ])
      pca <- summary(prcomp(x))$importance[2,1] # Variance explained
      svd <- summary(prcomp(x, center = FALSE))$importance[2,1]
      c(pca = pca, svd = svd)
    })
  }))
  adply(ref_var, c(3,1))  
})

# Cell-type reference signal based on PCA, SVD, and mean
ref_signal <- lapply(celltypes, function(s){ # For set of markers per study
  simplify2array(lapply(donorNames, function(d){
    simplify2array(lapply(s, function(ct){
      x <- t(brainExpr[[d]][ct, ])
      mean <- apply(x, 1, mean)
      pca <- prcomp(x)$x[, 1]# 1st PC (eigen gene expr)
      pca <- if (cor(pca, mean) > 0) pca else -pca # flip sign of eigen gene based on the data
      svd <- prcomp(x, center = FALSE)$x[, 1]# 1st PC (eigen gene expr)
      svd <- if (cor(svd, mean) > 0) svd else -svd # flip sign of eigen gene based on the data
      cbind(mean, pca, svd)
    })) # 3D-array: samples x method x cell-types
  }))
})

# Correlation PCA, SVD, and mean
cor_ref <- lapply(celltypes, function(s){
  lapply(donorNames, function(d){
    apply(ref_signal$Brainscope[[d]], 3, function(ct){
      r <- cor(ct)
      c(mean_vs_pca = r[1,2], mean_vs_svd = r[1,3], pca_vs_svd = r[2,3])
    })
  })
  cor_ref <- simplify2array(cor_ref)
  adply(cor_ref, c(3,1))
})

# Plot celltype reference expression
plot.ref.signal <- function(data, title, vline = NULL){
  colnames(data) <- make.names(colnames(data), unique = TRUE)
  df <- melt(data)
  colnames(df) <- c("celltype", "sample", "expr")
  df$sample <- factor(df$sample, levels = unique(df$sample))
  ggplot(df, aes(x=sample, y=expr, color=celltype)) + 
    geom_point() +
    geom_smooth(aes(group=celltype)) +
    geom_vline(xintercept = vline) +
    ggtitle(title) +c
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
}

pdf("celltype_markers_lineplot.pdf", 8, 6)
lapply(names(ref_signal), function(s){
  lapply(names(s), function(d){
    idx <- unlist(braak_idx[[d]]) # for samples
    intercepts <- match(c(1:6), braakLabels[[d]][idx])[-1]-1
    lapply(dimnames(ref)[[2]], function(m){
      df <- t(ref_signal[[s]][[d]][idx,m,])
      plot.ref.signal(df, paste(s, d, m, sep = "; "), intercepts)
    })
  })
})
dev.off()

# Correlation of cell-type markers with cell-type reference signal
pdf("celltype_markers_refcorrelation_hist_markersonly.pdf", 12, 8)
lapply(names(ref_signal), function(s){
  t <- lapply(donorNames, function(d){
    ref <- ref_signal[[s]][[d]]
    expr <- brainExpr[[d]]
    t <- melt(sapply(dimnames(ref)[[3]], function(ct){
      markers <- celltypes[[s]][[ct]]
      ref_expr <- ref[,,ct]
      data_expr <- t(expr[markers,])
      t(cor(ref_expr, data_expr))
    }, simplify = FALSE))
    cbind(t, d)
  })
  df <- Reduce(rbind, t)
  colnames(df) <- c("gene", "ref", "r", "celltype", "donor")
  df$celltype <- factor(df$celltype, levels = unique(df$celltype))
  df$donor <- factor(df$donor, levels = donorNames)
  p <- ggplot(df, aes(x = r, fill = ref)) +
    geom_histogram(binwidth = 0.05, position = "dodge") +
    geom_vline(xintercept = 0) +
    ggtitle(paste("Gene correlations with reference signal", s, d, sep = '; ')) + 
    facet_grid(donor~celltype) +
    theme_classic() +
    theme(strip.text = element_text(size=12),
          strip.background = element_blank(),
          axis.text = element_text(size = 10))
})
dev.off()

###############################################################################################

# Fit linear model Equation 1
lm_fit1 <- lapply(donorNames, function(d){
  braak <- braakLabels[[d]] # braak labels
  ct <- ref_signal[[d]][, "mean", ] # samples x cell-types
  expr <- t(brainExpr[[d]]) # samples x genes
  fit <- lm(expr ~ ct)
})

# Fit linear model Equation 2
lm_fit2 <- lapply(donorNames, function(d){
  idx <- unlist(braak_idx[[d]]) # idx for samples
  braak <- braakLabels[[d]][idx] # braak labels
  ct <- ref_signal[[d]][idx, "mean", ] # samples x cell-types
  expr <- t(brainExpr[[d]][, idx]) # samples x genes
  fit <- lm(expr ~ ct + braak)
})

########## Scatter plots ##########
g <- "85300"#name2EntrezId("SNCA")

pdf("EQ1and2_ATCAY_neuronal_marker.pdf", 6, 4)
lapply(names(celltypes), function(c){
  lapply(donorNames, function(d){
    samples <- unlist(braak_idx[[d]])
    label <- braakLabels[[d]][samples]
    ct <- mean_celltype[[d]][, samples]
    expr <- brainExpr[[d]][, samples]
    color <- sampleInfo[[d]][samples, "color_hex_triplet"]
    
    ref <- ct[c, ] # Cell-type expression
    gene <- unlist(expr[g, ])
    coef1 <- c(coef(lm_fit1[[d]])[, g], braak1 = 0)
    coef2 <- c(coef(lm_fit2[[d]])[, g], braak1 = 0)
    df <- data.frame(ref, gene, factor(label))
    # df$braakcoef <- coef1[paste0("braak", label)]
    ggplot(df) +
      geom_point(aes(x=ref, y=gene, shape = label), color = color) +
      geom_abline(intercept = coef1[1], slope = coef1[paste0("ct", c)], size = 1, linetype = 7) +# For cell-type estimate
      geom_abline(intercept = coef2[1], slope = coef2[paste0("ct", c)], size = 1, linetype = 7, color = "red") +# For cell-type estimate
      # geom_abline(aes(intercept = coef[1], slope = braakcoef, linetype = label), size = 1) + # braak region estimates
      # geom_smooth(aes(x=ref, y=gene), method = lm, se = FALSE, fullrange = TRUE) +
      scale_shape_manual(values=1:length(unique(label))) +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      labs(x = c, y = entrezId2Name(g)) +
      ggtitle(d) +
      expand_limits(x=0, y=min(coef1[1], coef2[1])-2) +
      theme
  })
})
dev.off()

# Plot concatenated expression data
concat.samples <- function(expr, idx = NULL){
  isVector <- is.vector(expr[[1]])
  res <- lapply(donorNames, function(d){
    e <- expr[[d]]
    if (!is.null(idx)){
      samples <- unlist(idx[[d]])
      if (isVector) e[samples] else e[, samples]
    } 
    else e
  })
  if (isVector) unlist(res) else Reduce(cbind, res)
}

braakExpr <- concat.samples(brainExpr, braak_idx)
ct <- concat.samples(mean_celltype, braak_idx)
braak <- concat.samples(braakLabels, braak_idx)
donor <- lapply(donorNames, function(d){
  samples <- unlist(braak_idx[[d]])
  rep(d, length(samples))
  # rep(d, ncol(brainExpr[[d]]))
})
donor <- unlist(donor)
color <- concat.samples(sapply(sampleInfo, function(x)x[, "color_hex_triplet"]))
color <- color[unlist(braak_idx)]
# color[label == "0"] <- "#808080"
gene <- unlist(braakExpr[g, ])
coefficients1 <- sapply(lm_fit1, function(fit) coef(fit)[, g])
coefficients2 <- sapply(lm_fit2, function(fit) coef(fit)[, g])

pdf("EQ1and2_ATCAY_neuronal_marker_concat.pdf", 6, 4)
lapply(names(celltypes), function(c){
  ref <- ct[c, ] # Cell-type expression
  df <- data.frame(ref, gene, factor(braak), donor)
  df$coef1 <- coefficients1[paste0("ct", c), df$donor]
  df$coef2 <- coefficients2[paste0("ct", c), df$donor]
  df$intercept1 <-coefficients1["(Intercept)", df$donor]
  df$intercept2 <-coefficients2["(Intercept)", df$donor]
  df$donor <- factor(df$donor, level = unique(df$donor))
  ggplot(df) +
    geom_point(aes(x=ref, y=gene, shape = braak), color = color, size = 1.5) +
    geom_abline(aes(intercept = intercept1, slope = coef1, linetype = donor), size = 1, color = "black") +
    geom_abline(aes(intercept = intercept2, slope = coef2, linetype = donor), size = 1, color = "red") +
    # geom_smooth(aes(x=ref, y=gene, color = donor), method = lm, se = FALSE, fullrange = TRUE) +
    scale_shape_manual(values=1:length(unique(label))) +
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) + 
    expand_limits(x=0, y=min(df$intercept1, df$intercept2)-2) +
    labs(x = c, y = entrezId2Name(g)) +
    theme
})
dev.off()

########## Use residuals as cell-corrected expression ##########
expr_celltype_corrected <- lapply(donorNames, function(d){
  fit <- lm_fit1[[d]]
  t(fit$residuals)
})
saveRDS(expr_celltype_corrected, file = "resources/expr_corrected_lm_eg_braaksamples.rds")

########## Heatmap Before correction ##########


colPal <- c("blue", "white", "red")
rampcols <- colorRampPalette(colors = colPal, space="Lab")(100)
rampcols <- c(rep(col2hex(colPal[1]), 50), rampcols, rep(col2hex(colPal[3]), 50))

lapply(donorNames[1], function(d){
  idx <- unlist(braak_idx[[d]]) # for samples
  info <- sampleInfo[[d]][idx, ]
  
  # Expression across all Braak samples
  expr <- brainExpr[[d]]
  
  file = paste0("heatmap_celltype_markers_", d, "_darmanis2015.pdf")
  pdf(file, 8, 6)
  t <- t(expr[ct_genes, idx])
  heatmap.2(t, col = rampcols,
            labRow = info$acronym, labCol = entrezId2Name(colnames(t)),
            Rowv=FALSE, Colv=FALSE,
            cexCol = .5, cexRow = .5,
            scale = "none", trace = "none", dendrogram = "none", #key = FALSE,
            RowSideColors = info$color_hex_triplet, ColSideColors = celltypeColors,
            main = paste0("Expression of cell types in ", d, " (uncorrected)"),
            margins = c(5, 5))
  dev.off()
  
  file = paste0("heatmap_celltype_markers2_", d, "_darmanis2015.pdf")
  pdf(file, 8, 6)
  lapply(names(celltypes), function(ct){
    t <- t(expr[celltypes[[ct]], idx])
    heatmap.2(t, col = rampcols,
              labRow = info$acronym, labCol = entrezId2Name(colnames(t)),
              Rowv=FALSE, Colv=FALSE,
              cexCol = .5, cexRow = .5,
              scale = "none", trace = "none", dendrogram = "none", #key = FALSE,
              RowSideColors = info$color_hex_triplet, #ColSideColors = celltypeColors,
              main = paste0("Expression of cell types; ", d, "; ", ct),
              margins = c(5, 5))
  })
  dev.off()
})