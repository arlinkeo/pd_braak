# Differential expression between regions in UKBEC
# First download all data files from http://www.braineac.org/ and store in UKBEC directory

# Set UKBEC data directory
ukbec_dir <- "C:/Users/dkeo/surfdrive/UKBEC"

##############################################################################################
# Map affy ID to entrez IDs
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version = 92)
affyID <- expr.map$exprID
system.time({
  ukbecGeneID <- getBM(c('ensembl_gene_id', 'entrezgene', 'hgnc_symbol', "affy_huex_1_0_st_v2"), 
                       filters=c("affy_huex_1_0_st_v2"), mart=ensembl, values=affyID)
})
ukbecGeneID$entrezgene <- as.character(ukbecGeneID$entrezgene)
ukbecGeneID$affy_huex_1_0_st_v2 <- as.character(ukbecGeneID$affy_huex_1_0_st_v2)
saveRDS(ukbecGeneID, file = "output/ukbecGeneID.rds")

##############################################################################################
# Probes to genes

# Read expression data for all brain regions
regions <- c("CRBL", "FCTX", "HIPP","MEDU", "OCTX", "PUTM", "SNIG", "TCTX", "THAL", "WHMT")
regionExpr <- sapply(regions, function(r){
  fName <- paste0(ukbec_dir, "/expr_", r, ".txt")
  x=read.csv(fName, header = TRUE, sep = " ", row.names = 1)
}, simplify = FALSE)
exprConcat <- Reduce(cbind, regionExpr) #Concatenate brain regions

# Mapping probes to genes
probes <- unique(ukbecGeneID$affy_huex_1_0_st_v2)
genes <- ukbecGeneID$entrezgene[match(probes, ukbecGeneID$affy_huex_1_0_st_v2)] # entrez IDs
nas <- is.na(genes)
probes <- probes[!nas]
genes <- genes[!nas]
exprConcat <- exprConcat[probes,]
probeSelection <- collapseRows(exprConcat, rowGroup = genes, rowID = probes, 
                                method = "maxRowVariance", connectivityBasedCollapsing = TRUE)

# Select probes for each data matrix
probes <- probes[probeSelection$selectedRow]
genes <- genes[probeSelection$selectedRow]
regionExpr <- lapply(regionExpr, function(x) {
  x <- x[probes, ]
  rownames(x) <- genes
  x
})
save(regionExpr, file = paste0(ukbec_dir, "/regionExpr.RData"))

##############################################################################################
# T-test in UKBEC

roi <- c('1' = "MEDU", '3' = "SNIG", '5' = "TCTX", '6' = "FCTX")
regionpairs <- combn(roi, 2)
colnames(regionpairs) <- apply(regionpairs, 2, function(x) paste0(x[1], "-",  x[2]))
ttest <- alply(regionpairs, 2, function(x){
  df1 <- regionExpr[[x[1]]]
  df2 <- regionExpr[[x[2]]]
  t.test.table(df1,df2)
}, .dims = TRUE)
ttest <- simplify2array(ttest) # 3D array: genes x measures x region pairs

# Number of diff. genes
apply(ttest, c(3), function(x){
  sum(x[, "BH"] < 0.05 & abs(x[, "meanDiff"]) > 1)
})
saveRDS(ttest, file = "output/ttest_ukbec.rds")

# Number of diff. genes
sum(abs(ttest[,"meanDiff",]) > 1 & ttest[,"BH",] < 0.05)

##############################################################################################
# Overlap with BRGs
genes_ukbec <- rownames(regionExpr$CRBL)

# Intersection with entrez IDs in UKBEC
bg_ukbec <- lapply(bg, function(g){
  intersect(g, genes_ukbec)
})

# BRGs differentially expressed in UKBEC
degs_num_ukbec <- apply(ttest, c(3), function(x){
  sapply(bg_ukbec, function(g){
    sum(x[g, "BH"] < 0.05 & abs(x[g, "meanDiff"]) > 1)
  })
})
degs_percentage_ukbec <- apply(degs_num_ukbec, 2, function(x) round(x/sapply(bg_ukbec,length)*100, digits = 2))
degs_num_ukbec
degs_percentage_ukbec

##############################################################################################
# Plotting functions

# prepare.data <- function(g){ # prepare ggplot dataframe for single genes
#   expr <- sapply(roi, function(r)  unlist(regionExpr[[r]][g, ] ))
#   df <- melt(expr)
#   colnames(df) <- c("sample", "region", "expr")
#   df$region <- paste0("R", df$region)
#   df$region <- factor(df$region, levels = sort(unique(df$region)))
#   df
# }

names(braakColors) <- gsub("braak", "R", names(braakColors))
theme <- theme(panel.background = element_blank(), panel.grid = element_blank(), 
               axis.line = element_line(colour = "black"))

box.plot <- function(df, title){
  ggplot(df) + 
    geom_boxplot(aes(y = expr, x = region, fill = region)) +
    labs(x = "Brain region", y = bquote("Expression ("*log[2]*"-transformed)")) +
    ggtitle(title) +
    scale_x_discrete(expand=c(0.2,0)) +
    scale_fill_manual(values = braakColors, guide = FALSE) +
    theme
}

# plot.pdf <- function(name, genes){ # For plots of single genes
#   pdf(name, 2, 3)
#   lapply(genes, function(g){
#     title <- paste0(entrezId2Name(g))
#     df <- prepare.data(g)
#     p <- box.plot(df, title)
#     print(p)
#   })
#   dev.off()
# }

##############################################################################################
# Box plots

# Mean expression across Braak genes within regions
df <- simplify2array(lapply(roi, function(r){
  sapply(bg_ukbec, function(g){
    expr <- regionExpr[[r]][g, ]
    apply(expr, 2, mean)
  })
}))
df <- melt(df)
colnames(df) <- c("sample", "dir", "region", "expr")
df$region <- paste0("R", df$region)
df$region <- factor(df$region, levels = unique(df$region))

pdf("output/boxplot_UKBEC.pdf", 2.5, 4)
box.plot(df, "BRGs in UKBEC") +
  facet_grid(.~dir, scales = 'free', space = 'free', switch = "y") +
  scale_y_continuous(limits = c(y_min, y_max))
dev.off()

########## Heatmap expression of BRGs ##########

expr <- sapply(roi, function(r){
    e <- regionExpr[[r]][unlist(bg_ukbec), ]
    apply(e, 1, function(x) mean(x, na.rm = TRUE))
})
colnames(expr) <- paste0("R", colnames(expr))
expr <- t(scale(t(expr))) # expr. is scaled across samples

pdf("output/heatmap_expr_BRGs_UKBEC_col1.2.pdf", 2.7, 10)
Heatmap(expr, name = 'Z-Score\nexpression',
        col = col_fun,
        row_split = rep(names(lengths(bg_ukbec)), lengths(bg_ukbec)),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        column_names_gp = gpar(fontsize = 10),
        column_names_rot = 0,
        column_names_centered = TRUE,
        row_title_rot = 0,
        width = unit(ncol(expr), "lines"),
        height = unit(nrow(expr)*.05, "lines")
)
dev.off()

# #boxplot of PD genes
# plot.pdf("boxplot_UKBEC_PD_variant_genes.pdf", 
#          name2EntrezId(c("SNCA", "ZNF184", "BAP1", "SH3GL2", "ELOVL7", "SCARB2")))
