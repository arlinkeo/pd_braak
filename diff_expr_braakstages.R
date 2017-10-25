# Differential expression in braak regions of PD

setwd("C:/Users/dkeo/surfdrive/Parkinson")
library("metap")

source("PD/base_script.R")
load("../ABA_Rdata/BrainExpr.RData")
nDonors <- length(donorNames)
load("resources/braakStages.RData")

#Differential expression for each braak stage compared to non-braak stage
diffGenesList <- lapply(braakNames, function(bs){
  braakInfo <- braakStages[[bs]]
  # Differential expression for each donor
  diffGenesPerBrain <- lapply(donorNames, function(d){
    expr <- brainExpr[[d]]
    nGenes <- nrow(expr)
    braakVec <- as.logical(braakInfo[[d]])
    nonBraakVec <- as.logical(nonBraak[[d]])
    # T-test for each gene
    pvalues <- as.data.frame(t(apply(expr, 1, function(g){
      exprB <- g[braakVec]
      exprNB <- g[nonBraakVec]
      test2tail <- t.test(exprB, exprNB, alternative = "two.sided")
      pval2tail <- test2tail$p.value
      meanB <- test2tail$estimate[1]
      meanNB <- test2tail$estimate[2]
      names(meanB) <- 'meanB'
      names(meanNB) <- 'meanNB'
      varB <- var(unlist(exprB))
      varNB <- var(unlist(exprNB))
      confidence95 <- test2tail$conf.int
      # fChange <- log2(2^test2tail$estimate[1] / 2^test2tail$estimate[2]) # estimate is mean
      # names(fChange) <- 'fold-change'
      c('pvalue' = pval2tail, meanB, 'varB' = varB, meanNB, 'varNB' = varNB, 'lower95' = confidence95[1], 'upper95' = confidence95[2])
    })))
    pvalues$'benjamini_hochberg' <- p.adjust(pvalues$'pvalue' , method = "BH", n = nGenes) #corrected p
    pvalues$'bonferroni' <- p.adjust(pvalues$'pvalue' , method = "bonferroni", n = nGenes) #corrected p
    pvalues
  })
})
save(diffGenesList, file = "resources/diffGenesBraak.RData")

# Number of diff. expr. genes in each braak stage en donor after multiple testing correction
nDiffgenes <- function(x){
  sapply(diffGenesList, function(bs){
    sapply(bs, function(d){
      sum(d[[x]] < 0.05)
    })
  })
}

nDiffgenes("bonferroni")
nDiffgenes("benjamini_hochberg")

nDonors <- 6

# convert lists to gene x donor table for each braak stage
dg.table <- function(x){
  lapply(braakNames, function(bs){
    brainList <- diffGenesList[[bs]]
    genes <- rownames(brainList[[1]])
    tab <- sapply(brainList, function(t){v <- t[ , x]})
    rownames(tab) <- genes
    tab
  })
}

dgTable1 <- dg.table("benjamini_hochberg")
dgTable2 <- dg.table("bonferroni")

# Diff. expressed in all donors
diff.all <- function(list, th) {
  braakTab <- sapply(list, function(tab){
    apply(tab < th, 1, sum) == nDonors
  })
  rownames(braakTab) <- rownames(list[[1]])
  braakTab
}

dgAll1 <- diff.all(dgTable1, 0.05)
dgAll2 <- diff.all(dgTable2, 0.05)
apply(dgAll1, 2, sum)
apply(dgAll2, 2, sum)

#Combine p-values
combine.p <- function(list){
  sapply(list, function(tab){
    apply(tab, 1, function(p){
      sumlog(p)$p#pchisq(-2*sum(log(p)), df = length(p)*2, lower.tail = FALSE)
    })
  })
}

fisherP1 <- combine.p(dgTable1)
fisherP2 <- combine.p(dgTable2)
apply(fisherP1 < 0.05, 2, sum)
apply(fisherP2 < 0.05, 2, sum)

#Get list of sorted significant genes
signif.genes <- function(tab){
  lapply(braakNames, function(bs){
    (head(tab[, bs])) < 0.05
  })
}