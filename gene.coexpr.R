# Gene coexpression
gene_coexpr <- function(exprList){
  
  # Correlations only without p-values
  sapply(exprList, function(expr){
    cor(t(expr), method = "pearson") #List with correlation-, size-, and p-value-matrix
  }, simplify = FALSE)
  
  # Correlations with p-values
  # require("Hmisc")
  # sapply(exprList, function(expr){
  #   m <- rcorr(t(expr), type = "pearson") #List with correlation-, size-, and p-value-matrix
  #   m$n <- m$n[1, 1]
  #   m
  # }, simplify = FALSE)
  
}