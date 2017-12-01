back.transform <- function(x){
  a <- exp(1)^(2*x)
  as.data.frame((a-1)/(a+1))
}