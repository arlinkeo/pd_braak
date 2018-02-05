# Differential co-expression

setwd("C:/Users/dkeo/surfdrive/pd_braak")
source("PD/base_script.R")

# Load co-expression matrices
b1_coexpr <- readRDS("resources/avgCoexpr_braak1.rds")
b6_coexpr <- readRDS("resources/avgCoexpr_braak6.rds")

# Differential co-expression
diff_coexpr_b1b6 <- b1_coexpr - b6_coexpr
