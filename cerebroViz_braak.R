#CerebroViz
setwd("C:/Users/dkeo/surfdrive/Parkinson/Images/CerebroViz")
library("cerebroViz")


data(regionMap)
regionMap

#braak regions
b1 <- "MED"
b2 <- "PON"
b3 <- "SN"
b4 <- "AMY"
b5 <- c("CNG", "TL")
b6 <- c("FL", "PL")
b<- c(b1, b2, b3, b4, b5, b6)

tab <- matrix(c(1:length(b)), length(b), 1)
rownames(tab) <- b

#AHBA region colors
ontology <- read.csv("../../../ABA_human_processed/Ontology_edited.csv")
ahbaAcronyms <- c("MY", "Pons", "SN", "Amg", "CgG", "TL", "FL", "PL")
colors <-  sapply(ontology$color_hex_triplet[match(ahbaAcronyms, ontology$acronym)], function (x)paste0("#", x))


cerebroViz(tab, regLabel = TRUE, palette = colors, legend = FALSE)
