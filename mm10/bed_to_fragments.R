gc()
rm(list=ls()) 
library(Biostrings)
library(dplyr)
data <- read.csv(file="/media/dimitris/3AEC27A3EC275881/Mouse_fragments_mapping/restsite_AAGCTT_test.bed", header=FALSE, sep="\t")
newdata <- data %>% 
  transmute(
    chr = V1,
    start = lag(V2, default = 1),
    end = V3)
newdata
