library(ggplot2)
library(plyr)
require(graphics)
require(grDevices)
library(Biostrings)
require(graphics)
#Create a table with the total normalized ranks from each sample 

data1 <- read.csv("ranksort_arab_1_1_qnorm.txt",header = F, sep=" ")
data1[1:10,]
data2 <- read.csv("ranksort_arab_2_1_qnorm.txt",header = F, sep=" ")
data2[1:10,]
data3 <- read.csv("ranksort_arab_3_1_qnorm.txt",header = F, sep=" ")
data4 <- read.csv("ranksort_arab_1_4_qnorm.txt",header = F, sep=" ")
data1[1:10,]
data5 <- read.csv("ranksort_arab_2_4_qnorm.txt",header = F, sep=" ")
data2[1:10,]
data6 <- read.csv("ranksort_arab_3_4_qnorm.txt",header = F, sep=" ")
table = cbind(as.character(data1$V1),data1$V2,data1$V3,data1$V5,data2$V5,data3$V5,data4$V5,data5$V5,data6$V5)
colnames(table) <- c("chr","start","end","arab_1_1","arab_2_1","arab_3_1","arab_1_4","arab_2_4","arab_3_4")

table[1:10,]
write.table(table, file = "totalnorm_all_arab_qnorm.txt" , quote = FALSE,  sep = "\t", row.names=FALSE , col.names=T)

newtab<- read.table(file= "totalnorm_all_arab_qnorm.txt", header=T, sep="\t")
newtab[1:10,]

data2 <- read.csv("G:/4CseqR-master/4CseqR/Library_Fragments/find_blind/AGATCT_blindNonBlind_frend_len_zero.txt",header = F, sep="\t")
data2[1:10,]

df1= newtab

df2=data2[1:3]

#add the missing zero fragments to the normalized results. 
index <- match(paste0(df2$V2, df2$V3), paste0(df1$start, df1$end))
df2$V4 <- df1$arab_1_1[index]
df2$V5 <- df1$arab_2_1[index]
df2$V6 <- df1$arab_3_1[index]
df2$V7 <- df1$arab_1_4[index]
df2$V8 <- df1$arab_2_4[index]
df2$V9 <- df1$arab_3_4[index]
df2[is.na(df2)] <- 0
df2

#Write the results in a new file
write.table(df2, file = "totalnorm_all_arab_qnorm_zero.txt" , quote = FALSE,  sep = " ", row.names=FALSE , col.names=FALSE)
