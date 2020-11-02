library(ggplot2)
library(Biostrings)
#Create a table with the total normalized ranks from each sample 

data1 <- read.csv("ESC_pcdhb19_1_salmon_binary_p_1.txt",header = F, sep="\t")
data1[1:10,]
data2 <- read.csv("ESC_pcdhb19_2_salmon_binary_p_1.txt",header = F, sep="\t")
data2[1:10,]
data3 <- read.csv("ESC_pcdhb19_3_salmon_binary_p_1.txt",header = F, sep="\t")
data3[1:10,]
data4 <- read.csv("iPSC_pcdhb19_1_salmon_binary_p_1.txt",header = F, sep="\t")
data4[1:10,]
data5 <- read.csv("iPSC_pcdhb19_2_salmon_binary_p_1.txt",header = F, sep="\t")
data5[1:10,]
data6 <- read.csv("iPSC_pcdhb19_3_salmon_binary_p_1.txt",header = F, sep="\t")
data6[1:10,]
table = cbind(as.character(data1$V1),data1$V2,data1$V3,data1$V9,data2$V8,data3$V8,data4$V8,data5$V8,data6$V9)
colnames(table) <- c("chr","start","end","ESC_pcdhb19_1","ESC_pcdhb19_2","ESC_pcdhb19_3","iPSC_pcdhb19_1","iPSC_pcdhb19_2","iPSC_pcdhb19_3")

table[1:10,]
write.table(table, file = "totalnorm_all_mouse_binary.txt" , quote = FALSE,  sep = "\t", row.names=FALSE , col.names=T)


newtab<- read.table(file= "totalnorm_all_mouse_binary.txt", header=T, sep="\t")
newtab[1:10,]

data2 <- read.csv("G:/4CseqR-master/4CseqR/Library_Fragments/find_blind/mm10_AAGCTT_blindNonBlind_frend_len_zero.txt",header = F, sep="\t")
data2[1:10,]

df1= newtab

df2=data2[1:3]

#add the missing zero fragments to the normalized results. 
index <- match(paste0(df2$V2, df2$V3), paste0(df1$start, df1$end))
df2$V4 <- df1$ESC_pcdhb19_1[index]
df2$V5 <- df1$ESC_pcdhb19_2[index]
df2$V6 <- df1$ESC_pcdhb19_3[index]
df2$V7 <- df1$iPSC_pcdhb19_1[index]
df2$V8 <- df1$iPSC_pcdhb19_2[index]
df2$V9 <- df1$iPSC_pcdhb19_3[index]
df2[is.na(df2)] <- 0
df2

#Write the results in a new file
write.table(df2, file = "totalnorm_all_mm10_binary_zero.txt" , quote = FALSE,  sep = " ", row.names=FALSE , col.names=FALSE)
