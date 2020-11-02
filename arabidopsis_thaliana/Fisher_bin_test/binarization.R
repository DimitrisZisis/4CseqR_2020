newdata <- read.csv(file="/4CseqR-master/4CseqR/Library_Fragments/arab_3_4_salmon.csv", header=T, sep="\t")

newdata$NumReads
bin_param <-1
newdata$binary <- ifelse(newdata$NumReads>bin_param, 1, 0)
newdata[1:10,]
write.table(newdata, file = "arab_3_4_salmon_binary_p_1.txt", append=TRUE , quote = FALSE,  sep = "\t", row.names=FALSE , col.names=FALSE)

