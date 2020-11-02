newtab<- read.table(file= "arab_2_1_nozero_qnorm_test.txt", header=FALSE, sep="\t")

chrsort=newtab[with(newtab, order(newtab$V1, newtab$V2)), ]
chrsort[1:10, ]
chrsort$V4
logCounts = log10(chrsort$V4 + 1)
rankedCounts = chrsort$V5
qnorm_ranked_counts = chrsort$V6

#density-dot plot 
ggplot(chrsort, aes(x=logCounts, y=rankedCounts) ) +
  geom_point(colour = "red", size = 1)
ggplot(chrsort, aes(x=logCounts, y=qnorm_ranked_counts) ) +
  geom_point(colour = "red", size = 1)


#Compare the normalization by different C 
newtab1<- read.table(file= "arab_1_1_nozero_qnorm_c_bloom.txt", header=FALSE, sep="\t")
newtab2<- read.table(file= "arab_1_1_nozero_qnorm_c_rankit.txt", header=FALSE, sep="\t")
newtab3<- read.table(file= "arab_1_1_nozero_qnorm_c_turkey.txt", header=FALSE, sep="\t")
newtab4<- read.table(file= "arab_1_1_nozero_qnorm_c_vanderWer.txt", header=FALSE, sep="\t")
newtab1[1:10, ]
newtab2[1:10, ]
newtab3[1:10, ]
newtab4[1:10, ]
total_ranks = cbind(newtab1$V1,newtab1$V2,newtab1$V3, newtab1$V4,newtab1$V5,newtab2$V5,newtab3$V5,newtab4$V5)
total_qnorm_ranks = cbind(newtab1$V1,newtab1$V2,newtab1$V3, newtab1$V4,newtab1$V6,newtab2$V6,newtab3$V6,newtab4$V6)
write.table(total_qnorm_ranks, file = "arab_1_1_total_c_qnormranks.txt", append=TRUE , quote = FALSE,  sep = "\t", row.names=FALSE , col.names=FALSE)
