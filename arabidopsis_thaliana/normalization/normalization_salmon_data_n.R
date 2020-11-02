gc()
rm(list=ls()) 
library(Biostrings)
library(Biostrings)
#setwd("/media/dimitris/3AEC27A3EC275881/4CseqR-master/4CseqR/Library_Fragments")

newdata <- read.csv(file="/4CseqR-master/4CseqR/Library_Fragments/arab_1_1_salmon.csv", header=T, sep="\t")
newdata[1:10,]

data2 <- read.csv("G:/4CseqR-master/4CseqR/Library_Fragments/find_blind/AGATCT_blindNonBlind_frend_len_zero.txt",header = F, sep="\t")
data2[1:10,]
combdata = cbind(newdata,data2[4:6])
#remove NA rows
countData1=combdata[complete.cases(combdata), ]
colnames(countData1) <- c("Name","start","end","Length","EffectiveLength","TPM","NumReads","Blind","5'frag_len","3'frag_len")

#remove zeros
countData = countData1[rowSums(countData1[6:7])>0, ]
countData[1:10,]
newdatacount=countData$NumReads
newdataChr=countData$Name
newdataend=countData$end
newdatastart=countData$start
fblind =as.character(countData$Blind)
newdata5flen=countData$`5'frag_len`
newdata3flen=countData$`3'frag_len`
#calculate the lenght and the middle point
lenght = countData$Length
lenght
middle = (newdataend + newdatastart)/2
middle[1:10]
ldist = countData$`5'frag_len`
rdist =countData$`3'frag_len`
#distributions of cov5 and cov3 are similar-fragment coverage total estimation
tot = newdatacount
tot[1:100]
tot = log10(tot + 1)

tot[1:100]
#create categories of distance and lenght 
head(lenght)
str(lenght)
str(factor(lenght))

#2nd way using fingInterval
fldist<-c("50", "100", "150", "200","250", "300", "350", "400",">400")[
  findInterval(ldist , c(-Inf, 51, 101, 151, 201, 251, 301, 351, 401, Inf) ) ]
fldist[1:100]
frdist <- c("50", "100", "150", "200","250", "300", "350", "400",">400")[
  findInterval(rdist , c(-Inf, 51, 101, 151, 201, 251, 301, 351, 401, Inf) ) ]
frdist[1:100]
flenght <- c("50", "100", "150", "200", "400","600", "800", ">800")[
  findInterval(lenght , c(-Inf, 51, 101, 151, 201, 401, 601, 801,  Inf) ) ]
flenght[1:100]

#categories for Normalization-combine all factors-3 ways 
f1 <- factor(fldist)
f2 <- factor(frdist)
f3 <- factor(flenght)
f4 <- factor(fblind)

#1st way
#kat5 <- cbind(fblind,fldist,flenght)
#kat5[1:10,]

kat3 <- cbind(fblind,frdist,fldist)
kat3[1:10,]

#Normalization for 3 prime ends
V3=paste(kat3[,1],kat3[,2],kat3[,3])
at3=sort(unique(V3))
at3
nlev=length(at3)
ach3=NULL
arst3=NULL
aren3=NULL
arbl3=NULL
ar5fln=NULL
ar3fln=NULL
ac3=NULL
ac3h=NULL
rank3=NULL
rankn=NULL
rank=NULL
c=0
for (kk in 1:nlev){
  ach3[[kk]]=newdataChr[which(V3==at3[kk])]
  arst3[[kk]]=newdatastart[which(V3==at3[kk])]
  aren3[[kk]]=newdataend[which(V3==at3[kk])]
  ac3[[kk]]=newdatacount[which(V3==at3[kk])]
  #arbl3[[kk]]=newdatablind[which(V3==at3[kk])]
  #ar5fln[[kk]]=newdata5flen[which(V3==at3[kk])]
  #ar3fln[[kk]]=newdata3flen[which(V3==at3[kk])]
  v3=ac3[[kk]]
  v3[which(v3==0)]=NA
  #rank only the numeric values withoun ranking NA's.Better results in ranking 
  #rank3[[kk]]=rank(v3,na.last="keep")
  #rank3[[kk]]=rank3[[kk]]/sum(!is.na(rank3[[kk]]))
  #rank3[[kk]][is.na(v3)] <- 0
  
  rank[[kk]]=rank(v3,na.last="keep")
  rank[[kk]]=( rank[[kk]] - c ) / (sum(!is.na(rank[[kk]])) - 2*c + 1  )
  rank[[kk]][is.na(v3)] <- 0
  rankn[[kk]]= qnorm(rank[[kk]], mean = 0, sd = 1, lower.tail = TRUE)
  #rank3[[kk]]=rank(v3,na.last=FALSE)
  #rank3[[kk]]=rank3[[kk]]/length(v3)
  #rank3[[kk]][is.na(v3)] <- 0
  # if ((all(is.na(ac3[[kk]]))) - make all 0
}


for (kk in 1:nlev) {
  U3=cbind(as.character(ach3[[kk]]),arst3[[kk]],aren3[[kk]],ac3[[kk]],rank[[kk]], rankn[[kk]],at3[[kk]])
  #print(U3)
  write.table(U3, file = "arab_1_1_nozero_qnorm_c_vanderWer.txt", append=TRUE , quote = FALSE,  sep = "\t", row.names=FALSE , col.names=FALSE)
}

#Sort the table and produce total normalized counts 
newtab<- read.table(file= "arab_3_1_nozero_qnorm_test.txt", header=FALSE, sep="\t")
newtab[1:10, ]
chrsort=newtab[with(newtab, order(newtab$V1, newtab$V2)), ]
#totnormal = chrsort$V5
#t=cbind(as.character(chrsort$V1),chrsort$V2,chrsort$V3,totnormal)
#totnormal[1:10]
write.table(chrsort, file = "ranksort_arab_3_1_qnorm.txt" , quote = FALSE,  sep = " ", row.names=FALSE , col.names=FALSE)
#write.table(t, file = "totalnorm_iPSC_pcdhb19_3_salmon_nozero.txt", append=TRUE , quote = FALSE,  sep = "\t", row.names=FALSE , col.names=FALSE)
