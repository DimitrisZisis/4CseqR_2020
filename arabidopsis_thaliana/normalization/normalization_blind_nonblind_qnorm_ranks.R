gc()
rm(list=ls()) 
library(Biostrings)
library(Biostrings)
#setwd("/media/dimitris/3AEC27A3EC275881/4CseqR-master/4CseqR/Library_Fragments")

newdata <- read.csv(file="/4CseqR-master/4CseqR/Library_Fragments/arab_1_1_salmon.csv", header=T, sep="\t")
newdata[1:10,]

data2 <- read.csv("G:/4CseqR-master/4CseqR/Library_Fragments/find_blind/AGATCT_blindNonBlind_frend_len_zero.txt",header = F, sep="\t")
data2[1:20,]
combdata = cbind(newdata[1:7],data2[4:6])

#select only blind 
blinddata=subset(combdata, V4=="Blind") 
blindData1=blinddata[1:8][complete.cases(blinddata[1:8]), ]
colnames(blindData1) <- c("Name","start","end","Length","EffectiveLength","TPM","NumReads","Blind")

#Select the rest of nonblind and remove NA rows
countData1=combdata[complete.cases(combdata), ]
colnames(countData1) <- c("Name","start","end","Length","EffectiveLength","TPM","NumReads","Blind","5'frag_len","3'frag_len")

#remove zeros
countData = countData1[rowSums(countData1[6:7])>0, ]
countData[1:10,]
blindData = blindData1[rowSums(blindData1[6:7])>0, ]
blindData[1:10,]

blindcount=blindData$NumReads
blindChr=blindData$Name
blindend=blindData$end
blindstart=blindData$start
fblind =as.character(blindData$Blind)

newdatacount=countData$NumReads
newdataChr=countData$Name
newdataend=countData$end
newdatastart=countData$start
fnonblind =as.character(countData$Blind)
newdata5flen=countData$`5'frag_len`
newdata3flen=countData$`3'frag_len`

#calculate the lenght and the middle point
lenght = countData$Length
lenght_blind = blindData$Length
middle = (newdataend + newdatastart)/2

ldist = countData$`5'frag_len`
rdist =countData$`3'frag_len`

#create categories of distance and lenght using fingInterval
fldist<-c("50", "100", "150", "200","250", "300", "350", "400",">400")[
  findInterval(ldist , c(-Inf, 51, 101, 151, 201, 251, 301, 351, 401, Inf) ) ]
fldist[1:100]
frdist <- c("50", "100", "150", "200","250", "300", "350", "400",">400")[
  findInterval(rdist , c(-Inf, 51, 101, 151, 201, 251, 301, 351, 401, Inf) ) ]
frdist[1:100]
flenght <- c("50", "100", "150", "200", "400","600", "800", ">800")[
  findInterval(lenght , c(-Inf, 51, 101, 151, 201, 401, 601, 801,  Inf) ) ]
flenght[1:100]
flenght_blind <- c("50", "100", "150", "200", "400","600", "800", ">800")[
  findInterval(lenght_blind , c(-Inf, 51, 101, 151, 201, 401, 601, 801,  Inf) ) ]
flenght_blind[1:100]


#1st way create categories
kat3 <- cbind(fnonblind,fldist,frdist)
kat3[1:10,]

kat_blind <- cbind(fblind)
kat_blind[1:10,]

#Normalization for nonblind combined categories
V3=paste(kat3[,1],kat3[,2],kat3[,3])
at3=sort(unique(V3))
nlev=length(at3)
ach3=NULL
arst3=NULL
aren3=NULL
ac3=NULL
ac3h=NULL
rank3=NULL
rankn=NULL
c=0.5
#Normalization for blind only
V=paste(fblind)
at=sort((V))
nlevb=length(at)
ach=NULL
arst=NULL
aren=NULL
arbl=NULL
ac=NULL
rank=NULL
rankm=NULL
c=0.5
for (kk in 1:nlev){
  ach3[[kk]]=newdataChr[which(V3==at3[kk])]
  arst3[[kk]]=newdatastart[which(V3==at3[kk])]
  aren3[[kk]]=newdataend[which(V3==at3[kk])]
  ac3[[kk]]=newdatacount[which(V3==at3[kk])]
  v3=ac3[[kk]]
  v3[which(v3==0)]=NA
  #rank only the numeric values withoun ranking NA's.Better results in ranking 
  #rank3[[kk]]=rank(v3,na.last="keep")
  #rank3[[kk]]=rank3[[kk]]/sum(!is.na(rank3[[kk]]))
  #rank3[[kk]][is.na(v3)] <- 0
  
  rank3[[kk]]=rank(v3,na.last="keep")
  rank3[[kk]]=( rank3[[kk]] - c ) / (sum(!is.na(rank3[[kk]])) - 2*c + 1  )
  rank3[[kk]][is.na(v3)] <- 0
  rankn[[kk]]= qnorm(rank3[[kk]], mean = 0, sd = 1, lower.tail = TRUE)

}

for (i in 1:nlevb){
  ach[[i]]=blindChr[which(V==at[i])]
  arst[[i]]=blindstart[which(V==at[i])]
  aren[[i]]=blindend[which(V==at[i])]
  ac[[i]]=blindcount[which(V==at[i])]
  v=ac[[i]]
  v[which(v==0)]=NA
  #rank only the numeric values withoun ranking NA's.Better results in ranking 
  #rank[[i]]=rank(v,na.last="keep")
  #rank[[i]]=rank[[i]]/sum(!is.na(rank[[i]]))
  #rank[[i]][is.na(v)] <- 0
  rank[[i]]=rank(v3,na.last="keep")
  rank[[i]]=( rank[[i]] - c ) / (sum(!is.na(rank[[i]])) - 2*c + 1  )
  rank[[i]][is.na(v3)] <- 0
  rankm[[i]]= qnorm(rank[[i]], mean = 0, sd = 1, lower.tail = TRUE)
}


for (kk in 1:nlev) {
  U3=cbind(as.character(ach3[[kk]]),arst3[[kk]],aren3[[kk]],ac3[[kk]],rankn[[kk]],at3[[kk]])
  #print(U3)
  write.table(U3, file = "arab_1_1_nonblind_qnorm.txt" ,append=TRUE, quote = FALSE,  sep = "\t", row.names=FALSE , col.names=FALSE)
}

Ub=cbind(as.character(ach[[i]]),arst[[i]],aren[[i]],ac[[i]],rankm[[i]],at[[i]])
print(Ub)
write.table(Ub, file = "arab_1_1_blind_qnorm.txt", append=TRUE , quote = FALSE,  sep = "\t", row.names=FALSE , col.names=FALSE)

#Sort the table and produce total normalized counts 
newtab1<- read.table(file= "arab_1_1_nonblind_qnorm.txt", header=FALSE, sep="\t")
newtab2<- read.table(file= "arab_1_1_blind_qnorm.txt", header=FALSE, sep="\t")
x=rbind(newtab1,newtab2)
chrsort=x[with(x, order(x$V1, x$V2)), ]
write.table(chrsort, file = "ranksort_arab_1_1_all_qnorm_blind.txt" , quote = FALSE,  sep = " ", row.names=FALSE , col.names=FALSE)
nrow(chrsort)
