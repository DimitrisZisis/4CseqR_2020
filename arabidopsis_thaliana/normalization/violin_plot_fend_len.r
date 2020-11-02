rm(list=ls()) 

library(ggplot2)
library(plyr)
require(graphics)
require(grDevices)

data1 <- read.csv("ranksort_arab_1_1_nozero_blind_fend_ln.txt",header = F, sep="\t")
data1[1:10,]
names(data1)
logcounts = log10(data1$V4 + 1)
logcounts[1:10]
newdata = cbind(data1, logcounts)
#Use only the fragments with counts in the graphs
countData = newdata[rowSums(newdata[4:5])>0, ]

#Create categories of length based on the values of 5 and 3 fragment end reordering levels
len5 = countData$V8
len3 = countData$V9
head(len5)
str(len5)
str(factor(len5))
str(factor(len3))
flen5<- c("50", "100", "150", "200", "400","600", "800", ">800")[
findInterval(len5 , c(-Inf, 51, 101, 151, 201, 401, 601, 801,  Inf) ) ]
flen5[1:100]
flen3<- c("50", "100", "150", "200", "400","600", "800", ">800")[
  findInterval(len3 , c(-Inf, 51, 101, 151, 201, 401, 601, 801,  Inf) ) ]
flen3[1:100]
countData2=cbind(countData,flen5,flen3)
countData2[1:10,]

#Re-order the levels
countData2$flen3 <- factor(countData2$flen3, levels(countData2$flen3)[c(6,2:5,7,8,1)])
str(countData2)
levels(countData2$flen3)
newdata2=revalue(countData2$flen3, c("'50'"="50","'100'"="100", "'150'"= "150","'200'"="200", "'400'"="400", "'600'"="600","'800'"="800","'>800'"=">800"))
levels(newdata2)


#time = 1

#data1[,time+5][which(data1[,time+5]==-5)]=NA

#data1[,1][which(data1[,2]==0)]="Non blind"
#data1[,1][which(data1[,1]==1)]="Blind"

#xx=factor(data1[,2],levels=c("50","100","150","200","250","300","350",">400"))
xx=newdata2
xx

col=c("#4c72b0","#55a868","#c44e52","#8172b2","#ccb974","#64b5cd")

yy = countData2$logcounts
yy
#alternative with  scale = "width"
p <- ggplot(data=countData2, aes(x=xx, y=yy))+ 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha=0.80, fill=col[1], lty=1, scale = "width") +
  #stat_summary(fun.y="median", geom="point", size=4,color="red")+
  #theme(panel.background=element_blank(),axis.line = element_line(),legend.key=element_rect(fill="white")) +
  theme(panel.background=element_blank(),legend.key=element_rect(fill="white"),panel.border=element_rect(fill=NA)) +
  labs(title="Violin_Plot_flen3'",x="",y = "Log10(Read count + 1)")

p 

pdf1=paste("violin_plot_arab_1_1_flen3_logcounts.pdf",sep="")
ggsave (pdf1, device="pdf", width=7, height=7, dpi=600)



