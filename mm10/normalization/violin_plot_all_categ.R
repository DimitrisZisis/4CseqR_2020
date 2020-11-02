library(ggplot2)
library(plyr)
require(graphics)
require(grDevices)

data1 <- read.csv("ranksort_ESC_pcdhb19_1_all_qnorm_blind.txt",header = F, sep=" ")
data1[1:10,]
names(data1)
logcounts = log10(data1$V4 + 1)
logcounts[1:10]
newdata1 = cbind(data1, logcounts)
#Use only the fragments with counts in the graphs
#countData = newdata[rowSums(newdata[4:5])>0, ]
#Re-order the levels
newdata=subset(newdata1, V6=="NonBlind") 

newdata$V7 <- factor(newdata$V7, levels(newdata$V7)[c(10,3:9,2)])
str(newdata)
levels(newdata$V7)
newdata_5end=revalue(newdata$V7, c("'50'"="50","'100'"="100", "'150'"= "150","'200'"="200","'250'"="250","'300'"="300", "'350'"="350", "'400'"="400","'>400'"=">400"))
levels(newdata_5end)

newdata$V8 <- factor(newdata$V8, levels(newdata$V8)[c(10,3:9,2)])
str(newdata)
levels(newdata$V8)
newdata_3end=revalue(newdata$V8, c("'50'"="50","'100'"="100", "'150'"= "150","'200'"="200","'250'"="250","'300'"="300", "'350'"="350", "'400'"="400","'>400'"=">400"))
levels(newdata_3end)

#violin plots for each category 

xx=newdata_5end
col=c("#4c72b0","#55a868","#c44e52","#8172b2","#ccb974","#64b5cd")
yy = newdata$logcounts
#alternative with  scale = "width"
  p <- ggplot(data=newdata, aes(x=xx, y=yy))+ 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha=0.80, fill=col[1], lty=1, scale = "width") +
  theme(panel.background=element_blank(),legend.key=element_rect(fill="white"),panel.border=element_rect(fill=NA)) +
  labs(title="violin_Plot_flen5'",x="",y = "Log10(Read count + 1)")

p 

pdf1=paste("violin_plot_ESC_pcdhb19_1_flen5_all_logcounts.pdf",sep="")
ggsave (pdf1, device="pdf", width=7, height=7, dpi=600)


#Normalized_ranks-qnorm
xx=newdata_5end
col=c("#4c72b0","#55a868","#c44e52","#8172b2","#ccb974","#64b5cd")
yy = newdata$V5
#alternative with  scale = "width"
p <- ggplot(data=newdata, aes(x=xx, y=yy))+ 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha=0.80, fill=col[1], lty=1, scale = "width") +
  theme(panel.background=element_blank(),legend.key=element_rect(fill="white"),panel.border=element_rect(fill=NA)) +
  labs(title="violin_plot_flen5'",x="",y = "QNorm_ranks")

p 

pdf1=paste("violin_plot_ESC_pcdhb19_1_flen5_qnorm.pdf",sep="")
ggsave (pdf1, device="pdf", width=7, height=7, dpi=600)


xx=newdata_3end
col=c("#4c72b0","#55a868","#c44e52","#8172b2","#ccb974","#64b5cd")
yy = newdata$logcounts
#alternative with  scale = "width"
p <- ggplot(data=newdata, aes(x=xx, y=yy))+ 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha=0.80, fill=col[1], lty=1, scale = "width") +
  theme(panel.background=element_blank(),legend.key=element_rect(fill="white"),panel.border=element_rect(fill=NA)) +
  labs(title="violin_plot_flen3'",x="",y = "Log10(Read count + 1)")

p 

pdf1=paste("violin_plot_ESC_pcdhb19_1_flen3_all_logcounts.pdf",sep="")
ggsave (pdf1, device="pdf", width=7, height=7, dpi=600)


#Normalized_ranks
xx=newdata_3end
col=c("#4c72b0","#55a868","#c44e52","#8172b2","#ccb974","#64b5cd")
yy = newdata$V5
#alternative with  scale = "width"
p <- ggplot(data=newdata, aes(x=xx, y=yy))+ 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha=0.80, fill=col[1], lty=1, scale = "width") +
  theme(panel.background=element_blank(),legend.key=element_rect(fill="white"),panel.border=element_rect(fill=NA)) +
  labs(title="violin_plot_flen3'",x="",y = "QNorm_ranks")

p 

pdf1=paste("violin_plot_ESC_pcdhb19_1_flen3_qnorm.pdf",sep="")
ggsave (pdf1, device="pdf", width=7, height=7, dpi=600)

#Blind non blind boxplots

xx=factor(newdata1$V6,levels=c("NonBlind","Blind"))
col=c("#4c72b0","#55a868","#c44e52","#8172b2","#ccb974","#64b5cd")
yy = newdata1$V5
#alternative with  scale = "width"
p <- ggplot(data=newdata1, aes(x=xx, y=yy))+ 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha=0.80, fill=col[1], lty=1, scale = "width") +
  theme(panel.background=element_blank(),legend.key=element_rect(fill="white"),panel.border=element_rect(fill=NA)) +
  labs(title="violin_plot blind/nonblind",x="",y = "QNorm_ranks")

p 
p + coord_trans(x= "log10")
pdf1=paste("violin_plot_ESC_pcdhb19_1_blind_qnorm.pdf",sep="")
ggsave (pdf1, device="pdf", width=7, height=7, dpi=600)

xx=factor(newdata1$V6,levels=c("NonBlind","Blind"))
col=c("#4c72b0","#55a868","#c44e52","#8172b2","#ccb974","#64b5cd")
yy = newdata1$logcounts
#alternative with  scale = "width"
p <- ggplot(data=newdata1, aes(x=xx, y=yy))+ 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha=0.80, fill=col[1], lty=1, scale = "width") +
  theme(panel.background=element_blank(),legend.key=element_rect(fill="white"),panel.border=element_rect(fill=NA)) +
  labs(title="violin_plot blind/nonblind",x="",y = "Log10(Read count + 1)")

p 
p + coord_trans(x= "log10")
pdf1=paste("violin_plot_ESC_pcdhb19_blind_logcounts.pdf",sep="")
ggsave (pdf1, device="pdf", width=7, height=7, dpi=600)

##########################################################################################################################
