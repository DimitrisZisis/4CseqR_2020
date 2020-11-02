library(ggplot2)
library(plyr)
require(graphics)
require(grDevices)
library(devtools)
require(reshape2) 

LMM_data <- read.csv("LMM_emmeans_200000_50_varf_qnorm_mm10.csv",header = T, sep="\t")
LMM_data[1:10,]
Fisher_data <- read.csv("Fisher_bintest_200000_50_all_fdr_mm10.csv",header = T, sep="\t")
Fisher_data[1:10,]

LMM_data$FDR
Fisher_data$Fisher_pvalue
bucket<-list(LMM=LMM_data$FDR,Fisher=Fisher_data$FDR)

ggplot(melt(bucket), aes(value, fill = L1)) + 
  #call geom_histogram with position="dodge" to offset the bars and manual binwidth of 2
  geom_histogram(position = "dodge")

#create scater plots for parameters in LMM and fiser 
Variant_A <- LMM_data$Est_meanCond1
Variant_B <- LMM_data$Est_meanCond2
# Create a new columns to see DCRs
LMM_data$DCW <- "NO"
# if LMM_FDR < 0.05  set as "YES" 
LMM_data$DCW[LMM_data$FDR < 0.05] <- "YES"
p <- ggplot(data=LMM_data, aes(x=Variant_A, y=Variant_B, col=DCW)) + geom_point()
print(p + ggtitle("Mean fragment coverage"))

#For fisher similar 
Variant_A <- Fisher_data$proba.in.group.con1
Variant_B <- Fisher_data$proba.in.group.con2

#plot(Variant_A, Variant_B, col=ifelse(Variant_A > 0.3, 2, 1)) +
#title(main="Probability of fragment coverage")
Fisher_data$DCW <- "NO"
# if LMM_FDR < 0.05  set as "YES" 
Fisher_data$DCW[Fisher_data$Fisher_pvalue < 0.05] <- "YES"
p <- ggplot(data=Fisher_data, aes(x=Variant_A, y=Variant_B, col=DCW)) + geom_point()
print(p + ggtitle("Probability of fragment coverage"))

#########################################################################################################


#Create volcano plot for LMM and Fisher
Variant_effect_LMM <- LMM_data$Est_meanCond2 - LMM_data$Est_meanCond1
# add a column of NAs
LMM_data$DCW <- "NO"
# if variant effect plays a role
#LMM_data$DCW[Variant_effect > 0.2 & LMM_data$FDR < 0.05] <- "NO"
# if FDR pvalue < 0.05, set as "YES"
LMM_data$DCW[ LMM_data$FDR < 0.05] <- "YES"
p1 <- ggplot(data=LMM_data, aes(x=Variant_effect_LMM, y=-log10(LMM_data$FDR), col=DCW)) + geom_point() + theme_minimal()
print(p1 + ggtitle("LMM analysis"))

#Create volcano plot for and Fisher
Variant_effect_fisher <- Fisher_data$proba.in.group.con2 - Fisher_data$proba.in.group.con1
# add a column of NAs
Fisher_data$DCW <- "NO"
# if Fisher_pvalue < 0.05, set as "YES"
Fisher_data$DCW[ Fisher_data$Fisher_pvalue < 0.05] <- "YES"
p2 <- ggplot(data=Fisher_data, aes(x=Variant_effect_fisher, y=-log10(Fisher_data$Fisher_pvalue), col=DCW)) + geom_point() + theme_minimal()
print(p2 + ggtitle("Fisher analysis"))

#Relationship between variant effects and significant scores 
new <- cbind(LMM_data$chr, LMM_data$FDR, Fisher_data$Fisher_pvalue, Fisher_data$FDR, Variant_effect_LMM,Variant_effect_fisher )
newdata=as.data.frame(new)
newdata$DCWs <- "NO_DCWs"
newdata$DCWs[ newdata$V2 < 0.05] <- "LMM_DCWs"
newdata$DCWs[ newdata$V3 < 0.05] <- "Fisher_DCWs"
newdata$DCWs[ newdata$V2 < 0.05 & newdata$V3 < 0.05] <- "BOTH_DCWs"

p3 <- ggplot(data=newdata, aes(x=Variant_effect_LMM, y=Variant_effect_fisher, col=DCWs)) + geom_point() + theme_minimal()
print(p3 + ggtitle("Comparison of results by two methods"))

#significant scores -log10(cor pvalue) for both methods 
newdata$V2
newdata$V4
LMM_test <- -log10(newdata$V2)
Fisher_test <- -log10(newdata$V3)
p4 <- ggplot(data=newdata, aes(x=LMM_test, y=Fisher_test, col=DCWs)) + geom_point() + theme_minimal()
print(p4 + ggtitle("Comparison of results by two methods/significant scores"))


