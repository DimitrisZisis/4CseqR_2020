
newdata.dt <- fread("Fisher_bintest_200000_50_alltable_mm10_new.csv")
newdata.dt[1:10,]
#remove NA rows
newdata.dt=newdata.dt[complete.cases(newdata.dt), ]
#remove rows with average probabilities <0,1
newdata.dt$mean <- rowMeans(newdata.dt[,c('proba in group con1', 'proba in group con1')], na.rm=TRUE)
newdata.dt$mean
clean.dt <- subset(newdata.dt, mean > 0.1)
trash.dt <- subset(newdata.dt, mean < 0.1)
#selecet if probabilies in condition 1 or 2 are > 0,1 
clean.dt2 <- newdata.dt[ which( newdata.dt$`proba in group con1` > 0.1 | newdata.dt$`proba in group con2` > 0.1) , ]

nrow(newdata.dt)
nrow(clean.dt)
nrow(trash.dt)

#calculate FDR for different cases
clean.dt[ , "Fisher_pvalue" := as.numeric(Fisher_pvalue)][ , "FDR_new" := p.adjust(Fisher_pvalue, method = "fdr")]
trash.dt[ , "Fisher_pvalue" := as.numeric(Fisher_pvalue)][ , "FDR_new" := p.adjust(Fisher_pvalue, method = "fdr")]
#fdr.dt[ , "Fisher_pvalue" := as.numeric(Fisher_pvalue)][ , "FDR_new" := p.adjust(Fisher_pvalue, method = "fdr")]
#check if there are sig contacts 
sig_cont=trash.dt[trash.dt$FDR_new<= 0.05,]
sig_cont[1:10,]


write.table(newdata.dt, file = "Fisher_bintest_200000_100_all_fdr_test.csv", col.names = T, row.names = F, quote = F, sep = "\t")


