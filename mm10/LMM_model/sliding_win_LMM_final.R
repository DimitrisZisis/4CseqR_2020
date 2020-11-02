library(parallel)
library(data.table)
library(GenomicRanges)
library(Biostrings)
library(lmerTest)
library(emmeans)
# make newdata as data.table for efficient operations
newdata.dt <- fread("totalnorm_all_mm10_qnorm_zero.txt")
#newdata.dt <- as.data.table(newdata)

# select sliding window length and step 
window <- 100000; step <- 20000

# hold the max possible value chromosomes have in the fragment table
chrEnds <- setnames(newdata.dt[, max(V3), by = "V1"], c("chr", "maxEnd"))
# make a data.table and then a GRanges object from those chr ends
chrLen <- data.table(chr = chrEnds$chr, start = 0, end = chrEnds$maxEnd, strand = ".", score = 0)
chrLen.gr <- makeGRangesFromDataFrame(chrLen)
# use GRanges object to make a sliding window. Use unlist() to ensure that you get back a GRanges and NOT a GRangeslist object.
slWin.gr <- unlist(slidingWindows(chrLen.gr, width = window, step = step))
# turn sliding window to data.table and we have all we need
slwin.dt <- as.data.table(slWin.gr)

#add column by chr with the ID of each fragment 
chrNums <- newdata.dt[, .N , by = V1]
newdata.dt$Fragm_chr <- unlist(lapply(1:nrow(chrNums), function(x) seq.int(1:chrNums[x, N])))
newdata.dt[, "Fragm_chr" := paste0(V1, "_", Fragm_chr)]

#newdata.dt[, Fragm_chr := nrow(newdata.dt), by = "V1"]

# function to find which fragments overlap with each sl window and run stats on them
LMMPerWin <- function(sliding_window, fragment_table, window_numb){
  
  # select fragments share the same chr as sliding window X AND that overlap in any amount with sliding window X
  newdatatable <- fragment_table[(V1 %in% sliding_window$seqnames)&((V2 <= sliding_window$start) & (V3 >= sliding_window$end)|(V2 > sliding_window$start) & (V2 < sliding_window$end)|(V3 > sliding_window$start) & (V3 < sliding_window$end)) , ]
  whichFragments <- newdatatable$Fragm_chr
  # perform operations to bring all values from the two conditions to format for stat testing
  flabel <- factor(c(rep("con1", nrow(newdatatable) * 3), rep("con2", nrow(newdatatable) * 3)), levels = c("con1", "con2"))
  fsample <- factor(c(rep("sample1", nrow(newdatatable)), rep("sample2", nrow(newdatatable)),rep("sample3", nrow(newdatatable)),rep("sample4", nrow(newdatatable)),rep("sample5", nrow(newdatatable)),rep("sample6", nrow(newdatatable))), levels = c("sample1", "sample2","sample3","sample4","sample5","sample6"))
  fragnum <- factor(nrow(newdatatable))
  data1 = unlist(c(newdatatable[ , 4], newdatatable[ , 5], newdatatable[ , 6], newdatatable[ , 7], newdatatable[ , 8], newdatatable[ , 9]))
  datav = data.table(flabel, fsample, whichFragments, data1, fragnum)
  #add column tp datav for window number
  datav$windownum <- window_numb
  
  #run LMM test in case data1 contains non-zero vals
  data1_check <- data1 != 0
  if (length( data1_check[data1_check == T]) >= 1) {
    lmodel<- lmer(formula = data1 ~ flabel + (1 | whichFragments), data = datav, control=lmerControl(check.nlev.gtr.1 = "ignore"))
    sfit <-summary(lmodel)
    serror <- sfit$coefficients[2,2] 
    pval <- sfit$coefficients[2,5]
    #Extract the variance of random effect component -fragments 
    sfit2 <- as.data.frame(VarCorr(lmodel))
    frag_variance <- sfit2$vcov[1]
    means <-  emmeans(lmodel, "flabel", mode = "satterth")
    test<-summary(means)
    con1mean <- test$emmean[1]
    con2mean <- test$emmean[2]
  } else {
    serror <- NA
    pval <- NA
    frag_variance <- NA
    con1mean <- NA
    con2mean <- NA
  }
  
  # calculate Estimated marginal means (Least-squares means)
  #means <-  emmeans(lmodel, "flabel", mode = "satterth")
 # test<-summary(means)
  #con1mean <- mean(datav[1:nrow(datav)/2, data1])
  #con2mean <- mean(datav[(nrow(datav)/2 + 1):nrow(datav), data1])
  #return(datav)
  return(c(as.vector(sliding_window$seqnames), sliding_window$start, sliding_window$end, fragnum, con1mean, con2mean, serror, pval, frag_variance))
}

# apply function for each row in slwin.dt (for each sliding window)

# Extract specific windows 
pvaluePerWin <- lapply(1:nrow(slwin.dt), function(idx) LMMPerWin(sliding_window = slwin.dt[idx, ], fragment_table = newdata.dt, window_numb = idx))
#cat(capture.output(print(pvaluePerWin), file="sliding_wind_1_20_arab.csv"))
sliding_window <- slwin.dt[1, ]; fragment_table <- newdata.dt; window_numb <- 1;
 #rm(sliding_window, fragment_table, window_numb)
# for parallel execution:
 #pvaluePerWin <- mclapply(1:nrow(slwin.dt), function(idx) kruskalPerWin(sliding_window = slwin.dt[idx, ], fragment_table = newdata.dt), mc.cores = 4 )


# make matrix with ncol = (number of things returned by each iteration), fill it by row
pval.mt <- matrix(unlist(pvaluePerWin), ncol = 9, byrow = T)
# make matrix to data table and set correct column names
output.dt <- setnames(as.data.table(pval.mt), c("chr", "start", "end","Fragments_No", "Est_meanCond1", "Est_meanCond2", "LMM_SE", "LMM_pvalue", "variance_fragm"))
# make kruskal wallis column numeric (it was character) and make a column with the pvalue adjustments
output.dt[ , "LMM_pvalue" := as.numeric(LMM_pvalue)][ , "FDR" := p.adjust(LMM_pvalue, method = "fdr")]

# order data.table by (CHARACTER) chr and start
setorderv(output.dt, c("chr", "start"))
# write file to getwd()
write.table(output.dt, file = "LMM_emmeans_100000_20_fram_results_qnorm_arab.csv", col.names = T, row.names = F, quote = F, sep = "\t")



##############################################################################################################################


