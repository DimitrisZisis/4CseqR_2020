library(parallel)
library(data.table)
library(GenomicRanges)
library(Biostrings)
library(lmerTest)
library(emmeans)
library(RVAideMemoire)
# make newdata as data.table for efficient operations
newdata.dt <- fread("totalnorm_all_mm10_binary_zero.txt")

# select sliding window length and step 
window <- 10000000; step <- 50000

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
FisherPerWin <- function(sliding_window, fragment_table, window_numb){
  
  # select fragments share the same chr as sliding window X AND that overlap in any amount with sliding window X
  newdatatable <- fragment_table[(V1 %in% sliding_window$seqnames)&((V2 <= sliding_window$start) & (V3 >= sliding_window$end)|(V2 > sliding_window$start) & (V2 < sliding_window$end)|(V3 > sliding_window$start) & (V3 < sliding_window$end)) , ]
  whichFragments <- newdatatable$Fragm_chr
  # perform operations to bring all values from the two conditions to format for stat testing
  flabel <- factor(c(rep("con1", nrow(newdatatable) * 3), rep("con2", nrow(newdatatable) * 3)), levels = c("con1", "con2"))
  fsample <- factor(c(rep("sample1", nrow(newdatatable)), rep("sample2", nrow(newdatatable)),rep("sample3", nrow(newdatatable)),rep("sample4", nrow(newdatatable)),rep("sample5", nrow(newdatatable)),rep("sample6", nrow(newdatatable))), levels = c("sample1", "sample2","sample3","sample4","sample5","sample6"))
  fragnum <- factor(nrow(newdatatable))
  data1 = unlist(c(newdatatable[ , 4], newdatatable[ , 5], newdatatable[ , 6], newdatatable[ , 7], newdatatable[ , 8], newdatatable[ , 9]))
  datan <- as.double(data1)
  #datab <- as.binary(data1)
  datav = data.table(flabel, fsample, whichFragments, datan, fragnum)

  #add column tp datav for window number
  datav$windownum <- window_numb
  #Create a table with covered uncovered fragments 
  
  datan_check <- length(datav[, unique(datan)])
  # if (length( data1_check[data1_check == T]) >= 1) {
  if (datan_check == 2) {
    fisher <- fisher.bintest(formula = datan ~ flabel, data = datav, alpha = 0.05, p.method = "fdr")
    pval <- fisher$p.value
    prob_con1 <- as.data.frame(fisher$estimate["proba in group con1"])[1,1]
    prob_con2 <- as.data.frame(fisher$estimate["proba in group con2"])[1,1]
    t1  <- rowSums(newdatatable[,c("V4", "V5", "V6")])>0
    t2  <- rowSums(newdatatable[,c("V7", "V8", "V9")])>0
    
    sum(t1) 
    sum(!t1) 
    sum(t2) 
    sum(!t2)
    #extract number of 1s and zeros out of the total of each condition 
    sum1s_con1 <- sum(newdatatable[,c("V4", "V5", "V6")])
    sum1s_con2 <- sum(newdatatable[,c("V7", "V8", "V9")])
    
    nr_con1=nrow(newdatatable[,c("V4", "V5", "V6")])
    nc_con1=ncol(newdatatable[,c("V4", "V5", "V6")])
    
    total_con1=nr_con1*nc_con1
    
    nr_con2=nrow(newdatatable[,c("V7", "V8", "V9")])
    nc_con2=ncol(newdatatable[,c("V7", "V8", "V9")])
    
    total_con2=nr_con2*nc_con2
    
  } else {
   
    pval <- NA
    prob_con1 <- NA
    prob_con2 <- NA
    sum1s_con1 <- NA
    sum1s_con2 <- NA
    total_con1 <- NA
    total_con2 <- NA
    t1 <- NA
    t2 <- NA
 }
  
 
  return(c(as.vector(sliding_window$seqnames), sliding_window$start, sliding_window$end, pval, prob_con1, prob_con2, as.numeric(as.vector(fragnum)), sum1s_con1, sum1s_con2, total_con1, total_con2))
  }

# apply function for each row in slwin.dt (for each sliding window)
# Extract specific windows nrow(slwin.dt) 
pvaluePerWin <- lapply(1:nrow(slwin.dt), function(idx) FisherPerWin(sliding_window = slwin.dt[idx, ], fragment_table = newdata.dt, window_numb = idx))
#cat(capture.output(print(pvaluePerWin), file="sliding_wind_1_20_arab.csv"))
#sliding_window <- slwin.dt[1,]; fragment_table <- newdata.dt; window_numb <- 10;
#rm(sliding_window, fragment_table, window_numb)
# for parallel execution:
#pvaluePerWin <- mclapply(1:nrow(slwin.dt), function(idx) kruskalPerWin(sliding_window = slwin.dt[idx, ], fragment_table = newdata.dt), mc.cores = 4 )


# make matrix with ncol = (number of things returned by each iteration), fill it by row
pval.mt <- matrix(unlist(pvaluePerWin), ncol = 11, byrow = T)
# make matrix to data table and set correct column names
output.dt <- setnames(as.data.table(pval.mt), c("chr", "start", "end","Fisher_pvalue", "proba in group con1", "proba in group con2", "number_of_fragments", "sums_covered_con1", 'sums_covered_con2', "total_elemets_con1","total_elemets_con2" ))

# order data.table by (CHARACTER) chr and start
setorderv(output.dt, c("chr", "start"))
# write file to getwd()
write.table(output.dt, file = "Fisher_bintest_10000000_50_mm10.csv", col.names = T, row.names = F, quote = F, sep = "\t")



##############################################################################################################################



