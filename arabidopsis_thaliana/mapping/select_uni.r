gc()
rm(list=ls()) 

# library(ShortRead)

#exptPath = c("arabidopsis_1_1_fragments","arabidopsis_2_1_fragments","arabidopsis_3_1_fragments","arabidopsis_1_4_fragments",
#	"arabidopsis_2_4_fragments","arabidopsis_3_4_fragments")

exptPath = c("arabidopsis_1_1_2mis_2nd","arabidopsis_2_1_2mis_2nd","arabidopsis_3_1_2mis_2nd","arabidopsis_1_4_2mis_2nd","arabidopsis_2_4_2mis_2nd","arabidopsis_3_4_2mis_2nd")

for (j in 1:6){

	file <- paste(exptPath[j] ,".sam",sep="")

	#exptPath

#	lines <- readLines(file , 1000)
	lines <- readLines(file)

	pos1 <- grepl("XS:i",lines)
	#pos1[1:100]

	lines = subset(lines, pos1==FALSE)
	#lines

	file <- paste(exptPath[j] ,"_un.sam",sep="")

	writeLines(lines, file)

}


###################### koniec

