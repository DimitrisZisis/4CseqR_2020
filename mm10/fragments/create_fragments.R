library(Biostrings)
library(Basic4Cseq)
library(BSgenome)
#file <- system.file("extdata", "TAIR10_chr_all.fas.gz", package="BSgenome")
fasta.seqlengths(file)

#create - ot call a BSgenome library for fragment preparation.
seed_files <- system.file("extdata", "GentlemanLab", package="BSgenome")
tail(list.files(seed_files, pattern="-seed$"))
Athaliana_seed <- list.files(seed_files, pattern="bsgenome.mmusculus.ucsc.mm10-seed$", full.names=TRUE)
cat(readLines(Athaliana_seed), sep="\n")
forgeBSgenomeDataPkg("/home/dimitris/R/x86_64-pc-linux-gnu-library/3.4/BSgenome/extdata/GentlemanLab/bsgenome.mmusculus.ucsc.mm10-seed")

library(bsgenome.mmusculus.ucsc.mm10)
genome <- bsgenome.mmusculus.ucsc.mm10

#Create the fragments for all chromosomes for mm10
libraryName2 = "fragments_mm10_2018_l82.csv"

require(BSgenome.Mmusculus.UCSC.mm10)
createVirtualFragmentLibrary(chosenGenome = Mmusculus, firstCutter = "AAGCTT", secondCutter = "GTAC", readLength = 82, onlyNonBlind = FALSE, libraryName = libraryName2)

