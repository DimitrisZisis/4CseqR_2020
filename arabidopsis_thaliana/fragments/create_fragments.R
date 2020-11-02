library(Biostrings)
library(Basic4Cseq)
library(BSgenome)
#file <- system.file("extdata", "TAIR10_chr_all.fas.gz", package="BSgenome")
fasta.seqlengths(file)

#create - ot call a BSgenome library for fragment preparation.
seed_files <- system.file("extdata", "GentlemanLab", package="BSgenome")
tail(list.files(seed_files, pattern="-seed$"))
Athaliana_seed <- list.files(seed_files, pattern="BSgenome.Athaliana.TAIR.TAIR10-seed$", full.names=TRUE)
cat(readLines(Athaliana_seed), sep="\n")
forgeBSgenomeDataPkg("/home/dimitris/R/x86_64-pc-linux-gnu-library/3.4/BSgenome/extdata/GentlemanLab/BSgenome.Athaliana.TAIR.TAIR10-seed")

library(BSgenome.Athaliana.TAIR.TAIR9)
genome <- BSgenome.Athaliana.TAIR.TAIR9
seqlengths(genome)
genome$Chr1
libraryName = "fragments_Athalian_2018_read30.csv"

#Create the fragments for all chromosomes for Athaliana TAIR9
require(BSgenome.Athaliana.TAIR.TAIR9)
createVirtualFragmentLibrary(chosenGenome = Athaliana, firstCutter = "AGATCT", secondCutter = "CATG", readLength = 30, onlyNonBlind = FALSE, chromosomeName = Chr, libraryName = libraryName)

#Create the fragments for all chromosomes for mm10
libraryName2 = "fragments_mm10_2018_l82.csv"

require(BSgenome.Mmusculus.UCSC.mm10)
createVirtualFragmentLibrary(chosenGenome = Mmusculus, firstCutter = "AAGCTT", secondCutter = "GTAC", readLength = 82, onlyNonBlind = FALSE, libraryName = libraryName2)

