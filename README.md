# 4CseqR_2020_version
4CseqR: A new pipeline to analyse 4C-seq data and investigate DNA contacts across the genome in different species 

4CSeqR is a new innovative method-pipeline based on common tools. The main goal of this pipeline-schema, which is based on existing tools and simple scripts in R statistical language, is to investigate DNA contacts across the genome based on 4C-seq data and to suggest some modifications of restriction fragment coverage estimation and normalization based on known algorithms and show advantages and possible weaknesses compared to other existing tools. The selected computational methods are studied and discussed on every step from the pre-processing of next-generation sequencing reads, through the mapping and treatment of mapped reads, the estimation and adjustment of read coverage, until the calling of significant contacts and the comparative statistical analysis of experimental variants. 4CSeqR main algorithm is developed anew in R. 

For reviews about the 4C-seq methodology refer to the following articles :

1.  Dimitrios Zisis , Paweł Krajewski, Maike Stam , Blaise Weber, Iris Hövel (2019) Analysis of 4C-seq data: A comparison of methods J Bioinform Comput Biol Feb;18(1):2050001. doi: 10.1142/S0219720020500018.

2.  Raviram R, Rocha PP, Bonneau R et al. (2014) Interpreting 4C-Seq data: how far can we go.  Epigenomics  6(5): 455-457. 

The required R packages to run 4CSeqR are:Bsgenome, ShortRead, Biostrings, GenomicRanges, ggplot2, stringr, dplyr, plyr, lmerTest, data.table, emmeans, RVAideMemoire, parallel

```
library(Biostrings)
library(BSgenome)
library(ShortRead)
showMethods(readFastq)
showMethods(writeFastq)
library(GenomicRanges)
library(ggplot2)
...
```
In order to be able to run all steps of 4CseqR the following tools are required : bowtie2, bedtools, salmon and an installation of Python.

The input to 4CseqR is the data from a 4C-seq experiment in fastq format. The 4CseqR pipeline provides script from pre-processing to mapping reads to the reference genome or to a library of genome, estimate coverage, and normalize it providing tables at each step. By using *run_4CseqR_pipeline* and *4CseqR_functions* a user can create normalized  reanked  counts of reads and collect them in a table with the script *collect_results.R*. 
Further to normalization, a new statistical approach is provided in R . The user can choose either LMM approach or Fisher binary test to analyze the counts. For each step of the analysis, 4CseqR provides tables and graphs with the use of simple R scripts, for further exploration of the data.  

**1. Library of fragments**

To create the library of fragments for the selected reference genome, in 4CseqR we are using the function *createVirtualFragmentLibrary* from Basic4Cseq tool. Two restriction enzymes have to be specified that cut the DNA, the read length is needed to check the fragment ends of corresponding length for uniqueness. Filtering options (minimum and maximum size) are provided on fragment level and on fragment end level. The following exapmle is used to create a library of fragments for Arabidopsis thaliana (TAIR10) : 
```
seed_files <- system.file("extdata", "GentlemanLab", package="BSgenome")
tail(list.files(seed_files, pattern="-seed$"))
Athaliana_seed <- list.files(seed_files, pattern="BSgenome.Athaliana.TAIR.TAIR9-seed$", full.names=TRUE)
cat(readLines(Athaliana_seed), sep="\n")

libraryName = "fragments_Athalian.csv"

require(BSgenome.Athaliana.TAIR.TAIR9)
createVirtualFragmentLibrary(chosenGenome = Athaliana, firstCutter = "AGATCT", secondCutter = "CATG", readLength = 50, onlyNonBlind = FALSE, chromosomeName = Chr, libraryName = libraryName)

```
In case of the example experiment only fragments longer than 100 nucleotides were used for mapping. For this reason there is  the *data.filtered* function to filter fragments:
```
fragmentLenght = 100
newfragments=data.filtered(libraryName, fragmentLenght)
filteredLibraryName = "fragments_Athaliana100.csv"
write.table(newfragments, file = filteredLibraryName,quote = FALSE,row.names=FALSE , sep="\t")
```

The resulting table with library of fragmens is saved in .bed format and is used in bedtools to create a fasta file which will be the input to the mapping procedure. A script (*getfasta.sh*) has been provided for this reason. 

**2. Data Pre processing**

4CseqR takes as input fastq files directly after sequencing. It filters the original reads to keep only those that start directly at a restriction enzyme cutting site. The NGS data obtained for all samples (fastq format) are at the beginning processed by filtering out truncated reads and the reads that did not contain the primary restriction site 
Main purpose of this preprocessing step is to find the reads from the original files which are legitimate 4C reads which means that they contain the sequence (primer + restr.site ). For filtering and trimming the primer sequences from the fastq reads the R function *select.reads* is provided. The user have to create two text files (primer_filename, restriction_filename) with the primer sequences and the enzyme recognition sequences for each experiment  and also the fastq read files (fastq_filename). 
The restiction_filename which includes the restriction enzyme is also used as pattern in the trimming step of this function for the specific experiment.
The following example is using the raw reads of Arabidopsis thaliana and the restriction enzyme AGATCT (BglII) as filtering parameter.
```
primer_filename = "/home/dimitris/4CseqR/Data/primer_seq.txt"
restriction_filename = "/home/dimitris/4CseqR/Data/restriction_seq.txt"
fastq_filename = "/home/dimitris/4CseqR/Data/INDEX_1_5.fq"
PrimerSeq <- readLines(primer_filename)
RestrEnzyme <- readLines(restriction_filename)
SelectReads=select.reads(fastq_filename, PrimerSeq[1], RestrEnzyme[1], pattern)
writeFastq(SelectReads, paste("INDEX_1_5", sep="_", "treatment_results.fq"))

```

Another additional process which can be done in preprocessing is the filtering of short reads in the output file which can be done by grep in the same r script. for example in arabidopsis thaliana data we want to filter the reads based on lenght 31. 
```
select_reads <- grep("^.{31}$", sread(SelectReads))
SelectReads <- SelectReads[select_reads]
writeFastq(SelectReads, paste("INDEX_1_5", "_treatment_results31.fq", sep=""))
```
**3. Mapping of reads**

The reads were mapped to the library of fragments created based on the genome of the selected species used as the reference.In order to create an index file the bowtie2-build script is used.  Mapping was conducted by Bowtie2 in all positions with up to 2 mismatches allowed in the contact region only and without mismatches in the restriction sequence. The reads mapped to regerence genome or library offragments  are combined for further analysis. A script (*mapping.sh*) is provided for this step. 
```
bowtie2-build --offrate 1 Library_fragments_AGATCT.fa  Library_fragments_AGATCT

bowtie2 -q -x AGATCT_fragments_new --no-unal --score-min L,0,-0.40 -U arabidopsis_1_1_Perfect31.fq -S arabidopsis_1_1_test_new.sam
```


**4. Estimate Coverage**

4CseqR provides a script (*run_salmon.sh*) which is taking care about the non-uniquely mapped reads  especially in the pericentromeric region of each chromosome.  The pipeline applies the mapper Salmon from R.Patro et, al. [Nat Methods 2017 ]. Salmon is a tool that is developed to perform quantification of gene expression in RNA-seq experiments. Salmon supports two modes either with a FASTA file containing a reference genome and a set of reads in a FASTQ file or a set of pre-computed alignments in a SAM/BAM file. We use the second mode of Salmon by taking the output of mapping with the alignments in SAM format and the library of fragments as reference in FASTA. 
A bash script to run salmon (salmon.sh) is provided by the 4CseqR pipeline. 
```
salmon quant -t Library_fragments_AGATCT.fa --libType A -a sample1.sam -o salmon_quant_sample1
```

Command line,bedtools and a python script (*merge_4cnew.py*) are used in order to create the proper input for normalization. An R script (*csv_to_bed.r*) is provided to create bed files in case of visalizations in a genome broswer and .csv files are used as input in the normalization function. 

**5. Normalization**

The counts of reads are normalized according to the procedure based on ranking within classes of fragments (blind/non blind, length of ends, total length). Ranks within each category are expressed in the scale from 0 to 1 and the total coverage is expressed as a value from 0 to 2. The input of the R function is the input_filename and the categories of distance and length (lefS_categories_limits, rigS_categories_limits,lenght_categories_limits) based to the needs of each experiment. To normalize the example data:
```
input_filename = "/home/dimitris/4CseqR/Data/finaltable1_5leng100.csv"
lefS_categories_limits = c(50, 100, 150, 200, 250, 300, 350, 400)
rigS_categories_limits = c(50, 100, 150, 200, 250, 300, 350, 400)
lenght_categories_limits = c(50, 100, 150, 200, 400, 600, 800)
NormalizationRanks= normalization_ranks(input_filename,lefS_categories_limits, rigS_categories_limits, lenght_categories_limits)```
```
**6. Linear mixed model (LMM)**
The normalized ranks with the coverage of restriction fragments are used as input to the linear mixed model (LMM) as a tool of data analysis in order to find genomic windows which can be inferred to be in contact with the bait to a different degree in different experimental variants. The 4CseqR performs a linear mixed model analysis through a function called *LMMPerWin*.
For a particular sliding window, with specific step in which there is a number of restriction fragments, an LMM model is created for the normalized coverage as the outcome variable and is fitted using the lmerTest package in R. 

An R script called *sliding_window_LMM.R* is provided for the LMM analysis.

**7. Binarization**

The binarization is performed on the estimated coverage obtained from Salmon. As it is described in Methods, the binary implementation is done with respect to a threshold (P). The proposed transformation takes the estimated coverage of fragments and converts it into 0 or 1 based on the threshold P = 1.

An R script called *binarization.R* is available for this step.

**8. Fisher exact test**

Following binarization of estimated fragment coverage, a Fisher exact test is used to identify DCWs and DCRs. As input to the Fisher exact test we use a table in .csv format with binary results for all replications under each treatment.The 4CseqR performs the fisher exacts test analysis through a function called *FisherPerWin*. 
For a particular sliding window, with specific step in which there is a number of restriction fragments, a fisher test is created for the binarized coverage as the outcome variable and is fitted using the fisher.test package in R. 

An R script called *sliding_window_fisher.R* is provided for the fisher analysis.



