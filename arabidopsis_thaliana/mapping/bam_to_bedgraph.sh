#!/bin/bash

##PBS -l nodes=1:ppn=1:haswell
##PBS -l mem=10GB
##PBS -l walltime=6:00:00
##PBS -N sambam

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=5GB
#SBATCH -p standard

# convert bam to bedgraph 

#samtools=/home/users/igrbiom2/tools/samtools-1.1/samtools

module load bedtools

DATA_PATH="/home/users/dzisis/4CseqR_new2018/arabidopsis/mapping"
cd $DATA_PATH

# bam to bedgraph

bedtools genomecov -ibam Ath_FLC4C_A_1_un_sorted_test.bam -bg > Ath_FLC4C_A_1_un_sorted.bedGraph
bedtools genomecov -ibam Ath_FLC4C_B_1_un_sorted_test.bam -bg > Ath_FLC4C_B_1_un_sorted.bedGraph
#bedtools genomecov -ibam arabidopsis_3_1_fragments_un_sorted.bam -bg > arabidopsis_3_1_fragments_un_sorted.bedGraph
#bedtools genomecov -ibam arabidopsis_1_4_fragments_un_sorted.bam -bg > arabidopsis_1_4_fragments_un_sorted.bedGraph
#bedtools genomecov -ibam arabidopsis_2_4_fragments_un_sorted.bam -bg > arabidopsis_2_4_fragments_un_sorted.bedGraph
#bedtools genomecov -ibam arabidopsis_3_4_fragments_un_sorted.bam -bg > arabidopsis_3_4_fragments_un_sorted.bedGraph


exit


