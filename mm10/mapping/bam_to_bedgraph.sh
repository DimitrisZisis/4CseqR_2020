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

DATA_PATH="/home/users/dzisis/4CseqR_new2018/mouse/mapping"
cd $DATA_PATH

# bam to bedgraph

bedtools genomecov -ibam ESC_pcdhb19_1_fragments_sorted.bam -bg > ESC_pcdhb19_1_fragments_sorted.bedGraph
bedtools genomecov -ibam iPSC_Pcdhb19_1_fragments_sorted.bam -bg > iPSC_Pcdhb19_1_fragments_sorted.bedGraph

exit


