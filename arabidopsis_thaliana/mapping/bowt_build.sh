#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=2GB
#PBS -l walltime=00:20:00
#PBS -N 4CSEQbuild

#BT=/home/users/dzisis/bedtools-2.17.0/bin/bedtools

DATA_PATH="/home/users/dzisis/4CseqR_new2018/arabidopsis/mapping"
cd $DATA_PATH

module load bowtie2

bowtie2-build --offrate 1 AGATCT_fragments_new.fa  AGATCT_fragments_new
