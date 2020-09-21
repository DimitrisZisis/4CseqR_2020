#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=4GB
#PBS -l walltime=01:00:00
#PBS -N btools

BT=/home/users/dzisis/bedtools-2.17.0/bin/fastaFromBed

DATA_PATH="/home/users/dzisis/PaperMethods/ArabidopsisData"
cd $DATA_PATH

$BT -fi Arabidopsis_thaliana.TAIR10.dna.chromosome.1.fa -bed AGATCT_fragments.bed -s -fo AGATCT_fragments.fa

exit






