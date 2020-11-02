#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=5
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=10GB
##SBATCH -p standard


BT=/home/users/dzisis/bedtools-2.17.0/bin/fastaFromBed

DATA_PATH="/home/users/dzisis/PaperMethods/GSE50029_MouseData"
cd $DATA_PATH

$BT -fi Mus_musculus.GRCm38.dna.toplevel.fa -bed valid_fragments_mm10.bed -s -fo valid_fragments_mm10.fa

exit






