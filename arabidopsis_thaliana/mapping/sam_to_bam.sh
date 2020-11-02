#!/bin/bash

##PBS -l nodes=1:ppn=1:haswell
##PBS -l mem=10GB
##PBS -l walltime=6:00:00
##PBS -N sambam

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --cpus-per-task=1
#SBATCH --time=00:45:00
#SBATCH --mem-per-cpu=5GB
#SBATCH -p fast

# przerob an wig wszystkich zmapowanych, nie tylko unique

#samtools=/home/users/igrbiom2/tools/samtools-1.1/samtools

module load samtools

#WigSoft=/home/users/dzisis/mapping/peakranger

DATA_PATH="/home/users/dzisis/4CseqR_new2018/arabidopsis/mapping"
cd $DATA_PATH

###### sam to bam

samtools view -bS arabidopsis_1_1_2mis_2nd_un.sam > Ath_FLC4C_A_1_un.bam
samtools view -bS arabidopsis_1_4_2mis_2nd_un.sam > Ath_FLC4C_B_1_un.bam
#samtools view -bS arabidopsis_2_1_2mis_2nd.sam > Ath_FLC4C_A_2_n.bam
#samtools view -bS arabidopsis_3_1_2mis_2nd.sam > Ath_FLC4C_A_3_n.bam
#samtools view -bS arabidopsis_1_4_2mis_2nd.sam > Ath_FLC4C_B_1_n.bam
#samtools view -bS arabidopsis_2_4_2mis_2nd.sam > Ath_FLC4C_B_2_n.bam
#samtools view -bS arabidopsis_3_4_2mis_2nd.sam > Ath_FLC4C_B_3_n.bam

# sort, index
samtools sort Ath_FLC4C_A_1_un.bam Ath_FLC4C_A_1_un_sorted_test
samtools sort Ath_FLC4C_B_1_un.bam Ath_FLC4C_B_1_un_sorted_test
#samtools sort Ath_FLC4C_A_2_n.bam Ath_FLC4C_A_2_sorted
#samtools sort Ath_FLC4C_A_3_n.bam Ath_FLC4C_A_3_sorted
#samtools sort Ath_FLC4C_B_1_n.bam Ath_FLC4C_B_1_sorted
#samtools sort Ath_FLC4C_B_2_n.bam Ath_FLC4C_B_2_sorted
#samtools sort Ath_FLC4C_B_3_n.bam Ath_FLC4C_B_3_sorted
 
samtools index Ath_FLC4C_A_1_un_sorted_test.bam
samtools index Ath_FLC4C_B_1_un_sorted_test.bam
#samtools index Ath_FLC4C_A_2_sorted.bam
#samtools index Ath_FLC4C_A_3_sorted.bam
#samtools index Ath_FLC4C_B_1_sorted.bam
#samtools index Ath_FLC4C_B_2_sorted.bam
#samtools index Ath_FLC4C_B_3_sorted.bam

exit


##Important terminal commands 

#samtools view -h Ath_FLC4C_A_1_n.bam | sed "s/()//" | samtools reheader - Ath_FLC4C_A_1_n.bam > Ath_FLC4C_A_1_new.bam
