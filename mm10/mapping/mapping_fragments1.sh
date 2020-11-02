#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --cpus-per-task=1
#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=20GB
##SBATCH -p standard



DATA_PATH="/home/users/dzisis/4CseqR_new2018/mouse/fragments"
cd $DATA_PATH

#Alignment with Bowtie 2
#chmod +x bowtie2-align

module load bowtie2

#max 2 mismatches

bowtie2 -q -x mm10_AAGCTT_fragments --no-unal --score-min L,0,-0.25 -U ESC_pcdhb19_1_filtered.fastq -S ESC_pcdhb19_1_fragments.sam
bowtie2 -q -x mm10_AAGCTT_fragments --no-unal --score-min L,0,-0.25 -U ESC_pcdhb19_2_filtered.fastq -S ESC_pcdhb19_2_fragments.sam
bowtie2 -q -x mm10_AAGCTT_fragments --no-unal --score-min L,0,-0.25 -U ESC_pcdhb19_3_filtered.fastq -S ESC_pcdhb19_3_fragments.sam
bowtie2 -q -x mm10_AAGCTT_fragments --no-unal --score-min L,0,-0.25 -U iPSC_Pcdhb19_1_filtered.fastq -S iPSC_Pcdhb19_1_fragments.sam
bowtie2 -q -x mm10_AAGCTT_fragments --no-unal --score-min L,0,-0.25 -U iPSC_Pcdhb19_2_filtered.fastq -S iPSC_Pcdhb19_2_fragments.sam
bowtie2 -q -x mm10_AAGCTT_fragments --no-unal --score-min L,0,-0.25 -U iPSC_Pcdhb19_3_filtered.fastq -S iPSC_Pcdhb19_3_fragments.sam

exit




