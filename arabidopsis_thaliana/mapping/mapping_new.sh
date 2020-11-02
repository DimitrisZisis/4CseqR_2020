#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --cpus-per-task=1
#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=20GB
##SBATCH -p standard

DATA_PATH="/home/users/dzisis/4CseqR_new2018/arabidopsis/mapping"
cd $DATA_PATH

#Alignment with Bowtie 2
#chmod +x bowtie2-align

module load bowtie2


#max 2 mismatches, min score = .40*31 = 12.4

bowtie2 -q -x AGATCT_fragments_new --no-unal --score-min L,0,-0.40 -U arabidopsis_1_1_Perfect31.fq -S arabidopsis_1_1_test_new.sam

#bowtie2 -q -x AGATCT_fragments --no-unal --score-min L,0,-0.40 -U arabidopsis_2_1_Perfect28.fq -S arabidopsis_2_1_2mis_2nd.sam
#bowtie2 -q -x AGATCT_fragments --no-unal --score-min L,0,-0.40 -U arabidopsis_3_1_Perfect28.fq -S arabidopsis_3_1_2mis_2nd.sam
#bowtie2 -q -x AGATCT_fragments --no-unal --score-min L,0,-0.40 -U arabidopsis_1_4_Perfect28.fq -S arabidopsis_1_4_2mis_2nd.sam
#bowtie2 -q -x AGATCT_fragments --no-unal --score-min L,0,-0.40 -U arabidopsis_2_4_Perfect31.fq -S arabidopsis_2_4_2mis_2nd.sam
#bowtie2 -q -x AGATCT_fragments --no-unal --score-min L,0,-0.40 -U arabidopsis_3_4_Perfect29.fq -S arabidopsis_3_4_2mis_2nd.sam


#./bowtie2-align -q -x fragment_endsAllChr_AGATCT100length -a --score-min L,0,-0.4 -U arabidopsis_2_1_Perfect28.fq -S arabidopsis_2_1_AllChrnew_2mismore100length.sam
#./bowtie2-align -q -x fragment_endsAllChr_AGATCT100length -a --score-min L,0,-0.4 -U arabidopsis_2_2_Perfect28.fq -S arabidopsis_2_2_AllChrnew_2mismore100length.sam

