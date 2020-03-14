#!/bin/bash

#SBATCH --partition=unlimited
#SBATCH --mem=20g

module load python/3.6
module load R/3.4

python CONETT_input.py \
/data/BMR_Genomics/data_for_gardner/exomeseq_activeDev/mutect_out

