#!/bin/bash

#SBATCH --partition=unlimited
#SBATCH --mem=16g

#module load python/3.6
module load R/4.0

Rscript GDC_vis.R
#python GDC_preproc.py
