#!/bin/bash
#
# Script for the Comparison of both BwS approximations
#
#SBATCH --job-name=VerbAnalysis
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --array=1-4
#SBATCH --mem=10M
#
########################################################################
 
# Decide which input file this Slurm task will process.
# We use the special $SLURM_ARRAY_TASK_ID variable, which tells the
# task which one it is among the whole job array.
input_file=LangData${SLURM_ARRAY_TASK_ID}.txt
output_file=output${SLURM_ARRAY_TASK_ID}.txt
filename="LangData${SLURM_ARRAY_TASK_ID}.txt"

./VerbAnalysis.exe $filename
