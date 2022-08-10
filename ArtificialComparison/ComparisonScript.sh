#!/bin/bash
#
# Script for the Comparison of both BwS approximations
#
#SBATCH --job-name=BwSAccuracyComparison
#SBATCH --partition=short
#SBATCH --time=2:00:00
#SBATCH --array=1-8
#
########################################################################
 
# Decide which input file this Slurm task will process.
# We use the special $SLURM_ARRAY_TASK_ID variable, which tells the
# task which one it is among the whole job array.
input_file=Parameters${SLURM_ARRAY_TASK_ID}.txt
output_file=output${SLURM_ARRAY_TASK_ID}.txt
 
./LRAndresArtificialComparison1 $ "Parameters${SLURM_ARRAY_TASK_ID}.txt"
