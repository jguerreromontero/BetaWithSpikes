#!/bin/bash
#
# Script for the Comparison of both BwS approximations
#
#SBATCH --job-name=BwSAccuracyComparison
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --array=1-3
#
########################################################################
 
# Decide which input file this Slurm task will process.
# We use the special $SLURM_ARRAY_TASK_ID variable, which tells the
# task which one it is among the whole job array.
input_file=Parameters${SLURM_ARRAY_TASK_ID}
output_file=output${SLURM_ARRAY_TASK_ID}
filename="Parameters${SLURM_ARRAY_TASK_ID}"

./ArtificialComparisonMoments.exe $filename
