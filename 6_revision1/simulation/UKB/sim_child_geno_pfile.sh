#!/bin/bash
#SBATCH --job-name=child_gen
#SBATCH --mem=200G
#SBATCH --nodes=1
#SBATCH --array=1-1000
#SBATCH --cpus-per-task=1


####################

module load R

R CMD BATCH --no-save --no-restore sim_child_geno_pfile.R script_$SLURM_ARRAY_TASK_ID