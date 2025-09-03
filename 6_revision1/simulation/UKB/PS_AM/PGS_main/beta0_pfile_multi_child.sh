#!/bin/bash
#SBATCH --job-name=family_child
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1



module load R

R CMD BATCH --no-save --no-restore beta0_negbmi_pfile_multi_child.R 