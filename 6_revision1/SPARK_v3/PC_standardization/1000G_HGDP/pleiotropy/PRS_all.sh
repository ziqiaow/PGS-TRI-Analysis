#!/bin/bash
#SBATCH --job-name=pleiotropy
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --array=1-22
#SBATCH --cpus-per-task=5


####################
readarray -t a < /dcs04/nilanjan/data/zwang/family/pleiotropy/redo/pgs_trait.txt

for i in {1..10}
do
/dcs04/nilanjan/data/zwang/tools/20240105/plink2 \
--threads 5 \
--mach-r2-filter \
--rm-dup exclude-all \
--bfile /dcs04/nilanjan/data/ydun/HGDP+1000G/GRCh38/original_files/chr${SLURM_ARRAY_TASK_ID} \
--remove /dcs04/nilanjan/data/ydun/HGDP+1000G/GRCh38/related_low_quality.id \
--score /dcs04/nilanjan/data/zwang/family/pleiotropy/redo/${a[$(($i-1))]}_score.txt 1 2 3 list-variants cols=+scoresums,-scoreavgs \
--out /dcs04/nilanjan/data/zwang/family/1000G_HGDP/pleiotropy/${a[$(($i-1))]}/score_chr$SLURM_ARRAY_TASK_ID  
done


