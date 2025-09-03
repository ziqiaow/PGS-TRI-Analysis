#!/bin/bash
#SBATCH --job-name=family_child
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --partition=chatterjee

simulate.py 1000 0.5 /dcs04/nilanjan/data/zwang/family/simulation/snipar/simulated_data/ --nfam 100000 --n_am 20 --r_par 0.5 --save_par_gts


/dcs04/nilanjan/data/zwang/tools/20240105/plink2 \
--bfile /dcs04/nilanjan/data/zwang/family/simulation/snipar/simulated_data/chr_1 \
--score /dcs04/nilanjan/data/zwang/family/simulation/snipar/simulated_data/causal_effects.txt 1 2 4 cols=+scoresums,-scoreavgs \
--out /dcs04/nilanjan/data/zwang/family/simulation/snipar/simulated_data/PGS_children


/dcs04/nilanjan/data/zwang/tools/20240105/plink2 \
--bfile /dcs04/nilanjan/data/zwang/family/simulation/snipar/simulated_data/chr_1_par \
--score /dcs04/nilanjan/data/zwang/family/simulation/snipar/simulated_data/causal_effects.txt 1 2 4 cols=+scoresums,-scoreavgs \
--out /dcs04/nilanjan/data/zwang/family/simulation/snipar/simulated_data/PGS_parents