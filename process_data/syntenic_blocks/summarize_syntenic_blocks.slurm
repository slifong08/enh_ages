#!/bin/bash
#SBATCH --mail-user=sarah.l.fong@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1    # comments allowed
#SBATCH --ntasks=1
#SBATCH --time=06:00:00
#SBATCH --mem=60GB
#SBATCH --output=/data/hodges_lab/ATAC-STARR_V2/bin_human-evolution/SF_hu-rhe

#SBATCH --job-name=hu_rhe
#SBATCH --account=capra_lab
#SBATCH --error=%J.stdout
#SBATCH --output=%J.stderr


ml load Anaconda3
source activate /home/fongsl/.conda/envs/sfenv/


python ${1}