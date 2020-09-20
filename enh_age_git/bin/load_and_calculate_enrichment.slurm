#!/bin/bash
#SBATCH --mail-user=sarah.l.fong@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1    # comments allowed 
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
#SBATCH --mem=4GB
#SBATCH --output=/dors/capra_lab/users/fongsl/enh_age/

#SBATCH --job-name=load_and_calculate_enrichment
#SBATCH --account=capra_lab
#SBATCH --error=%J.stdout
#SBATCH --output=%J.stderr

ml load Anaconda3

python /dors/capra_lab/users/fongsl/enh_age/bin/load_and_calculate_enrichment.py ${1} ${2} ${3}