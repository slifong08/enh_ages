#!/bin/bash
#SBATCH --mail-user=sarah.l.fong@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1    # comments allowed 
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mem=30GB
#SBATCH --output=/dors/capra_lab/users/fongsl/enh_age/

#SBATCH --job-name=get_oldest_age_shuf-roadmap
#SBATCH --account=capra_lab
#SBATCH --error=%A_%a.stdout
#SBATCH --output=%A_%a.stderr

ml load Anaconda3

# run all the shuffle files

python /dors/capra_lab/users/fongsl/enh_age/bin/get_oldest_age.py ${1} -a 0