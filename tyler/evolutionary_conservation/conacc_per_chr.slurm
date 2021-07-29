#!/bin/bash
#SBATCH --mail-user=sarah.l.fong@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1    # comments allowed
#SBATCH --ntasks=1
#SBATCH --time=03:00:00
#SBATCH --mem=8GB
#SBATCH --output=/dors/capra_lab/users/fongsl/tyler/data/

#SBATCH --job-name=CON_ACC
#SBATCH --account=capra_lab
#SBATCH --error=%J.stdout
#SBATCH --output=%J.stderr


ml load Anaconda3


python /dors/capra_lab/users/fongsl/enh_ages/tyler/evolutionary_conservation/conacc_per_chr.py ${1} -br ${2} -msa ${3} -mod ${4}
