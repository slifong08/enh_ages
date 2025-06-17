#!/bin/bash
#SBATCH --mail-user=sarah.l.fong@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1    # comments allowed
#SBATCH --ntasks=1
#SBATCH --time=06:00:00
#SBATCH --mem=8GB
#SBATCH --output=/dors/capra_lab/users/fongsl/enh_age/

#SBATCH --job-name=hepg2_ELS
#SBATCH --account=capra_lab
#SBATCH --error=%J.stdout
#SBATCH --output=%J.stderr


ml load Anaconda3


python /dors/capra_lab/users/fongsl/enh_ages/landscape/ENCODE3_cCRE/HepG2/HepG2_ELS.py 