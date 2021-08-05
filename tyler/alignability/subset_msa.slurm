#!/bin/bash
#SBATCH --mail-user=sarah.l.fong@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1    # comments allowed
#SBATCH --ntasks=1
#SBATCH --time=01:30:00
#SBATCH --mem=8GB
#SBATCH --output=/dors/capra_lab/users/fongsl/tyler/data/alignability/

#SBATCH --job-name=CON_ACC
#SBATCH --account=capra_lab
#SBATCH --error=%J.stdout
#SBATCH --output=%J.stderr


ml load Anaconda3
ml load  GCC/8.2.0  OpenMPI/3.1.4 Biopython/1.74
source /dors/capra_lab/users/fongsl/sfenv/bin/activate

python /dors/capra_lab/users/fongsl/enh_ages/tyler/alignability/subset_msa.py 
