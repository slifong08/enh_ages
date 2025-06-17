#!/bin/bash
#SBATCH --mail-user=sarah.l.fong@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1    # comments allowed
#SBATCH --ntasks=1
#SBATCH --time=15:00:00
#SBATCH --mem=64GB
#SBATCH --output=/dors/capra_lab/users/fongsl/enh_age/

#SBATCH --job-name=age_enh
#SBATCH --account=capra_lab
#SBATCH --error=%J.stdout
#SBATCH --output=%J.stderr

ml load Anaconda3

python /dors/capra_lab/users/fongsl/enh_ages/core_v_derived/seqid/2_sequence_identity-parse_msa-slurm.py ${1}