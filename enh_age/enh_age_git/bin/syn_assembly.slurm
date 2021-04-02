#!/bin/bash
#SBATCH --mail-user=sarah.l.fong@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1    # comments allowed
#SBATCH --ntasks=1
#SBATCH --time=00:20:00
#SBATCH --mem=4GB
#SBATCH --output=/dors/capra_lab/users/fongsl/enh_age/

#SBATCH --job-name=age_enh
#SBATCH --account=capra_lab
#SBATCH --error=%J.stdout
#SBATCH --output=%J.stderr


ml load Anaconda3

# ${1} = /dors/capra_lab/users/fongsl/enh_age/data/<<FILE_TO_AGE>>.bed # str

echo $(pwd)/${1}
python /dors/capra_lab/users/fongsl/enh_age/enh_age_git/bin/syntenic_assembly.py  $(pwd)/${1}
