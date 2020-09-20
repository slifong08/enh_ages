#!/bin/bash
#SBATCH --mail-user=sarah.l.fong@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1    # comments allowed 
#SBATCH --ntasks=1
#SBATCH --time=01:30:00
#SBATCH --mem=2GB
#SBATCH --output=/dors/capra_lab/users/fongsl/enh_age/

#SBATCH --job-name=get_oldest_age_shuf-roadmap
#SBATCH --account=capra_lab
#SBATCH --error=%A_%a.stdout
#SBATCH --output=%A_%a.stderr
#SBATCH --array=[0-100]%10

ml load Anaconda3

# run all the shuffle files

FILES=(/dors/capra_lab/projects/enhancer_ages/fantom/data/shuffle/shuf-all_fantom_enh_112_tissues-*.bed)

echo $SLURM_ARRAY_TASK_ID
echo ${FILES[$SLURM_ARRAY_TASK_ID]}
python /dors/capra_lab/users/fongsl/enh_age/bin/get_oldest_age.py ${FILES[$SLURM_ARRAY_TASK_ID]}