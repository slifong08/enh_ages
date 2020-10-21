#!/bin/bash
#SBATCH --mail-user=sarah.l.fong@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1    # comments allowed
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --mem=2GB
#SBATCH --output=/dors/capra_lab/users/fongsl/enh_age/

#SBATCH --job-name=get_oldest_age_shuf-roadmap
#SBATCH --account=capra_lab
#SBATCH --error=%A_%a.stdout
#SBATCH --output=%A_%a.stderr
#SBATCH --array=[0-229]%10

ml load GCC
ml load Anaconda3

# run all the shuffle files


FILES=(${1})
echo $SLURM_ARRAY_TASK_ID
echo ${FILES[$SLURM_ARRAY_TASK_ID]}
python /dors/capra_lab/users/fongsl/enh_age/bin/run_parallel_breaks.py ${FILES[$SLURM_ARRAY_TASK_ID]}
