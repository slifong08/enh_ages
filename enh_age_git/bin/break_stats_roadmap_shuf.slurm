#!/bin/bash
#SBATCH --mail-user=sarah.l.fong@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1    # comments allowed 
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
#SBATCH --mem=2GB
#SBATCH --output=/dors/capra_lab/users/fongsl/enh_age/

#SBATCH --job-name=age_enh
#SBATCH --account=capra_lab
#SBATCH --error=%A_%a.stdout
#SBATCH --output=%A_%a.stderr
#SBATCH --array=[0-819]%20

ml load Anaconda3
FILES=(/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/shuffle/breaks/shuf-Hsap_H3K27ac_plus_H3K4me3_minus_E*.bed)

echo $SLURM_ARRAY_TASK_ID
echo ${FILES[$SLURM_ARRAY_TASK_ID]}
python /dors/capra_lab/users/fongsl/enh_age/bin/break_stats.py ${FILES[$SLURM_ARRAY_TASK_ID]}
