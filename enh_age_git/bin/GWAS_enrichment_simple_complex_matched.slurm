#!/bin/bash
#SBATCH --mail-user=sarah.l.fong@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1    # comments allowed 
#SBATCH --ntasks=1
#SBATCH --time=03:30:00
#SBATCH --mem=2GB
#SBATCH --output=/dors/capra_lab/users/fongsl/enh_age/

#SBATCH --job-name=GWAS_enrichment
#SBATCH --account=capra_lab
#SBATCH --error=%A_%a.stdout
#SBATCH --output=%A_%a.stderr
#SBATCH --array=[0-99]%5

ml load Anaconda3

# run all the roadmap files

FILES=(/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/architecture_coordinates/trimmed/*.bed)

echo $SLURM_ARRAY_TASK_ID
echo ${FILES[$SLURM_ARRAY_TASK_ID]}
python /dors/capra_lab/users/fongsl/enh_age/bin/GWAS_enrichment_simple_complex_matched.py ${FILES[$SLURM_ARRAY_TASK_ID]} -g ${1}