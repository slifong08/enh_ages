#!/bin/bash
#SBATCH --mail-user=sarah.l.fong@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1    # comments allowed
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
#SBATCH --mem=64GB
#SBATCH --output=/dors/capra_lab/users/fongsl/enh_age/

#SBATCH --job-name=1kg_intersection
#SBATCH --account=capra_lab
#SBATCH --error=%J.stdout
#SBATCH --output=%J.stderr


ml load Anaconda3

# ${1} = bed file
# ${2} = outdir
# ${3} = "-a <<age enhancers>>" # str
# ${4} = "-b <<assemble breaks enhancers>>" # str
# ${5} = "-t <<tfbs density enhancers>>" # str
# ${6} = "-sh <<shuffle, age, breaks, tfbs_den enhancers>>" # str
# ${7} = "-rb <<runbed, age, breaks, tfbs_den enhancers>>" # str

echo python /dors/capra_lab/users/fongsl/enh_ages/core_v_derived/1000G/1kg_intersection.py -b ${1} -o ${2} 

python /dors/capra_lab/users/fongsl/enh_ages/core_v_derived/1000G/1kg_intersection.py -b ${1} -o ${2} 
