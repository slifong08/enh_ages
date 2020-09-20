#!/bin/bash
#SBATCH --mail-user=sarah.l.fong@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1    # comments allowed 
#SBATCH --ntasks=1
#SBATCH --time=04:00:00
#SBATCH --mem=20GB
#SBATCH --output=/dors/capra_lab/users/fongsl/enh_age/

#SBATCH --job-name=age_enh
#SBATCH --account=capra_lab
#SBATCH --error=%J.stdout
#SBATCH --output=%J.stderr


ml load Anaconda3

# ${1} = /dors/capra_lab/users/fongsl/enh_age/data/<<FILE_TO_AGE>>.bed # str
# ${2} = "-i <<num shuffle interations>>" # str
# ${3} = "-a <<age enhancers>>" # str
# ${4} = "-b <<assemble breaks enhancers>>" # str
# ${5} = "-t <<tfbs density enhancers>>" # str
# ${6} = "-sh <<shuffle, age, breaks, tfbs_den enhancers>>" # str
# ${7} = "-rb <<runbed, age, breaks, tfbs_den enhancers>>" # str

python /dors/capra_lab/users/fongsl/enh_age/bin/age_enhancers.py ${1} -i ${2} -a ${3} -b ${4} -t ${5} -sh ${6} -rt ${7}
