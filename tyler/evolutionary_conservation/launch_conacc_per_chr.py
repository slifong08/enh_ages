import glob
import numpy as np
import os, sys
import subprocess


# add the path with the file containing all the paths to other files.
CONFIG_PATH = '/dors/capra_lab/users/fongsl/enh_ages/tyler/evolutionary_conservation'
sys.path.append(CONFIG_PATH)

import config # import the config.py file with all the files for tyler's project

all_ = config.all
all_ = "/dors/capra_lab/users/fongsl/tyler/data/CON_ACC/prune/GG-LL_all_OCRs.bed"


#%% function

def make_chr_list():
    n = list(np.arange(1, 23))
    #n.append("X")

    chr_list = []
    for num in n:
        chrn = "chr" + str(num)
        chr_list.append(chrn)

    return chr_list


def run_conacc_slurm(config_path, file, branches, msaway, mod):

    # make the command
    cmd = f"sbatch {config_path}/conacc_per_chr.slurm {file} {branches} {msaway} {mod}"

    # tell us what is being run
    print("running CONACC slurm job on", file,
    "\nbranches:", branches, "\nmsa:", msaway,
    "\nmod:", mod)
    print(cmd)
    # run it
    subprocess.call(cmd, shell = True)

def run_conacc_py(config_path, file, branches, msaway, mod):

    # make the command
    cmd = f"python {config_path}/conacc_per_chr.py {file} -br {branches} -msa {msaway} -mod {mod}"

    # tell us what is being run
    print("running CONACC python job on", file, "\nbranches:", branches, "\nmsa:", msaway, "\nmod:", mod)
    print(cmd)
    # run it
    #subprocess.call(cmd, shell = True)

#%%

BRANCHES = ["hg38", "rheMac8", 'hg38-rheMac8']
MSAWAY = "30"
MODEL = 'hg38-rheMac8'

PATH = "/".join(all_.split("/")[:-1]) + "/" # the path
for BRANCH in BRANCHES:
    run_conacc_slurm(CONFIG_PATH, all_, BRANCH, MSAWAY, MODEL)
#run_conacc_py(CONFIG_PATH, file, BRANCHES, MSAWAY, MODEL)
