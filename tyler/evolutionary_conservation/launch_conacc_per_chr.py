import glob
import numpy as np
import os, sys
import subprocess


# add the path with the file containing all the paths to other files.
CONFIG_PATH = '/dors/capra_lab/users/fongsl/enh_ages/tyler/evolutionary_conservation'
sys.path.append(CONFIG_PATH)

import config # import the config.py file with all the files for tyler's project

all_, hars_, phast_ = config.all, config.hars, config.phastCons
FS = [hars_, all_, phast_] # a list of the files

#%% function

def make_chr_list():
    n = list(np.arange(1, 23))
    #n.append("X")

    chr_list = []
    for num in n:
        chrn = "chr" + str(num)
        chr_list.append(chrn)

    return chr_list


def run_conacc_slurm(config_path, file, branches, msaway):

    # make the command
    cmd = f"sbatch {config_path}/conacc_per_chr.slurm {file} {branches} {msaway}"

    # tell us what is being run
    print("running CONACC slurm job on", file, "\nbranches:", branches, "\nmsa:", msaway)
    print(cmd)
    # run it
    subprocess.call(cmd, shell = True)

def run_conacc_py(config_path, file, branches, msaway):

    # make the command
    cmd = f"python {config_path}/conacc_per_chr.py {file} -br {branches} -m {msaway}"

    # tell us what is being run
    print("running CONACC python job on", file, "\nbranches:", branches, "\nmsa:", msaway)
    print(cmd)
    # run it
    #subprocess.call(cmd, shell = True)

#%%

BRANCHES = 'hg38-rheMac8'
MSAWAY = "30"
chr_list = make_chr_list()
F = phast_

PATH = "/".join(F.split("/")[:-1]) + "/" # the path
for chr_f in chr_list[2:]:
    file = f"{PATH}{chr_f}.bed"
    print(file)
    run_conacc_slurm(CONFIG_PATH, file, BRANCHES, MSAWAY)
    #run_conacc_py(CONFIG_PATH, file, BRANCHES, MSAWAY)
