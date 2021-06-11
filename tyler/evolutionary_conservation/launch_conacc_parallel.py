import os, sys
import subprocess

# add the path with the file containing all the paths to other files.
PATH = '/dors/capra_lab/users/fongsl/enh_ages/tyler/evolutionary_conservation'
sys.path.append(PATH)

import config # import the config.py file with all the files for tyler's project

all_, hars_, phast_ = config.all, config.hars, config.phastCons
hu_specific_noTE = config.hu_specific_noTE

#%% function


def run_conacc_slurm(path, file, branches, msaway):

    # make the command
    cmd = f"sbatch {path}/conacc_parallel.slurm {file} {branches} {msaway}"

    # tell us what is being run
    print("running CONACC slurm job on", file, "\nbranches:", branches, "\nmsa:", msaway)
    print(cmd)
    # run it
    subprocess.call(cmd, shell = True)

def run_conacc_py(path, file, branches, msaway):

    # make the command
    cmd = f"python {path}/conacc_parallel.py {file} -br {branches} -m {msaway}"

    # tell us what is being run
    print("running CONACC python job on", file, "\nbranches:", branches, "\nmsa:", msaway)
    print(cmd)
    # run it
    #subprocess.call(cmd, shell = True)

#%%

BRANCHES = 'hg38-rheMac8'
MSAWAY = "30"

run_conacc_slurm(PATH, all_, BRANCHES, MSAWAY)
run_conacc_slurm(PATH, hu_specific_noTE, BRANCHES, MSAWAY)
#%%
run_conacc_py(PATH, hu_specific_noTE, BRANCHES, MSAWAY)
