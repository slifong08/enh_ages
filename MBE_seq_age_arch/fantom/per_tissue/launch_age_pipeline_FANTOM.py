

import argparse
from chr_splitter import chr_splitter
import datetime
import glob
from functools import reduce
from functools import partial
from itertools import groupby
from multiprocessing import Pool
import os
from pybedtools import BedTool
from pybedtools.helpers import BEDToolsError, cleanup, get_tempdir, set_tempdir
import subprocess
import sys, traceback
import time
print("last run", datetime.datetime.now())



#%% FUNÇÕES
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def sbatch_age_sequence(file, run_any, iterations,\
 run_age, run_break, run_tfbs, run_shuffle, run_testbed):
    path = "/dors/capra_lab/users/fongsl/enh_age/enh_age_git/bin/"

    filelen = file_len(file)
    print("FILE_LEN", filelen)


    time = "02:00:00"



    cmd = "sbatch %sage_enhancers_w_parallelbreaks.slurm %s %s %s %s %s %s %s --time=%s" \
    % (path, file, iterations, run_age,\
     run_break, run_tfbs, run_shuffle, run_testbed, time)

    print("SLURM", file, "\n")
    print(cmd)
    if run_any ==1:
        subprocess.check_call(cmd, shell=True)

def python_age_sequence(file, run_any, iterations, run_age,\
 run_break, run_tfbs, run_shuffle, run_testbed):
    path = "/dors/capra_lab/users/fongsl/enh_age/enh_age_git/bin/"
    cmd = '''python %sage_enhancers_w_parallelbreaks.py %s -i %s -a %s -b %s -t %s -sh %s -rt %s '''\
    % (path, file, iterations, run_age, run_break, run_tfbs, run_shuffle, run_testbed)
    print(cmd)
    if run_any ==1:
        subprocess.check_call(cmd, shell=True)

def breaks_array(search_str):
    slurm_script = "/dors/capra_lab/users/fongsl/enh_age/enh_age_git/bin/age_enhancers_w_parallelbreaks_array.slurm"
    cmd = "sbatch %s %s" % (slurm_script, search_str)
    subprocess.check_call(cmd, shell=True)


# %% select the file(s) to run
os.chdir("/dors/capra_lab/users/fongsl/enh_age/enh_age_git/bin/")

# each individual tissue sample
source_path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/download/"
samples = glob.glob("%sshuffle/ages/*enh_ages.bed"%source_path)

source_path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/"
samples = glob.glob("%sall_fantom_enh.bed"%source_path)
len(samples)

#%%
already_done = []

doneFs = glob.glob("%sshuffle/breaks/*matrix.bed"% source_path)
for f in doneFs:
    fid = "_".join(((f.split("/")[-1]).split("huf-")[1]).split("_")[0:2])
    already_done.append(fid)

print(len(already_done))
already_done
#%%

sample_dict = {}

ITERATIONS = 200
AGE_VAL = 1
BREAK_VAL = 1
TFBS_VAL = 0
SHUF_VAL = 1
RUN_BED = 0

RUN = 1 # Launch command or dont
SBATCH =1  # Sbatch or run w/ python interpreter
RUN_BREAKS = 0
val = 0

for sample in samples:
    sample_id = (sample.split("/")[-1]).split(".")[0]
    fid = "_".join(sample_id.split("_")[0:2])


    if fid not in already_done:
        sample_dict[sample_id] = sample
print(len(sample_dict.keys()))
#%%
print("start", datetime.datetime.now(), "\n")

for sample_id, file in sample_dict.items():
    if "_".join(sample_id.split("_")[:2]) not in already_done:
        print(sample_id,  "\n")

        if SBATCH ==1:
            sbatch_age_sequence(file, RUN, ITERATIONS, AGE_VAL, BREAK_VAL, TFBS_VAL, SHUF_VAL, RUN_BED)
        else:
            python_age_sequence(file, RUN, ITERATIONS, AGE_VAL, BREAK_VAL, TFBS_VAL, SHUF_VAL, RUN_BED)

        val +=1

    #if val %10 == 0:

    #    time.sleep(300)

print("\nend", datetime.datetime.now())

#%%
for sample_id, file in sample_dict.items():
    chr_splitter(file)
#%%
if RUN_BREAKS ==1:
    search_str = "/dors/capra_lab/projects/enhancer_ages/emera16/data/shuffle/shuf*.bed"
    breaks_array(search_str)
