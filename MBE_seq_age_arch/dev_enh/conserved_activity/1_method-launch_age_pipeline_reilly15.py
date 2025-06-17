import os
import sys, traceback
import argparse
import datetime
import glob
from functools import reduce
from itertools import groupby
from functools import partial
from multiprocessing import Pool
from pybedtools import BedTool
from pybedtools.helpers import BEDToolsError, cleanup, get_tempdir, set_tempdir
import subprocess

print("last run", datetime.datetime.now())

#%% FUNÇÕES


def sbatch_age_sequence(file, run_any, iterations,\
 run_age, run_break, run_tfbs, run_shuffle, run_testbed):
    path = "/dors/capra_lab/users/fongsl/enh_age/enh_age_git/bin/"
    cmd = "sbatch %sage_enhancers_w_parallelbreaks.slurm %s %s %s %s %s %s %s" \
    % (path, file, iterations, run_age,\
     run_break, run_tfbs, run_shuffle, run_testbed)

    print("SLURM", file, datetime.datetime.now())
    print(cmd)
    if run_any ==1:
        subprocess.check_call(cmd, shell=True)

def python_age_sequence(file, run_any, iterations, run_age,\
 run_break, run_tfbs, run_shuffle, run_testbed):
    path = "/dors/capra_lab/users/fongsl/enh_age/enh_age_git/bin/"
    cmd = '''python %sage_enhancers_w_parallelbreaks.py %s -i %s -a %s -b %s -t %s -sh %s -rt %s'''\
    % (path, file, iterations, run_age, run_break, run_tfbs, run_shuffle, run_testbed)
    print(cmd)
    if run_any ==1:
        subprocess.check_call(cmd, shell=True)


# %% markdown
# # launch REILLY15 w/ python in cmd line


#%%  set paths

# scripts for pipeline
os.chdir("/dors/capra_lab/users/fongsl/enh_age/enh_age_git/bin/")

# data path
source_path = "/dors/capra_lab/projects/enhancer_ages/reilly15/data/"
samples = glob.glob("%sMmu*_brain_enhancers_reilly15_gapexcluded.liftOver.to.Hg19.bed" % source_path)


#%%


sample_dict = {}

ITERATIONS = 10
AGE_VAL = 1
BREAK_VAL = 1
TFBS_VAL = 0
SHUF_VAL = 0
RUN_BED = 1

RUN = 0 # Launch command or dont
SBATCH = 0  # Sbatch or run w/ python interpreter
val = 0
for sample in samples:
    sample_id = (sample.split("/")[-1]).split(".")[0]
    if "chrX" not in sample_id and "chrY" not in sample_id:
        sample_dict[sample_id] = sample
print(sample_dict.keys(), len(sample_dict.keys()))
#%%
print("start", datetime.datetime.now())

for sample_id, file in sample_dict.items():

    if SBATCH ==1:
         sbatch_age_sequence(file, RUN, ITERATIONS, AGE_VAL, BREAK_VAL, TFBS_VAL, SHUF_VAL, RUN_BED)
    else:
        python_age_sequence(file, RUN, ITERATIONS, AGE_VAL, BREAK_VAL, TFBS_VAL, SHUF_VAL, RUN_BED)

    val +=1

print("end", datetime.datetime.now())

#%%
