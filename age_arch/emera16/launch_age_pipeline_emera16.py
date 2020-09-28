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

def sbatch_age_sequence(file, run_any, iterations,\
 run_age, run_break, run_tfbs, run_shuffle, run_testbed):
    path = "/dors/capra_lab/users/fongsl/enh_age/enh_age_git/bin/"
    cmd = "sbatch %sage_enhancers.slurm %s %s %s %s %s %s %s" \
    % (path, file, iterations, run_age,\
     run_break, run_tfbs, run_shuffle, run_testbed)

    print("SLURM", file, datetime.datetime.now())

    if run_any ==1:
        subprocess.check_call(cmd, shell=True)

def python_age_sequence(file, run_any, iterations, run_age,\
 run_break, run_tfbs, run_shuffle, run_testbed):
    path = "/dors/capra_lab/users/fongsl/enh_age/enh_age_git/bin/"
    cmd = '''python %sage_enhancers.py %s -i %s -a %s -b %s -t %s -sh %s -rt %s'''\
    % (path, file, iterations, run_age, run_break, run_tfbs, run_shuffle, run_testbed)
    print(cmd)
    if run_any ==1:
        subprocess.check_call(cmd, shell=True)

# %% markdown
# # launch FANTOM w/ python in cmd line
# %% codecell
os.chdir("/dors/capra_lab/users/fongsl/enh_age/bin/")
source_path = "/dors/capra_lab/data/enhancers/emera16/"

samples = glob.glob("%semera_2016_neocortical_dev_enhancers_hu_ms.bed"%source_path)

sample_dict = {}

ITERATIONS = 100
AGE_VAL = 1
BREAK_VAL = 1
TFBS_VAL = 0
SHUF_VAL = 1
RUN_BED = 1

RUN = 1 # Launch command or dont
SBATCH = 1  # Sbatch or run w/ python interpreter
val = 0
for sample in samples:
    sample_id = (sample.split("/")[-1]).split(".")[0]
    sample_dict[sample_id] = sample

print("start", datetime.datetime.now())

for sample_id, file in sample_dict.items():

    if SBATCH ==1:
         sbatch_age_sequence(file, RUN, ITERATIONS, AGE_VAL, BREAK_VAL, TFBS_VAL, SHUF_VAL, RUN_BED)
    else:
        python_age_sequence(file, RUN, ITERATIONS, AGE_VAL, BREAK_VAL, TFBS_VAL, SHUF_VAL, RUN_BED)

    val +=1

print("end", datetime.datetime.now())
