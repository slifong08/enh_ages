

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
print("last run", datetime.datetime.now())

need_shuffles = ["E055", "E111", "E078", "E102", "E103", "E008", "E098", "E012",\
 "E114", "E056", "E105", "E034", "E066", "E029", "E003", "E063"]
print(len(need_shuffles))


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

def breaks_array(search_str):
    slurm_script = "/dors/capra_lab/users/fongsl/enh_age/enh_age_git/bin/age_enhancers_w_parallelbreaks_array.slurm"
    cmd = "sbatch %s %s" % (slurm_script, search_str)
    subprocess.check_call(cmd, shell=True)
# %% markdown
# # launch emera16 w/ python in cmd line


# %% select the file(s) to run
os.chdir("/dors/capra_lab/users/fongsl/enh_age/enh_age_git/bin/")
source_path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/non-genic/"
samples = glob.glob("%sno-exon_E*.bed"%source_path)
#%%

#%%
sample_dict = {}
ITERATIONS = 10
AGE_VAL = 1
BREAK_VAL = 1
TFBS_VAL = 0
SHUF_VAL = 1
RUN_BED = 0

RUN = 1 # Launch command or dont
SBATCH = 1  # Sbatch or run w/ python interpreter
RUN_BREAKS = 0
val = 0
for sample in samples:
    sample_id = (sample.split("/")[-1]).split(".")[0]
    if sample_id.split("_")[1] in need_shuffles:
        sample_dict[sample_id] = sample
print(len(sample_dict.keys()))
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
for sample_id, file in sample_dict.items():
    chr_splitter(file)
#%%
if RUN_BREAKS ==1:
    search_str = "/dors/capra_lab/projects/enhancer_ages/emera16/data/shuffle/shuf*.bed"
    breaks_array(search_str)
