

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

missing = [ 'E050', 'E108', 'E074', 'E071', 'E019','E041', 'E042', 'E062',
 'E047', 'E029', 'E045', 'E127', 'E069', 'E115', 'E123', 'E032', 'E058', 'E055',
 'E039', 'E044', 'E072', 'E126', 'E048', 'E040', 'E124', 'E122', 'E046', 'E118', 'E116']

# missing = [ 'E050', 'E108', 'E074', 'E071', 'E029',  'E127', 'E069', 'E115', 'E058', 'E055',
 # 'E072',  'E126', 'E124', 'E122',]


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


# %% select the file(s) to run
os.chdir("/dors/capra_lab/users/fongsl/enh_age/enh_age_git/bin/")

source_path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/trimmed/"

samples = glob.glob("%sno-exon_trimmed-310-E*.bed"%source_path)
print(len(samples))
#%%
sample_dict = {}
ITERATIONS = 10
AGE_VAL = 1
BREAK_VAL = 1
TFBS_VAL = 0
SHUF_VAL = 1
RUN_BED = 0

RUN = 1 # Launch command or dont
SBATCH = 0  # Sbatch or run w/ python interpreter
RUN_BREAKS = 0
val = 0
for sample in samples:

    sample_id = ((sample.split("/")[-1]).split(".")[0]).split("-")[-1]
    if sample_id in missing:
        print(sample_id)

        sample_dict[sample_id] = sample
print(len(sample_dict.keys()))
#%%
print("start", datetime.datetime.now())

for sample_id, file in sample_dict.items():

    uncap = file.split(".")[0]
    new_f = uncap +"_enh_ages.bed"
    cmd = "mv %s %s" % (file, new_f)
    #subprocess.call(cmd, shell = True)
    print(sample_id)


    if SBATCH ==1:
        sbatch_age_sequence(file, RUN, ITERATIONS, AGE_VAL, BREAK_VAL, TFBS_VAL, SHUF_VAL, RUN_BED)
    else:
        python_age_sequence(file, RUN, ITERATIONS, AGE_VAL, BREAK_VAL, TFBS_VAL, SHUF_VAL, RUN_BED)
    break
    val +=1

print("end", datetime.datetime.now())

#%%
file
#%%
if RUN_BREAKS ==1:
    search_str = "/dors/capra_lab/projects/enhancer_ages/emera16/data/shuffle/shuf*.bed"
    breaks_array(search_str)
