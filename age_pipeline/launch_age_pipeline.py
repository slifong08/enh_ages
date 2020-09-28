# %% codecell
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
source_path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/download/"

samples = glob.glob("%s*.bed"%source_path)

sample_dict = {}

ITERATIONS = 0
AGE_VAL = 1
BREAK_VAL = 1
TFBS_VAL = 0
SHUF_VAL = 0
RUN_BED = 1

RUN = 1 # Launch command or dont
SBATCH =0  # Sbatch or run w/ python interpreter
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
# %% markdown
# # launch ROADMAP w/ python in cmd line
# %% codecell
os.chdir("/dors/capra_lab/users/fongsl/enh_age/bin/")
source_path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/trimmed/"
#samples = glob.glob("%sHsap_H3K27ac_plus_H3K4me3_minus_E003.bed"%source_path)
samples = glob.glob("%s*.bed"%source_path)
sample_dict = {}

ITERATIONS = 0
AGE_VAL = 1
BREAK_VAL = 1
TFBS_VAL = 0
SHUF_VAL = 0
RUN_BED = 1

RUN = 1 # Launch command or dont
SBATCH = 1 # Sbatch or run w/ python interpreter

for sample in samples:
    sample_id = (sample.split("/")[-1]).split(".")[0]
    sample_dict[sample_id] = sample

print("start", datetime.datetime.now())
val = 0
for sample_id, file in sample_dict.items():
    sid = sample_id.split("_")[-1]

    print(sid, val)

    if SBATCH ==1:
         sbatch_age_sequence(file, RUN, ITERATIONS, AGE_VAL, BREAK_VAL, TFBS_VAL, SHUF_VAL, RUN_BED)
    else:
        python_age_sequence(file, RUN, ITERATIONS, AGE_VAL, BREAK_VAL, TFBS_VAL, SHUF_VAL, RUN_BED)
        print(sample_id)
    val +=1
print("end", datetime.datetime.now())
# %% markdown
# # trimmed_shuffle ROADMAP
# %% codecell
os.chdir("/dors/capra_lab/users/fongsl/enh_age/bin/")
source_path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/trimmed/shuffle/"
#samples = glob.glob("%sHsap_H3K27ac_plus_H3K4me3_minus_E003.bed"%source_path)
samples = glob.glob("%sshuf-trimmed-310*.bed"%source_path)
sample_dict = {}

ITERATIONS = 0
AGE_VAL = 1
BREAK_VAL = 1
TFBS_VAL = 0
SHUF_VAL = 0
RUN_BED = 1

RUN = 1 # Launch command or dont
SBATCH = 1 # Sbatch or run w/ python interpreter

for sample in samples[1:]:
    sample_id = (sample.split("/")[-1]).split(".")[0]
    sample_dict[sample_id] = sample

print("start", datetime.datetime.now())
val = 0
for sample_id, file in sample_dict.items():
    sid = sample_id.split("_")[-1]

    print(sid, val)

    if SBATCH ==1:
         sbatch_age_sequence(file, RUN, ITERATIONS, AGE_VAL, BREAK_VAL, TFBS_VAL, SHUF_VAL, RUN_BED)
    else:
        python_age_sequence(file, RUN, ITERATIONS, AGE_VAL, BREAK_VAL, TFBS_VAL, SHUF_VAL, RUN_BED)
        print(sample_id)
    val +=1

print("end", datetime.datetime.now())

# %% markdown
# # launch Ernst MPRA w/ python in cmd line
# %% codecell
########################
# get MPRA enhancers
########################
enh_path ="/dors/capra_lab/data/mpra/ernst16/"

os.chdir(enh_path)
samples = glob.glob("%s**/tiles_*_ScaleUpDesign1and2_combinedP.bed" % enh_path, recursive = True)

sample_dict = {}
for sample in samples:
    sample_id = (sample.split("/")[-1]).split("_")[1]
    sample_dict[sample_id] = sample

print("start", datetime.datetime.now())
for sample_id, file in sample_dict.items():
    cmd = "python /dors/capra_lab/users/fongsl/enh_age/bin/age_enhancers.py %s %s" % (file, sample_id)
    print(sample_id, datetime.datetime.now())
    os.system(cmd)
print("end", datetime.datetime.now())
# %% markdown
# # launch Klein MPRA w/ python in cmd line
# %% codecell
enh_path ="/dors/capra_lab/projects/enhancer_ages/klein2018/data/"
samples = glob.glob("%smerged.bed"%enh_path) # published and not published enhancers tiles merged together

# %% codecell
enh_path ="/dors/capra_lab/projects/enhancer_ages/klein2018/data/"
samples = glob.glob("%smerged.bed"%enh_path) # published and not published enhancers tiles merged together
samples = glob.glob("%sklein_tile_activity.bed" % enh_path)
print("start", datetime.datetime.now())

ITERATIONS = 100
AGE_VAL = 1
BREAK_VAL = 1
TFBS_VAL = 1
SHUF_VAL = 1
RUN_BED = 0

RUN = 1 # Launch command or dont
SBATCH = 0 # Sbatch or run w/ python interpreter

val = 0
sample_dict={}
for sample in samples:
    sample_id = (sample.split("/")[-1]).split(".")[0]
    sample_dict[sample_id] = sample



for sample_id, file in sample_dict.items():

    if SBATCH ==1:
         sbatch_age_sequence(file, RUN, ITERATIONS, AGE_VAL, BREAK_VAL, TFBS_VAL, SHUF_VAL, RUN_BED)
    else:
        python_age_sequence(file, RUN, ITERATIONS, AGE_VAL, BREAK_VAL, TFBS_VAL, SHUF_VAL, RUN_BED)


        #%run '/dors/capra_lab/users/fongsl/enh_age/bin/age_enhancers.py'\
        #$file '-i' $ITERATIONS "-a" $AGE_VAL "-b" $BREAK_VAL "-t" $TFBS_VAL "-sh" $SHUF_VAL

    val +=1
# %% markdown
# # launch villar
# %% codecell
########################
# get MPRA enhancers
########################
enh_path ="/dors/capra_lab/projects/enhancer_ages/villar15/Hsap/"

os.chdir(enh_path)
file = enh_path ="/dors/capra_lab/projects/enhancer_ages/villar15/data/HSap_H3K4me3.H3K27Ac_overlap_H3K27Aconly.bed"

# %% markdown
# # fimo fantom
# %% codecell
source_path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/download/"
samples = glob.glob("%s*enhancers.bed"%source_path)
sample_dict = {}

RUN = 0 # Launch command?
SBATCH = 1 # Sbatch or run w/ python interpreter

for sample in samples:
    sample_id = ''.join(((sample.split("/")[-1]).split(".")[0]).split("_")[:2])
    sample_dict[sample_id] = sample

for sample_id, file in sample_dict.items():
    if SBATCH ==1:
        cmd = "sbatch /dors/capra_lab/users/fongsl/enh_age/bin/fimo_tfbs.slurm %s %s"\
        % (file, sample_id)

        print("SLURM", sample_id, datetime.datetime.now())
        print(cmd)
        if RUN ==1:
            os.system(cmd)
    else:
        cmd = "python /dors/capra_lab/users/fongsl/enh_age/bin/fimo_tfbs.py %s %s" \
        % (file, sample_id)

        print("PYTHON", sample_id, datetime.datetime.now())
        print(cmd)
        if RUN ==1:
            os.system(cmd)
    break
# %% codecell
