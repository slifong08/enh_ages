import os, sys
import subprocess
import glob

def run_python(infile, iters, age, breaks, tfbs, shuf, run_test, run):

    cmd = "python /dors/capra_lab/users/fongsl/enh_age/enh_age_git/bin/awk_breaks.py\
    %s -i %d -a %d -b %d -t %d -sh %d -rt %d" % (infile, iters, age, breaks, tfbs,
    shuf, run_test)
    print(cmd)

    if run ==1:
        subprocess.call(cmd, shell = True)

def run_sbatch(infile, iters, age, breaks, tfbs, shuf, run_test, run):
    cmd = "sbatch /dors/capra_lab/users/fongsl/enh_age/enh_age_git/bin/awk_breaks.slurm %s %d %d %d %d %d %d" \
    % (infile, iters, age, breaks, tfbs, shuf, run_test)

    print("SLURM")
    print(cmd)
    if run ==1:
        subprocess.call(cmd, shell=True)
#%%
path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/"
fs = glob.glob("%s**/cleaned_Hsap_H3K27ac_plus_H3K4me3_minus_*.bed" % path)


print(len(fs))

iters, age, breaks, tfbs, shuf, run_test = 100, 1, 1, 0, 1, 1
run = 1
sbatch = 1


missing = ['E013', 'E011', 'E015', 'E012', 'E004', 'E008', 'E073']

#%%


for f in fs:

    fid = ((f.split("/")[-1]).split("_")[-1]).split(".")[0]

    if fid in missing:

        print(fid)

        if sbatch ==0:
            run_python(f, iters, age, breaks, tfbs, shuf, run_test, run)

        elif sbatch ==1:
            run_sbatch(f, iters, age, breaks, tfbs, shuf, run_test, run)
