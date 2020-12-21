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
#%%
path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/download/"
fs = glob.glob("%s*/*.bed" % path)
print(len(fs))

iters, age, breaks, tfbs, shuf, run_test = 100, 1, 1, 0, 1, 1
run = 1


#%%

missing = ["UBERON_0000955", "CL_0000767", "CL_0000451", "CL_0000359", "CL_0000576", "UBERON_0000945"]

for f in fs:
    fid = ("_".join((f.split("/")[-1]).split("_")[:2]))
    if fid in missing:
        print(f)

        run_python(f, iters, age, breaks, tfbs, shuf, run_test, run)
