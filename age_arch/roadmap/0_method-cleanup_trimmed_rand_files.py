#%%
# 20200813
# sarahfong
# clean up the rand_files where age breaks have been assembled.

import os
import glob
import subprocess

path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/\
hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks\
/trimmed/shuffle/ages/"

#%%
rand_files = glob.glob("%srand_file_shuf-E*_*.bed" %path)

rand_dict = {} # dict for rand_sid:rand_file

#%%
shuf_dict = {}
shuf_files = glob.glob("%sshuf-trimmed-310_E*_*_parallel_breaks.bed.gz" %path)

for f in shuf_files:
    fid = "_".join(((f.split("/")[-1]).split(".")[0]).split("_")[1:3])
    shuf_dict[fid] = f

#%%
for f in rand_files:
    fid = ((f.split("/")[-1]).split("-")[1]).split(".")[0]
    rand_dict[fid] = f
    if fid not in shuf_dict.keys():
        rename = "%sshuf-trimmed-310_%s_parallel_breaks.bed" %(path, fid)
        cmd = "mv %s %s | gzip" % (f, rename)
        subprocess.call(cmd, shell = True)
