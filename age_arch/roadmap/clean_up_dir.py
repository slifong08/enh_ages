import os
import glob
import subprocess

path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19/download/h3k27ac_plus_h3k4me3_minus_peaks/"
fs = glob.glob("%sHsap_H3K27ac_plus_H3K4me3_minus_E*/cleaned_Hsap_H3K27ac_plus_H3K4me3_minus_E*/shuffle/breakstrimmed/" % path)

print(len(fs))
#%%
fs
#%%
for f in fs:
    #path = "/".join(f.split("/")[:-2])
    cmd = "rm -r %s" % (f)

    subprocess.call(cmd, shell = True)

    
