import glob
import os, sys
import subprocess

PATH = "/dors/capra_lab/data/tf_motif/vierstra_motif_archetypes/"

FS = glob.glob(f"{PATH}*.bed.gz" )
FS
#%%
for f in FS:
    new_f = "".join(f.split("_temp")[:])
    cmd = f"mv {f} {new_f}"
    subprocess.call(cmd, shell = True)

    
