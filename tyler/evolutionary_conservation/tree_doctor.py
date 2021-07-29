import glob
import os
import subprocess
import numpy as np

#%%
msa = 30


msaway = str(msa) + "way"

# the neutral tree
mod_path = f"/dors/capra_lab/data/ucsc/hg38/multiz{msaway}"
mod = f"{mod_path}/hg38.phastCons{msaway}.mod"
new_mod = f"{mod_path}/hg38.phastCons{msaway}_hg38-rheMac8.mod"

PHAST_PATH = "/dors/capra_lab/bin/"
treedoctor = f"{PHAST_PATH}./tree_doctor"

ape_list = [
"panTro5","panPan2", "gorGor5","ponAbe2","nomLeu3", # ape branches
"macFas5", "macNem1", "cerAty1", "papAnu3", "chlSab2", "manLeu1",
"nasLar1", "colAng1", "rhiRox1", "rhiBie1"
]

prunelist = ",".join(ape_list)
cmd = f"{treedoctor} -p {prunelist} {mod} > {new_mod}"
subprocess.call(cmd, shell = True)


#%%
new_mod
mod
