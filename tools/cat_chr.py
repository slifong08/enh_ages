import glob
import os, sys
import numpy as np
import subprocess

def make_chr_list():
    n = list(np.arange(1, 23))
    n.append("X")

    chr_list = []
    for num in n:
        chrn = "chr" + str(num)
        chr_list.append(chrn)

    return chr_list

def cat_cmd(chrn):

    path = os.getcwd() # get the path.

    outfile = f"{path}/syn_{chrn}.bed" # outfile
    fs = glob.glob(f"{chrn}:*.bed")

    if len(fs)>0:
        print(f"concatenating {chrn}")
        cmd = f"cat {chrn}:*.bed > {outfile}" # cat by chromosome
        subprocess.call(cmd, shell = True)

        if os.path.getsize(outfile) >0:
            rm = f"{chrn}:*.bed" # clean up
            os.remove(rm)
            print(f"removed files {chrn}")

            return outfile

        else:
            return None
    else:
        print(f"nothing to concatenate {chrn}")
        return None

def final_cat():
    outfile = "syn_age_breaks.bed"
    cmd = f"cat syn_chr*.bed > {outfile}" # cat by chromosome
    subprocess.call(cmd, shell = True)

    if os.path.getsize(outfile) >0:
        rm = f"syn_chr*.bed" # clean up
        os.remove(rm)
        print(f"removed syn files")

        return outfile


#%%
chr_list = make_chr_list()

chr_flist = []

for chrn in chr_list:
    chr_f = cat_cmd(chrn)
    chr_flist.append(chr_f)
