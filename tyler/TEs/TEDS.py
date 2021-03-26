import matplotlib.pyplot as plt
import os, sys
import pandas as pd
import pybedtools as pb
from scipy import stats
import seaborn as sns
import subprocess


PATH = "/dors/capra_lab/users/fongsl/tyler/data/"
FILE = "GG-LL_species-specific_OCRs_rank.bed"
F = os.path.join(PATH, FILE)


TEDSPATH = "/dors/capra_lab/data/transposable_elements/repeatmasker/"
TEDSFILE = "hg38.noheader.bed"
TEDSF = os.path.join(TEDSPATH, TEDSFILE)


#%%

def bed_intersect(f, tedsf, path):

    sid = (f.split("/")[-1]).split(".bed")[0]
    tedsid = (tedsf.split("/")[-1]).split(".")[0]
    outf = f"{path}{sid}_x_{tedsid}.bed"

    cmd = f"bedtools intersect -a {f} -b {tedsf} -wao > {outf}" # write all the overlapping and non-overlapping hars
    print(cmd)
    if os.path.exists(outf) == False:

        subprocess.call(cmd, shell = True)

    return outf

def get_df(outf):
    df = pd.read_csv(outf, sep = '\t', header = None)
    return df
#%%

outf = bed_intersect(F, TEDSF, PATH)
df = get_df(outf)
print(outf)
