import glob
import os, sys
import pandas as pd
from scipy import stats

PATH = "/dors/capra_lab/projects/enhancer_ages/encode/data/"
FS = glob.glob(f"{PATH}*genomecov*.bed")
FS
F = FS[0]

def get_genomecov(f):
    cols = ["#chr", "cov_bin", "start", "stop", "percent"]
    df = pd.read_csv(f, sep = '\t', header = None, names = cols)
    print(df.head())

    cov = df.loc[(df.cov_bin ==1) & (df["#chr"]== "genome"), "percent"]
    print(cov)
#%%
FS[0]
get_genomecov(FS[2])
