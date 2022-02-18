import glob
import matplotlib.pyplot as plt
import numpy as np
import os, sys
import pandas as pd

import seaborn as sns
from scipy import stats

#%%

ENHPATH = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19/download/h3k27ac_plus_h3k4me3_minus_peaks/"
ENHFILE = "no-exon_all_roadmap_enh_breaks.bed"
ENHF = os.path.join(ENHPATH, "all_roadmap_enh/breaks", ENHFILE)


ENHFS = glob.glob("/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19/download/h3k27ac_plus_h3k4me3_minus_peaks/Hsap_H3K27ac_plus_H3K4me3_minus_E*/non-genic/no-exon_E*_parallel_breaks_enh_age_arch_summary_matrix.bed")



SHUFFILE = ""
SHUFF = os.path.join(ENHPATH, "shuffle/breaks", SHUFFILE)


#%%


def format_df(df):
    cols = ["chr", "start", "end", "enh_id",
    #"reads",
    "id",
    'core_remodeling', 'arch',
    "seg_index",
    #"mrca",
     "enh_len",
     #"taxon",
     "mrca_2", "taxon2"]
     #"mya", "mya2", "density", "id_again"]

    df.columns = cols

    return df



#%%

len_dict ={}
mrca_dict = {}

for enhf in ENHFS:
    enh = pd.read_csv(enhf, sep = '\t', header = None, usecols = [0,1,2,3,5,6,7,8,10,12,13])
    df = format_df(enh)
    lens = df.groupby(["id", "arch"])["enh_len"].median().reset_index()
    mrcas = df.groupby(["id", "arch"])["mrca_2"].median().reset_index()
    sid = df.id.iloc[0]
    len_dict[sid] = lens
    mrca_dict[sid] = mrcas


#%%
lens = pd.concat(len_dict.values())
lens.groupby("arch")["enh_len"].median()
sLens = lens.loc[lens.arch == "simple", "enh_len"]
cLens = lens.loc[lens.arch != "simple", "enh_len"]

stats.mannwhitneyu(sLens, cLens)
"""
arch
complexenh    1960.0 median bp long (median of medians)
simple         795.5

MannwhitneyuResult(statistic=55.0, pvalue=3.101001445036452e-33)
"""


#%%
mrcas = pd.concat(mrca_dict.values())
mrcas.groupby("arch")["mrca_2"].median()

sMrcas = mrcas.loc[mrcas.arch == "simple", "mrca_2"]
cMrcas = mrcas.loc[mrcas.arch != "simple", "mrca_2"]

stats.mannwhitneyu(sMrcas, cMrcas)


'''
arch
complexenh    0.308 median of median complex MRCAs
simple        0.175  median of median simple MRCAs


MannwhitneyuResult(statistic=686.0, pvalue=5.826731641379351e-34)
'''
