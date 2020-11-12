import glob
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

RE ="/dors/capra_lab/projects/enhancer_ages/fantom/results/for_publication/non-genic/"

#%%
path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/non-genic/"
fs = glob.glob("%sexonOverlap_SHUFFLE_FANTOM_shuf-all_fantom_enh_112_tissues-*_age_breaks_summary_matrix.bed" %path)

f = fs[0]
f
#%%
genic_dict = {}
for i, f in enumerate(fs):
    df = pd.read_csv(f, sep = '\t', header = None, usecols = [3, 4,5,6, 10, 11]).drop_duplicates()
    cols = ["enh_id",  "core_remodeling", "arch", "seg_index", "mrca_2", "taxon"]
    df.columns = cols
    genic_dict[i] = df

df = pd.concat(genic_dict.values())
#%%
total = len(df)
arch_counts = df.groupby("arch")["enh_id"].count().reset_index()
arch_counts.columns = ["arch", "arch_total"]

arch_counts ["freq"] = arch_counts.arch_total.divide(total)
arch_counts
#%%
age_counts = df.groupby("mrca_2")["enh_id"].count().reset_index()
age_counts ["freq"] = age_counts.enh_id.divide(total)
age_counts
#%%
freq = df.groupby(["mrca_2","arch"])["enh_id"].count().reset_index()
freq.columns = ["mrca_2","arch", "age_arch_count"]

freq["total"] = total

freq = pd.merge(freq, arch_counts, how = "left")
freq.head()
#%%
freq["percent_arch"] = freq.age_arch_count.divide(freq.arch_total)
freq["percent_total"] = freq.age_arch_count.divide(freq.total)
freq
#%%
fs_noex = glob.glob("%sno-exon_SHUFFLE_FANTOM_shuf-all_fantom_enh_112_tissues-*_age_breaks_summary_matrix.bed" %path)
nongenic_dict = {}
for i, f in enumerate(fs_noex):
    df = pd.read_csv(f, sep = '\t', header = None, usecols = [3, 4,5,6, 10, 11]).drop_duplicates()
    cols = ["enh_id",  "core_remodeling", "arch", "seg_index", "mrca_2", "taxon"]
    df.columns = cols
    nongenic_dict[i] = df
#%%
df = pd.concat(nongenic_dict.values())
#%%
total = len(df)
arch_counts = df.groupby("arch")["enh_id"].count().reset_index()
arch_counts.columns = ["arch", "arch_total"]

arch_counts ["freq"] = arch_counts.arch_total.divide(total)
arch_counts
#%%
age_counts = df.groupby("mrca_2")["enh_id"].count().reset_index()
age_counts ["freq"] = age_counts.enh_id.divide(total)
age_counts
#%%
freq = df.groupby(["mrca_2","arch"])["enh_id"].count().reset_index()
freq.columns = ["mrca_2","arch", "age_arch_count"]

freq["total"] = total

freq = pd.merge(freq, arch_counts, how = "left")
freq.head()
#%%
freq["percent_arch"] = freq.age_arch_count.divide(freq.arch_total)
freq["percent_total"] = freq.age_arch_count.divide(freq.total)
freq
#%%
100714/(2626874+100714)
#%%
(2626874+100714)
