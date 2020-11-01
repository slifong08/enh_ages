#!/usr/bin/env python
# coding: utf-8

# In[1]:


import glob
import numpy as np
import pandas as pd
import os, sys
import datetime
import seaborn as sns
last_run = datetime.datetime.now()
today = (datetime.date.today())
RE = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/results/for_publication/relative_simple/"

print("last run", datetime.datetime.now())

GENIC = 0

base_path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/"

if GENIC == 1:
    path = "%sh3k27ac_plus_h3k4me3_minus_peaks/breaks/" % base_path
    samples = glob.glob("%sE*_age_breaks.bed" % path)
elif GENIC == 0:
    path = "%sh3k27ac_plus_h3k4me3_minus_peaks/breaks/non-genic/" % base_path
    samples = glob.glob("%sno-exon_ROADMAP_E*_enh_and_shuf_age_arch_summary_matrix.bed" % path)
#%%
len(samples)
#%%


desc_file = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/roadmap_hg19_sample_id_desc.csv"
desc_df= pd.read_csv(desc_file, sep = ',', header = None)
desc_df.columns = ["sid", "desc"]


### FUNCTION ###


def get_percentile(df):


    # reduce dataframe to enhancer architecture information
    if "sid" in list(df):
        enh = df.groupby(["enh_id", "core_remodeling",\
     "sid", "enh_len"])[["seg_index", "mrca"]].max().reset_index()
        enh["sid2"] = df.sid.apply(lambda x: x.split("_")[0])

    else:
        enh = df.groupby(["enh_id", "core_remodeling",\
     "shuf_id", "enh_len"])[["seg_index", "mrca"]].max().reset_index()
        enh["sid2"] = df.shuf_id.apply(lambda x: x.split("_")[0])

    # add sample_id


        # calculate number of breaks as percentile
    enh["pct"] = enh["seg_index"].rank(pct = True)
    enh.pct = enh.pct.round(1) # round percentile to nearest 10%
    # reduce the dataframe
    pctBreaks = enh.groupby(["sid2", "pct"])["seg_index"].min().reset_index()

    pctBreaks.sort_values(by = "seg_index") # sort the new dataframe.

    return enh, pctBreaks

#%% ## RUN FUNCTION ##


pct_dict = {} # dictionary to collect percentile info


for sample in samples:

    if GENIC ==1:
        sid = (sample.split("/")[-1]).split("_")[0] # get sid
        outenh = "%s%s_enh_breaks.bed" %(path, sid)
        df = pd.read_csv(sample, sep = '\t',  low_memory=False)
    elif GENIC ==0:
        sid = (sample.split("/")[-1]).split("_")[2] # get sid
        outenh = "%sno_exon_%s_enh_breaks.bed" %(path, sid)
        df = pd.read_csv(sample, sep = '\t',  low_memory=False, header = None)
        df.columns = ["chr_enh","start_enh","end_enh", "shuf_id", "enh_id",
        "core_remodeling", "arch", "seg_index", "mrca", "enh_len", "taxon",
        "mrca_2", "taxon2", "mya", "mya2", "seg_den", "datatype"]
    print(sid)

    if sid != "E003": # DEAL WITH THIS PROBLEM CHILD LATER.

        enh, pctBreaks = get_percentile(df) # get percentile info

        pct_dict[sid] = pctBreaks

        enh = pd.merge(enh, pctBreaks)

        enh.to_csv(outenh, sep = '\t', index = False, header = True)


pct = pd.concat(pct_dict.values()) # concat all the percentiles

if GENIC ==1:
    outpct = "%sall_ROADMAP_breaks_percentiles.bed" % path
elif GENIC ==0:
    outpct = "%sall_noexon_ROADMAP_breaks_percentiles.bed" % path
#%%
pct.to_csv(outpct, sep ='\t', header = True, index = False) # expore all percentiles.
#%%
len(pct)
#%%
pct = pct.loc[~pct.sid2.str.contains("shuf")]

#%%
len(pct.sid2.unique())
fig, ax = plt.subplots()
sns.distplot(pct.loc[pct.pct == 0.5, "seg_index"])

ax.set(title = "median break distribution\n96/98 ROADMAP samples",\
xlabel = "median breaks in dataset", ylabel = "frequency of datasets",
xlim = (0,10))

if GENIC ==1:
    plt.savefig("%srelative_simple_break_dist.pdf" %RE)
elif GENIC ==0:
    plt.savefig("%srelative_simple_no_exon_break_dist.pdf" %RE)


#%%
pct.loc[pct.sid2 == "E100"]
