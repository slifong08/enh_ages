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
RE = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/"
print("last run", datetime.datetime.now())

path = "%sh3k27ac_plus_h3k4me3_minus_peaks/breaks/" % RE
samples = glob.glob("%sE*_age_breaks.bed" % path)


# In[3]:


desc_file = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/roadmap_hg19_sample_id_desc.csv"
desc_df= pd.read_csv(desc_file, sep = ',', header = None)
desc_df.columns = ["sid", "desc"]


### FUNCTION ###


def get_percentile(samplef):
    df = pd.read_csv(samplef, sep = '\t',  low_memory=False)

    # reduce dataframe to enhancer architecture information
    if "sid" in list(df):
        enh = df.groupby(["enh_id", "core_remodeling",\
     "sid", "enh_len"])["seg_index", "mrca"].max().reset_index()
        enh["sid2"] = df.sid.apply(lambda x: x.split("_")[0])

    else:
        enh = df.groupby(["enh_id", "core_remodeling",\
     "shuf_id", "enh_len"])["seg_index", "mrca"].max().reset_index()
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


    sid = (sample.split("/")[-1]).split("_")[0] # get sid
    print(sid)

    if sid != "E003": # DEAL WITH THIS PROBLEM CHILD LATER.

        enh, pctBreaks = get_percentile(sample) # get percentile info


        pct_dict[sid] = pctBreaks

        enh = pd.merge(enh, pctBreaks)
        outenh = "%s%s_enh_breaks.bed" %(path, sid)
        enh.to_csv(outenh, sep = '\t', index = False, header = True)


pct = pd.concat(pct_dict.values()) # concat all the percentiles

outpct = "%sall_ROADMAP_breaks_percentiles.bed" % path

pct.to_csv(outpct, sep ='\t', header = True, index = False) # expore all percentiles.
#%%
print(sample)
df = pd.read_csv(sample, sep = '\t',  low_memory=False)
df.head()
#%%
len(pct.sid2.unique())
fig, ax = plt.subplots()
sns.distplot(pct.loc[pct.pct == 0.5, "seg_index"])

ax.set(title = "median distribution")

#%%
pct.loc[pct.sid2 == "E100"]
