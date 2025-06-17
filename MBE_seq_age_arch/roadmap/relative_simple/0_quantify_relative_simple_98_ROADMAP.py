#!/usr/bin/env python
# coding: utf-8

# In[1]:


import glob
import numpy as np
import pandas as pd
import os, sys
import datetime
import seaborn as sns
import matplotlib.pyplot as plt
last_run = datetime.datetime.now()
today = (datetime.date.today())
RE = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/results/for_publication/relative_simple/"

print("last run", datetime.datetime.now())

GENIC = 1

base_path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/"

if GENIC == 1:
    path = "%sh3k27ac_plus_h3k4me3_minus_peaks/breaks/" % base_path
    samples = glob.glob("%sROADMAP_E*_enh_summary_matrix.bed" % path)
elif GENIC == 0:
    path = "%sh3k27ac_plus_h3k4me3_minus_peaks/non-genic/breaks/" % base_path
    samples = glob.glob("%sno-exon_E*_parallel_breaks_enh_age_arch_summary_matrix.bed" % path)
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
     "sid", "enh_len"])[["seg_index"]].max().reset_index()
        enh["sid2"] = df.sid.apply(lambda x: x.split("_")[0])

    else:
        enh = df.groupby(["enh_id", "core_remodeling",\
     "shuf_id", "enh_len"])[["seg_index"]].max().reset_index()
        enh["sid2"] = df.shuf_id.apply(lambda x: x.split("_")[0])

    # add sample_id


        # calculate number of breaks as percentile
    enh["pct"] = enh["seg_index"].rank(pct = True)
    enh.pct = enh.pct.round(2) # round percentile to nearest 10%
    # reduce the dataframe
    pctBreaks = enh.groupby(["sid2", "pct"])["seg_index"].min().reset_index()

    pctBreaks.sort_values(by = "seg_index") # sort the new dataframe.

    return enh, pctBreaks

#%% ## RUN FUNCTION ##


pct_dict = {} # dictionary to collect percentile info
cols = ["enh_id", "shuf_id", "core_remodeling","seg_index",  "enh_len",]

for sample in samples:

    if GENIC ==1:
        sid = (sample.split("/")[-1]).split("_")[1] # get sid

        col_check = pd.read_csv(sample, sep = '\t', header =None, nrows = 1)
        if ":" in col_check[3]:

            df = pd.read_csv(sample, sep = '\t',  header =None, low_memory=False, usecols = [3,4,5,7,9])
            df = df[[4,3,5,7,9]]
            df.columns =cols
            print(df.shuf_id.unique())
        else:
            df = pd.read_csv(sample, sep = '\t',  header =None, low_memory=False, usecols = [16,4,5,7,9])
            df = df[[4,16,5,7,9]]

            df.columns =cols

    elif GENIC ==0:
        sid = (sample.split("/")[-1]).split("_")[1] # get sid
        #outenh = "%sno_exon_%s_enh_percentile_breaks.bed" %(path, sid)
        df = pd.read_csv(sample, sep = '\t',
        low_memory=False, header = None, usecols = [3, 5, 6, 8, 10 ])
        df.columns =cols
    print(sid)

    enh, pctBreaks = get_percentile(df) # get percentile info

    pct_dict[sid] = pctBreaks

    #enh = pd.merge(enh, pctBreaks)

    #enh.to_csv(outenh, sep = '\t', index = False, header = True)


#%%
pct_dict.values()
#%%
pct = pd.concat(pct_dict.values()) # concat all the percentiles

if GENIC ==1:
    outpct = "%sall_ROADMAP_breaks_percentiles.bed" % path
elif GENIC ==0:
    outpct = "%sall_noexon_ROADMAP_breaks_percentiles.bed" % path
#%%
len(pct)
pct.head(50)
#%%
pct.to_csv(outpct, sep ='\t', header = True, index = False) # expore all percentiles.
#%%
len(pct)
#%%
pct = pct.loc[~pct.sid2.str.contains("shuf")]
#%%
pct.head()
#%%
pct.loc[pct.pct<0.5].groupby("sid2")["pct", "seg_index"].max().reset_index()
#%%
len(pct.sid2.unique())
#%%
plot = pct.loc[pct.pct<=0.5].groupby("sid2")["pct", "seg_index"].max().reset_index()
fig, (ax, ax2) = plt.subplots(ncols =2, figsize = (14,6))
sns.distplot(plot.seg_index, ax = ax)

ax.set(title = "median break distribution\n98 ROADMAP samples",\
xlabel = "median breaks in dataset", ylabel = "frequency of datasets",
xlim = (0,10))

sns.boxplot(x = "seg_index", y= "pct", data = plot, ax = ax2)
counts = plot.groupby("seg_index")["pct"].count()
ax2.set(title = "median break distribution\n98 ROADMAP samples",\
xlabel = "median breaks in dataset\n(agebreaks, n)\n%s" % counts, ylabel = "dataset percentile",
)
if GENIC ==1:
    plt.savefig("%srelative_simple_break_pct.pdf" %RE)
elif GENIC ==0:
    plt.savefig("%srelative_simple_no_exon_break_pct.pdf" %RE, bbox_inches = 'tight')

#%%
counts = plot.groupby("seg_index")["pct"].count()

#%%
plot.seg_index.mean()
#%%
sns.lineplot(x = "seg_index", y = "pct", data = pct, hue = "sid2")
