#!/usr/bin/env python
# coding: utf-8

# In[1]:


import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os, sys
from scipy import stats
import seaborn as sns
import statsmodels
get_ipython().run_line_magic('matplotlib', 'inline')

from functools import partial
from multiprocessing import Pool
from pybedtools import BedTool
import subprocess
get_ipython().run_line_magic('matplotlib', 'inline')
# sns colors
arch_colors = ["amber", "dusty purple", "windows blue","greyish"]
arch_palette = sns.xkcd_palette(arch_colors)
sns.palplot(arch_palette)

colors = [ "amber", "faded green", "dusty purple", "windows blue","greyish"]
palette = sns.xkcd_palette(colors)
sns.palplot(palette)

cs = ["faded green", "greyish"]
cs_pal = sns.xkcd_palette(cs)
sns.palplot(cs_pal)

es_colors = [ "slate grey","greyish"]
es_pal = sns.xkcd_palette(es_colors)
sns.palplot(es_pal)

# sns graphing preferences
sns.set(color_codes=True)
sns.set(font_scale=1.5)
sns.set_style("white")
sns.despine(bottom=True, left=True)

import datetime
last_run = datetime.datetime.now()
today = (datetime.date.today())
print("last run", datetime.datetime.now())

RE ="/dors/capra_lab/projects/enhancer_ages/fantom/results/arch/"


# In[2]:


syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t') # read the file
syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2", "mya", "mya2"]] # whittle down the df
syn_gen_bkgd


# In[3]:


path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/te/"
f = "%sall_unique_fantom_erna_112_tissue_x_te.bed" % path


# In[4]:


df = pd.read_csv(f, sep ='\t', header = None)

df.columns = ["chr_s", "start_s", "end_s", "enh_id", "chr_e", "start_e", "end_e",
             "seg_index", "core_remodeling", "core", "mrca", "mrca_2", "code", "syn_id",
             "chr_t", "start_t", "end_t", "te", "te_overlap", "mystery"]
df = df.loc[df.chr_s != "chrX"]

df = df.drop( ["mrca_2", "code", "mystery"], axis = 1) # bad columns
df["syn_len"] = df.end_s -df.start_s
df = df.loc[df.syn_len>5]
df.head()


# In[5]:


core = df.groupby("enh_id")["mrca"].max().reset_index()
core.columns = ["enh_id", "mrcacore"]

enh_ids = core.sample(frac = 1) # randomly sample 1% of dataframe.

enh_wTE = df.groupby("enh_id")["te_overlap"].max().reset_index()

enh_wTE = enh_wTE.loc[enh_wTE.te_overlap !="."] # select only enh_ids w/ TE overlap
enh_wTE = enh_wTE.drop(["te_overlap"], axis = 1) # drop the TE overlap column to merge w/ main dataframe

dfTE = pd.merge(enh_wTE, df, how = "left", on = "enh_id") # merge main dataframe against only enh w/ TE
dfTE.head()


# # FIGURE OUT HOW TO EVALUATE TE LOCATION

# In[6]:


def bin_enh(df, bins, id_dict, iteration):
# split enhancer topology up into bins, return ages.
    enh_id = id_dict[iteration]

    test = df.loc[df.enh_id == enh_id].copy()

    test["len"] =test.end_e - test.start_e # enhancer length
    test["midpoint"] = (test.start_e + (test.len)/2).round(0).map(int) # calculate the midpoint
    test.head()

    start, stop = test.start_e.unique().item(), test.end_e.unique().item()
    percent = (test.len/bins).unique() # estimate how many bp per bin

    one_percent_series = np.arange(start, stop, step = percent, dtype = float) #create series of bin coordinates
    binpercent = np.arange(start = 0, stop = bins) # create a series of bins
    one_percent_series = one_percent_series[:bins] # map bins onto enhancer
    binpercent = binpercent[:bins]  # make sure there are only bins 1-99 for 100-binnings

    test = test[["enh_id", "syn_id", "start_s", "end_s","mrca", "start_t", "end_t", "te", "te_overlap", "core_remodeling"]].drop_duplicates()
    newdf = pd.DataFrame({"binCoor": one_percent_series, "bin":binpercent}) # create a dataframe of the percents

    newdf["enh_id"] = test.enh_id.unique().item()

    newdf = pd.merge(newdf, test, how = "left") # merge the binned coordinates w/ dataframe

    # assign mrca to bin
    newdf["mrcapercent"] = "none"
    newdf.loc[(newdf.start_s <= newdf.binCoor)                           & (newdf.end_s > newdf.binCoor), "mrcapercent"] = newdf.mrca
    # ID TE overlaps in bins

    newdf["te_bin"] = -1
    newdf.loc[(newdf.start_t <= newdf.binCoor)                      & (newdf.end_t > newdf.binCoor), "te_bin"] = newdf.bin

    TEdf = newdf.loc[(newdf.te_bin!=-1)&(newdf.mrcapercent.map(str) !="none" )] # make df w/ only the TEs
    TEdf = TEdf.loc[TEdf.te_overlap.map(int)>=6] # erase TFs with less than 6bp overlap
    TEdf=TEdf[['bin', 'enh_id', 'start_t', 'end_t', 'te', 'te_overlap',
 'mrcapercent', 'te_bin']] # reduce df

    return_df = newdf.loc[newdf["mrcapercent"].astype(str) != "none", ["enh_id", "binCoor","bin", "mrcapercent", "core_remodeling"]]

    return_df.mrcapercent = return_df.mrcapercent.astype(float)
    return_df = pd.merge(return_df, TEdf, how = "left").drop_duplicates() # merge bin df + TE df

    return return_df


# In[7]:


# bin enhancer landscapes

RUN_BINNING = 0

if RUN_BINNING ==1:
    id_dict = {}

    for n, enh_id in enumerate(dfTE.enh_id.unique()):
        id_dict[n] = enh_id


    collection_dict = {}
    te_dict = {}
    val =0

    end_val = len(dfTE.enh_id.unique())
    num_threads = 100
    bins = 100

    print(val, end_val)

    while val <  end_val:

        if end_val - val > num_threads:

            new_range = np.arange(num_threads) + val # figure out the number of threads to assign

            print(val)

        else: # unless this is the last batch of threads, which will be less than a thread request.
            num_threads = abs(end_val - val)

            new_range = np.arange(num_threads) + val

            print(val)
            print(end_val,"minus", val, "equals", num_threads)


        pool = Pool(num_threads)
        partial_bin = partial(bin_enh, dfTE, bins, id_dict)
        results = pool.map(partial_bin, [i for i in new_range])
        pool.close()
        pool.join()

        temp = pd.concat(results)

        collection_dict[val] = temp

        val +=num_threads
    col = pd.concat(collection_dict.values())
    col.to_csv("%sfantom_age_complexity_arch_te_only.tsv" % path, sep = '\t', index = False)

else:
    col = pd.read_csv("%sfantom_age_complexity_arch_te_only.tsv" % path, sep = '\t')


# In[8]:


col.mrcapercent = col.mrcapercent.map(float).round(3)
core.mrcacore = core.mrcacore.round(3)

col = pd.merge(col, core, how = "left")
col2 = pd.merge(col, syn_gen_bkgd, how = "left", right_on = "mrca", left_on = "mrcacore")

col.head()


# # complex landscape

# In[9]:


# COMPLEX ENH groupby core_age and #te in bin
complx = col2.loc[col2.core_remodeling ==1].copy()
complx_mrca_te_counts = complx.groupby(["mrca_2","te_bin"])["te"].count().reset_index()
complx_mrca_te_counts.columns =["mrca_2", "bin", "te_count"]


# In[10]:


for mrca in complx.mrca_2.unique():
    print(mrca)

    test = complx.loc[(complx.mrca_2 == mrca) & (complx.te_bin>-1)] # mrca_2 is based off the mrca_core age

    test_te = complx_mrca_te_counts.loc[(complx_mrca_te_counts.mrca_2 == mrca)&(complx_mrca_te_counts.bin>-1)]
    test_te["te_count_z"] = stats.zscore(test_te.te_count)

    fig, ax = plt.subplots(figsize = (8,8))

    plot_ = sns.barplot(x = "bin", y = "te_count_z", data = test_te, color= "green")
    sns.set(font_scale=2)
    ax.set_ylabel("composite TE count z-score")

    ax.set_xlabel("normalized enhancer coordinates\nn = %s" % (len(test.enh_id.unique())))

    taxon = test.taxon2.unique().item()
    ax.set_title(taxon)

    plt.axhline(0, color = "k")

    # show only some of the xticklabels

    for ind, label in enumerate(plot_.get_xticklabels()):
        if ind % 20 == 0:  # every 10th label is kept
            label.set_visible(True)
        else:
            label.set_visible(False)

    plt.savefig("%sfigS5-2_complex_mrca2_%s_fantom_enh_TE_overlap.pdf" % (RE, mrca), bbox_inches = "tight")


# In[11]:


sns.set("poster")
complx_mrca_te_counts = complx.groupby(["te_bin"])["te"].count().reset_index()
complx_mrca_te_counts.columns =[ "bin", "te_count"]

test = complx.loc[(complx.te_bin>-1)]

test_te = complx_mrca_te_counts.loc[(complx_mrca_te_counts.bin>-1)]
test_te["te_count_z"] = stats.zscore(test_te.te_count)

fig, ax = plt.subplots(figsize = (8,8))

sns.lineplot(x = "bin", y = "mrcapercent", data = test)

ax.set_ylabel("mean composite mrca")

ax.set_xlabel("n = %s" % len(test.enh_id.unique()))

ax2 = ax.twinx()

sns.lineplot(x = "bin", y = "te_count_z", data = test_te, color= "green")

ax2.set_ylabel("composite TE count z-score")
# show only some of the xticklabels

for ind, label in enumerate(plot_.get_xticklabels()):
    if ind % 20 == 0:  # every 10th label is kept
        label.set_visible(True)
        print(ind)
    else:
        label.set_visible(False)

plt.savefig("%scomplex_mrca2_ALL_fantom_enh_TE_overlap.pdf" % (RE), bbox_inches = "tight")


# In[12]:


test_te.head()


# In[18]:


m, p = stats.mannwhitneyu(test_te.loc[(test_te.bin<=25)| (test_te.bin> 75), "te_count_z"],
                  test_te.loc[(test_te.bin>25) & (test_te.bin <=75), "te_count_z"])
print(m, p)


# In[19]:


test_te.loc[(test_te.bin<=25)| (test_te.bin> 75), "te_count_z"].median()

test_te.loc[(test_te.bin>25) & (test_te.bin <=75), "te_count_z"].median()


# In[13]:


test_te.bin = test_te.bin.astype(int)

fig, ax = plt.subplots(figsize = (8,8))

plot_ = sns.barplot(x = "bin", y = "te_count_z", data = test_te, color= "green")

ax.set_xlabel("Normalized Enhancer Coordinates\nn = %s" % len(test.enh_id.unique()))

ax.set_ylabel("composite TE count z-score")

plt.axhline(0, color = "k")
ax.set_ylim(-2, 3.7)
# show only some of the xticklabels

for ind, label in enumerate(plot_.get_xticklabels()):
    if ind % 20 == 0:  # every 10th label is kept
        label.set_visible(True)
        print(ind)
    else:
        label.set_visible(False)
plt.savefig("%sfig5b-complex_mrca2_ALL_fantom_enh_TE_ONLY_barplot.pdf" % (RE), bbox_inches = "tight")


# # SINE/ALU in complex

# In[17]:


col2.head()


# In[21]:


# COMPLEX ENH groupby core_age and #te in bin
SINEdf = col2.loc[(col2.core_remodeling ==1) & col.te.str.contains("SINE/Alu")].copy()
alu_complx_mrca_te_counts = SINEdf.groupby(["mrca_2","te_bin"])["te"].count().reset_index()
alu_complx_mrca_te_counts.columns =["mrca_2", "bin", "te_count"]

alu_complx_mrca_te_counts.head()


# In[23]:


alu_complx_mrca_te_counts.mrca_2.unique()


# In[24]:


for mrca in alu_complx_mrca_te_counts.mrca_2.unique():
    test = SINEdf.loc[(SINEdf.mrca_2 == mrca) & (SINEdf.te_bin>-1)] # mrca_2 is based off the mrca_core age

    test_te = alu_complx_mrca_te_counts.loc[(alu_complx_mrca_te_counts.mrca_2 == mrca)&(alu_complx_mrca_te_counts.bin>-1)]
    test_te["te_count_z"] = stats.zscore(test_te.te_count)

    fig, ax = plt.subplots(figsize = (8,8))

    plot_ = sns.barplot(x = "bin", y = "te_count_z", data = test_te, color= "green")
    sns.set(font_scale=2)
    ax.set_ylabel("composite TE count z-score")

    ax.set_xlabel("normalized enhancer coordinates\nn = %s" % (len(test.enh_id.unique())))

    taxon = test.taxon2.unique().item()
    ax.set_title(taxon)

    plt.axhline(0, color = "k")

    # show only some of the xticklabels

    for ind, label in enumerate(plot_.get_xticklabels()):
        if ind % 20 == 0:  # every 10th label is kept
            label.set_visible(True)
        else:
            label.set_visible(False)

#    plt.savefig("%sfigS5-2_complex_mrca2_%s_fantom_enh_TE_overlap.pdf" % (RE, mrca), bbox_inches = "tight")


# # simple landscape

# In[21]:


# SIMPLE ENH groupby core_age and #te in bin
simpl = col2.loc[col2.core_remodeling ==0].copy()
simpl_mrca_te_counts = simpl.groupby(["mrca_2","te_bin"])["te"].count().reset_index()
simpl_mrca_te_counts.columns =["mrca_2", "bin", "te_count"]

for mrca in simpl.mrca_2.unique():
    print(mrca)

    test = simpl.loc[(simpl.mrca_2 == mrca) & (simpl.te_bin>-1)]

    test_te = simpl_mrca_te_counts.loc[(simpl_mrca_te_counts.mrca_2 == mrca)&(simpl_mrca_te_counts.bin>-1)]
    test_te["te_count_z"] = stats.zscore(test_te.te_count)

    fig, ax = plt.subplots(figsize = (8,8))
    plot_ = sns.barplot(x = "bin", y = "te_count_z", data = test_te, color= "gold")
    sns.set(font_scale=2)
    ax.set_ylabel("composite TE count z-score")

    ax.set_xlabel("normalized enhancer coordinates\nn = %s" % (len(test.enh_id.unique())))

    taxon = test.taxon2.unique().item()
    ax.set_title(taxon)
    plt.axhline(0, color = "k")

    # show only some of the xticklabels

    for ind, label in enumerate(plot_.get_xticklabels()):
        if ind % 20 == 0:  # every 10th label is kept
            label.set_visible(True)
            print(ind)
        else:
            label.set_visible(False)
    plt.savefig("%ssimple_mrca2_%s_fantom_enh_TE_overlap.pdf" % (RE, mrca), bbox_inches = "tight")


# In[15]:


sns.set("poster")
simpl_mrca_te_counts = simpl.groupby(["te_bin"])["te"].count().reset_index()
simpl_mrca_te_counts.columns =[ "bin", "te_count"]

test = simpl.loc[(simpl.te_bin>-1)]

test_te = simpl_mrca_te_counts.loc[(simpl_mrca_te_counts.bin>-1)]
test_te["te_count_z"] = stats.zscore(test_te.te_count)

fig, ax = plt.subplots(figsize = (8,8))

plot_ = sns.barplot(x = "bin", y = "te_count_z", data = test_te, color= "gold")

ax.set_xlabel("Normalized Enhancer Coordinates\nn = %s" % len(test.enh_id.unique()))

ax.set_ylabel("composite TE count z-score")
ax.set_ylim(-2, 3.7)
plt.axhline(0, color = "k")

# show only some of the xticklabels

for ind, label in enumerate(plot_.get_xticklabels()):
    if ind % 20 == 0:  # every 10th label is kept
        label.set_visible(True)
        print(ind)
    else:
        label.set_visible(False)
plt.savefig("%sfig5b-simple_mrca2_ALL_fantom_enh_TE_ONLY_barplot.pdf" % (RE), bbox_inches = "tight")


# In[16]:


fig, ax = plt.subplots(figsize = (8,8))

plot_ = sns.barplot(x = "bin", y = "te_count", data = test_te, color= "gold")

ax.set_xlabel("Normalized Enhancer Coordinates\nn = %s" % len(test.enh_id.unique()))

ax.set_ylabel("composite TE counts in bins")
plt.axhline(0)

# show only some of the xticklabels

for ind, label in enumerate(plot_.get_xticklabels()):
    if ind % 20 == 0:  # every 10th label is kept
        label.set_visible(True)
        print(ind)
    else:
        label.set_visible(False)
plt.savefig("%sfig5b-complex_mrca2_ALL_fantom_enh_TE_ONLY_counts_bar.pdf" % (RE), bbox_inches = "tight")


# In[22]:


m, p = stats.mannwhitneyu(test_te.loc[(test_te.bin<=25)| (test_te.bin> 75), "te_count_z"],
                  test_te.loc[(test_te.bin>25) & (test_te.bin <=75), "te_count_z"])
print(m, p)


# In[28]:


test_te.loc[(test_te.bin<=25)| (test_te.bin> 75), "te_count_z"].median()


# In[29]:


test_te.loc[(test_te.bin>25) & (test_te.bin <=75), "te_count_z"].median()


# In[ ]:
