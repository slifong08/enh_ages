#!/usr/bin/env python
# coding: utf-8

# In[2]:


import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os, sys
from scipy import stats
import seaborn as sns
import statsmodels

from functools import partial
from multiprocessing import Pool
from pybedtools import BedTool
import subprocess
get_ipython().run_line_magic('matplotlib', 'inline')


import datetime
last_run = datetime.datetime.now()
today = (datetime.date.today())
print("last run", datetime.datetime.now())

RE ="/dors/capra_lab/projects/enhancer_ages/roadmap_encode/results/arch/"

if os.path.exists(RE) == False:
    os.system("mkdir %s" % RE)


# In[3]:


syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t') # read the file
syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2", "mya", "mya2"]] # whittle down the df
syn_gen_bkgd


# In[4]:


path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19/download/h3k27ac_plus_h3k4me3_minus_peaks/Hsap_H3K27ac_plus_H3K4me3_minus_E072/breaks/"
f = "%sE072_uniq_age_breaks.bed" % path # in file

df = pd.read_csv(f, sep ='\t', header = None) # open the dataframe


col_names = ["chr_s", "start_s", "end_s", "enh_id", "chr_e", "start_e", "end_e",
             "seg_index", "core_remodeling", "core", "mrca", "id", "enh_len"] # name columns

df.columns = col_names # name df columns

df = df.loc[df.chr_s != "chrX"] # autosomes only
df = df.loc[df.enh_len < 10000] # enhancers shorter than 10kb


df["syn_len"] = df.end_s -df.start_s # calculate syntenic length

df.head()


# # subtract exon overlapping enhancers.

# In[5]:


noexon_f = "%sno-exon_E072.bed" % path # file w/ non-exon overlapping enhancers

noexon = pd.read_csv(noexon_f, sep = '\t', header = None, usecols = [3]) # open only the enhancer_id column

noexon.columns = ["enh_id"] # rename column
noexon["noexon"] = 1 # add binary for filtering in the next cell.

noexon.head()


# In[6]:


# merge and filter the original dataframe for enhancers that do not overlap exons.

noexdf = pd.merge(df, noexon, how = 'left') # merge

df = noexdf.loc[noexdf.noexon ==1] # filter for non-exon overlapping enhancers.


# # identify relative simple label and relabel

# In[7]:


seg_index = df.groupby("enh_id")["seg_index"].max().reset_index() # groupby enhancer, get max number of segments

relative_simple = seg_index.seg_index.median()

print(relative_simple, "> age segments are defined as relatively simple enhancer architectures" )

seg_index["relative_core_remodeling"] = 0 # create a new binary for simple, complex definition

seg_index.loc[seg_index.seg_index >= relative_simple, "relative_core_remodeling"] = 1 # complex architecture is >= relative_simple age segments

seg_index.head()


# In[8]:


df = pd.merge(df, seg_index[["enh_id", "relative_core_remodeling"]], how = "left") # merge the relative simple definition


# # sample enhancers
# ## running all enhancers would take a way-way long time

# In[9]:


core = df.loc[df.syn_len>5].groupby("enh_id")["mrca", "relative_core_remodeling"].max().reset_index()
core.columns = ["enh_id", "mrcacore", "core_remodeling"]

enh_ids = core.sample(frac = 0.05) # randomly sample 5% of dataframe.


# In[10]:


def bin_enh(df, bins, id_dict, iteration):


    # split enhancer topology up into bins, return ages.

    enh_id = id_dict[iteration] # find the enumerated enh_id

    test = df.loc[df.enh_id == enh_id].copy() # get dataframe for one enhancer


    test["len"] = test.end_e - test.start_e # calculate the length of the enhancer
    test["midpoint"] = (test.start_e + (test.len)/2).round(0).map(int) # find the midpoint genomic coordinate

    start, stop = test.start_e.iloc[0], test.end_e.iloc[0] # get genomic enhancer boundaries
    percent = (test.len/bins).unique() # calculate the genomic length reflecting 1 percent of enhancer sequence.

    one_percent_series = np.arange(start, stop, step = percent, dtype = float) # make as series of 1% genomic coordinates
    binpercent = np.arange(start = 0, stop = bins) # label percent.
    one_percent_series = one_percent_series[:bins]
    binpercent = binpercent[:bins]

    newdf = pd.DataFrame({"coorPercent": one_percent_series,
                          "bin": binpercent}) # create a dataframe of the percents

    newdf["enh_id"] = test.enh_id.iloc[0]

    newdf = pd.merge(newdf, test, how = "left")

    if test.seg_index.max()>0:
        newdf["mrcapercent"] = "none"
        newdf.loc[(newdf.start_s <= newdf.coorPercent)                               & (newdf.end_s >= newdf.coorPercent),"mrcapercent"] = newdf.mrca
    else:
        newdf["mrcapercent"] = test.mrca.iloc[0]


    return_df = newdf.loc[(newdf.mrcapercent!= "none"), ["enh_id", "coorPercent","bin", "mrcapercent"]]

    return_df.mrcapercent = return_df.mrcapercent.astype(float)
    return_df = return_df.groupby(["enh_id", "bin"])[["coorPercent","mrcapercent"]].max().reset_index() # prevent duplicates()
    return return_df


# In[11]:


id_dict = {}
for n, enh_id in enumerate(enh_ids.enh_id.unique()):
    id_dict[n] = enh_id


# In[12]:


NEW_RUN = 1
if NEW_RUN ==1:

    # bin the complex enhancers

    id_dict = {}

    num_enh = len(enh_ids) # how many enhancers are we mapping?

    for n, enh_id in enumerate(enh_ids.enh_id.unique()):
        id_dict[n] = enh_id

    collection_dict = {}

    val = 0
    end_val = len(enh_ids.enh_id.unique())
    num_threads = 100
    bins = 100
    print(val, end_val)

    while val <  end_val:
        if end_val - val > num_threads:
            new_range = np.arange(num_threads) + val
            print(val)
        else:
            num_threads = abs(end_val - val)

            new_range = np.arange(num_threads) + val
            print(val)
            print(end_val, "minus", val, "equals", num_threads)


        pool = Pool(num_threads)
        partial_bin = partial(bin_enh, df, bins, id_dict)
        results = pool.map(partial_bin, [i for i in new_range])
        pool.close()
        pool.join()

        temp = pd.concat(results)
        collection_dict[val] = temp
        val +=num_threads



    col = pd.concat(collection_dict.values())

    col.mrcapercent = col.mrcapercent.round(3)
    core.mrcacore = core.mrcacore.round(3)

    col = pd.merge(col, core, how = "left").drop_duplicates()
    col = col.sort_values(by = "mrcacore")
    col2 = pd.merge(col, syn_gen_bkgd[["mrca", "mrca_2", "taxon2"]], left_on = "mrcapercent", right_on = "mrca")
    outF = "%s%s_E072_roadmap_enh_landscapes.csv" % (path, num_enh)

    #col2.to_csv(outF, sep = '\t', header = True, index = False) # save the new run


else:
    outF = "%s5932_E072_roadmap_enh_landscapes.csv" %path
    col2 = pd.read_csv(outF, sep = '\t')


# In[8]:


#outF = "%s5932_E072_roadmap_enh_landscapes.csv" %path
#col2.to_csv(outF, sep = '\t', header = True, index = False)


# In[13]:


def plot_cm(table, df,  sample):
    enh_count = len(table)
    # setup palette
    pal = sns.color_palette("Blues", len(df.mrca_2.unique()))
    pal_lut = dict(zip(set(df["bin_val"].sort_values(ascending = True)), pal))
    sns.set(font_scale=2)

    cm = sns.clustermap(table, col_cluster = False, cmap = pal, figsize = (10,10),
                    yticklabels=0, xticklabels=10)

    cm.cax.remove()
    for label in df.sort_values(by = "bin_val")["bin_val"].unique():
        cm.ax_col_dendrogram.bar(0, 0, color = pal_lut[label],
                                label=une[label], linewidth=0)
    cm.ax_col_dendrogram.legend(loc="lower center", ncol=2)
    plt.xlabel("normalized enhancer coordinates")
    plt.ylabel("N = %s enhancers" % enh_count)
    plt.savefig("%sSupFig2.3_Roadmap_arch_%s_mrca_%s_enh.pdf"% (RE, sample, enh_count), bbox_inches = "tight")

enu = {}
une = {}
ms = {}
for n, mrca_2 in enumerate(col2.mrca_2.sort_values(ascending = True).unique()):

    enu[mrca_2] = n
    taxon = syn_gen_bkgd.loc[syn_gen_bkgd.mrca_2 == round(mrca_2,3), "taxon2"].iloc[0]

    une[n] = taxon
    ms[mrca_2] = taxon

col2["bin_val"] = 0
for mrca in col2.mrca_2.unique():

    col2.loc[col2.mrca_2 == mrca, "bin_val"] = enu[mrca]

complex_table = col2.loc[col2.core_remodeling ==1].sort_values(by = ["enh_id", "bin"]).pivot(index = "enh_id", columns = "bin", values = "bin_val")

simple_table = col2.loc[col2.core_remodeling ==0].sort_values(by = ["enh_id", "bin"]).pivot(index = "enh_id", columns = "bin", values = "bin_val")
complex_df = col2.loc[col2.core_remodeling ==1]
simple_df = col2.loc[col2.core_remodeling ==0]

plot_cm(complex_table, complex_df, "complex")

plot_cm(simple_table, simple_df, "simple")


# In[20]:


une


# In[ ]:
