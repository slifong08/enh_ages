#!/usr/bin/env python
# coding: utf-8

#%% In[1]:


import glob
from itertools import product
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os, sys
from scipy import stats
import seaborn as sns
import statsmodels
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

RE ="/dors/capra_lab/projects/enhancer_ages/roadmap_encode/results/for_publication/arch/"

if os.path.exists(RE) == False:
    os.system("mkdir %s" % RE)


#%% In[3]:


syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t') # read the file
syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2", "mya", "mya2"]] # whittle down the df


#%%
# # load the file

path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19/download/h3k27ac_plus_h3k4me3_minus_peaks/all_roadmap_enh/trimmed/trimmed-310-all_roadmap_enh/"


#%% In[6]: format the dataframe
 def format_df(df):

    enhdf = df.loc[df.chr_s != "chrX"] # autosomes only

    enhdf = enhdf.drop( ["mrca_2"], axis = 1) # bad columns

    com = enhdf.loc[enhdf.core_remodeling == 1] # complex enhancer df
    sim = enhdf.loc[enhdf.core_remodeling == 0] # simple enhancer df


    core = com.groupby("enh_id")[["mrca", "seg_index"]].max().reset_index() # get core ages, maxseg number
    core.columns = ["enh_id", "mrcacore", "max_seg"] # rename cols

    com = pd.merge(com, core, how = "left") # merge w/ copmlex enhancer dataframe

    return com, sim, core

# Function to bin enhancers into 100 bins

def bin_enh(test, enh_id, bins):

# split enhancer topology up into bins, return ages.

    test["len"] = test.end_e - test.start_e
    test["midpoint"] = (test.start_e + (test.len)/2).round(0).map(int)

    start, stop = test.start_e.iloc[0], test.end_e.iloc[0]
    one_percent_step_size = (test.len/bins).unique()

    one_percent_series = np.arange(start, stop, step = one_percent_step_size, dtype = float)
    binpercent = np.arange(start = 0, stop = bins)
    one_percent_series = one_percent_series[:bins] # map the coordinates as vector of step_size_bins
    binpercent = binpercent[:bins]# map the bins as vector of bins

    # create a dataframe of the percent coordinates
    newdf = pd.DataFrame({"coorPercent": one_percent_series, "bin":binpercent})

    newdf["enh_id"] = test.enh_id.unique().item()

    newdf = pd.merge(newdf, test, how = "left", on = "enh_id")


    # map MRCA to bins
    newdf["mrcapercent"] = "none"

    newdf.loc[(newdf.start_s <= newdf.coorPercent.astype(int))\
    & (newdf.end_s >= newdf.coorPercent.astype(int)), "mrcapercent"] = newdf.mrca

    # keep only bins w/mrca assigned to them
    return_df = newdf[["enh_id", "coorPercent","bin",\
     "mrcapercent"]].loc[(newdf.mrcapercent!= "none")] #

    return_df.mrcapercent = return_df.mrcapercent.astype(float)
    return_df.sort_values(by = "bin")

    midpoint = len(return_df)/2
    print(midpoint)
    if midpoint ==50:
        return_df["bin2"] = np.arange(-50,50,1) # for plotting

    else:
        midpoint = round(midpoint,0)
        min, max = int(midpoint-50), int(midpoint+50)
        print(min, max)
        return_df = return_df.iloc[min:max]
        return_df["bin2"] = np.arange(-50, 50,1) # for plotting

    return return_df


#%%
# function to plot mrca landscape

def plot_mrca(col2, Scol2):

    for n, mrca in enumerate(col2.mrca_2.unique()):
        fig, ax = plt.subplots(figsize=(6,6))
        # WORK HERE TO  assign axis FOR FUNCTION
        plot = col2.loc[col2.mrca_2 == mrca]
        Splot = Scol2.loc[Scol2.mrca_2 == mrca] # compare with the shuffle df
        sns.lineplot(x = "bin", y = "mrcapercent", data = plot.loc[plot.bin2>-51],
        color = "green", label = "enh")
        sns.lineplot(x = "bin", y = "mrcapercent", data = Splot.loc[Splot.bin2>-51],
        color = "grey", label = "shuf")
        ax.axvline(75, linewidth=2, color='k')
        ax.axvline(25, linewidth=2, color='k')
        ax.set_ylabel("mean mrca")
        sns.set(font_scale=2)
        sns.set_style("white")
        sns.set("poster")

    # MWU of outer v. inner 50%
        m, mp = stats.mannwhitneyu(plot.mrcapercent.loc[(plot.bin>=25)&(plot.bin<75)],
                                   plot.mrcapercent.loc[(plot.bin>=75)|(plot.bin<25)])
        sm, smp = stats.mannwhitneyu(Splot.mrcapercent.loc[(Splot.bin>=25)&(Splot.bin<75)],
                                  Splot.mrcapercent.loc[(Splot.bin>=75)|(Splot.bin<25)])
        # annotate the graph w/ quartile ages
        first_quartile = plot.mrcapercent.loc[(plot.bin<25)].mean().round(3)
        inter_quartile =  plot.mrcapercent.loc[(plot.bin>=25)&(plot.bin<75)].mean().round(3)
        last_quartile = plot.mrcapercent.loc[(plot.bin>=75)].mean().round(3)
        ax.text(0.15, 0.05, first_quartile, size=24,\
                va="center", ha="center",transform=ax.transAxes)
        ax.text(0.5, 0.05, inter_quartile, size=24,\
                va="center", ha="center",transform=ax.transAxes)
        ax.text(0.85, 0.05, last_quartile, size=24,\
                va="center", ha="center",transform=ax.transAxes)

        sfirst_quartile = Splot.mrcapercent.loc[(Splot.bin<25)].mean().round(3)
        sinter_quartile =  Splot.mrcapercent.loc[(Splot.bin>=25)&(Splot.bin<75)].mean().round(3)
        slast_quartile = Splot.mrcapercent.loc[(Splot.bin>=75)].mean().round(3)
        ax.text(0.15, 0.12, sfirst_quartile, size=24,\
                va="center", ha="center",transform=ax.transAxes,  color = "grey")
        ax.text(0.5, 0.12, sinter_quartile, size=24,\
                va="center", ha="center",transform=ax.transAxes,  color = "grey")
        ax.text(0.85, 0.12, slast_quartile, size=24,\
                va="center", ha="center",transform=ax.transAxes, color = "grey")


        n = len(plot.enh_id.unique())
        Sn = len(Splot.enh_id.unique())

        ax.set(xlabel = "Normalized Enhancer Bin\nn = %s\nEnh mwu-p=%s\
        \nShuf n= %s, Shuf mwu-p=%s"\
         % (n, round(mp, 3), Sn, round(smp, 3)), title = plot.taxon2.unique().item())

        ax.legend().remove()

        plt.savefig("%sSupFig2.15B_complex_mrca_2_%s_fantom_enh_shuf.pdf" % (RE, mrca), bbox_inches = "tight")

    return fig

#%% evaluate enhancers

columns = ["chr_enh", "start_enh", "end_enh", "enh_id", "shuf_id", "seg_index",
"core_remodeling", "arch", "mrca", "taxon", "mrca_2", "taxon2", "mya", "mya2"]

enh = "%sbreaks/trimmed-310-all_roadmap_enh_enh_age_arch_summary_matrix.bed" % path

shuf = "%sshuffle/breaks/shuf-trimmed-310-all_roadmap_enh_enh_age_arch_summary_matrix.bed" % path


shuffle = pd.read_csv(shuf, sep = '\t', header = None) # open shuffle file

shuffle.columns = columns # rename columns

print("shuffle relative simple cut off", shuffle.seg_index.median())

shuffle.loc[shuffle.seg_index == 1, "core_remodeling"] = 0 # redo the core_remodeling
shuffle.loc[shuffle.seg_index == 1, "arch"] = "simple"
shuffle["datatype"] = "SHUFFLE"
shuffle["enh_len"] = shuffle.end_enh - shuffle.start_enh

shuffle.mrca_2 = shuffle.mrca_2.round(3) # round values


df = pd.read_csv(enh, sep = '\t', header = None) # open enhancer file

df.columns = columns # rename columns

print("enhancer relative simple cut off", df.seg_index.median())

df.loc[df.seg_index == 1, "core_remodeling"] = 0

df.loc[df.seg_index == 1, "arch"] = "simple"

df.mrca_2 = df.mrca_2.round(3) # round values

df["enh_len"] = df.end_enh - df.start_enh
df["datatype"] = "roadmap_310"

#%%


com, sim, core = format_df(df)


#%%
collection_dict = {}

for enh_id in com.enh_id.unique():
    test = com.loc[com.enh_id == enh_id]
    newdf  = bin_enh(test, enh_id, 100)
    collection_dict[enh_id] = newdf

#%%
col = pd.concat(collection_dict.values())

col.head()

#%% In[14]:

col.mrcapercent = col.mrcapercent.round(3)
core.mrcacore = core.mrcacore.round(3)

col = pd.merge(col, core, how = "left")
col2 = pd.merge(col, syn_gen_bkgd, how = "left", right_on = "mrca", left_on = "mrcacore")
col = col.sort_values(by = "mrcacore")
col2 = col2.sort_values(by = "mrcacore")

col2.head()

#%%
fig, ax = plt.subplots(figsize = (8,8))
sns.set_context("poster")
#sns.set_style("white")
sns.lineplot(x = "bin", y = "mrcapercent", data = col2, )#hue = "mrca_2")
ax.set(ylabel= "mean mrca", xlabel = "Normalized Enhancer Bin",
 title = "Composite Complex Enhancer per MRCA")

ax.legend(bbox_to_anchor=(1.3, 1))
#plt.savefig("%sfig2_complex_ALLmrca_2_fantom_enh_3plus.pdf" % RE, bbox_inches = "tight")

#%%
fig, ax = plt.subplots(figsize = (8,8))
sns.set_context("poster")
#sns.set_style("white")
sns.lineplot(x = "bin", y = "mrcapercent", data = col2, hue = "mrca_2")
ax.set(ylabel= "mean mrca", xlabel = "Normalized Enhancer Bin",
 title = "Composite Complex Enhancer per MRCA")

ax.legend(bbox_to_anchor=(1.3, 1))
plt.savefig("%sfig2_complex_ALLmrca_2_fantom_enh.pdf" % RE, bbox_inches = "tight")

#%% evaluate shuffled landscape
# "S" for Shuffle
shuf_dict = {}
Spath = "/dors/capra_lab/projects/enhancer_ages/fantom/data/shuffle/breaks/"
Sfs = glob.glob("%s*_age_breaks.bed" % Spath)
print(len(Sfs))

#%%
for Sf in Sfs:
    Sdf = pd.read_csv(Sf, sep ='\t', header = None) # read file

    # name columns
    Sdf.columns = ["chr_s", "start_s", "end_s", "enh_id", "chr_e", "start_e", "end_e",
                 "seg_index", "core_remodeling", "core", "mrca"]
    ids = Sdf.enh_id.drop_duplicates().sample(frac = 0.05).to_list()
    Sdf = Sdf.loc[Sdf.enh_id.isin(ids)]
    Sdf.mrca = Sdf.mrca.round(3)
    Sdf = pd.merge(Sdf, syn_gen_bkgd[["mrca", "mrca_2", "taxon2"]], how = "left")

    Sdf["code"] = "complex_core"
    Sdf.loc[(Sdf.core_remodeling ==0), "code"] = "simple"
    Sdf.loc[(Sdf.core == 0), "code"] = "derived"
    shuf_dict[Sf] = Sdf
#%%
Sdf = pd.concat(shuf_dict.values())
Scom, Ssim, Score = format_df(Sdf)
print(len(Scom))
#%%

#%%
Scom.sort_values(by = ["chr_s", "start_s"]).head()
#%% sample some enhancer ids

ids = Scom.enh_id.drop_duplicates().sample(frac = 0.30)
print(len(ids))
Scollection_dict = {}
#%%
len(col2.enh_id.unique())
#%%

for enh_id in ids:
    test = Scom.loc[Scom.enh_id == enh_id].copy()
    newdf = bin_enh(test, enh_id, 100)
    Scollection_dict[enh_id] = newdf


#%%
Scol = pd.concat(Scollection_dict.values())
#%%

Scol.mrcapercent = Scol.mrcapercent.round(3)
Score.mrcacore = Score.mrcacore.round(3)

Scol = pd.merge(Scol, Score, how = "left")
Scol2 = pd.merge(Scol, syn_gen_bkgd, how = "left", right_on = "mrca", left_on = "mrcacore")
Scol = Scol.sort_values(by = "mrcacore")
Scol2 = Scol2.sort_values(by = "mrcacore")
#%%
Scol2.sort_values(by =["enh_id", "bin"]).head()
#%%

fig, ax = plt.subplots(figsize = (8,8))

sns.lineplot(x = "bin", y = "mrcapercent", data = col.loc[col.bin2>-51],\
color = "green", label = "FANTOM")
sns.lineplot(x = "bin", y = "mrcapercent", data = Scol.loc[Scol.bin2>-51], \
color = "grey", label = "Shuffle")

sns.set(font_scale=2)
sns.set_style("white")

ax.axvline(75, linewidth=2, color='k')
ax.axvline(25, linewidth=2, color='k')
ax.set_ylabel("mean mrca")

# MWU of outer v. inner 50%
m, mp = stats.mannwhitneyu(col.mrcapercent.loc[(col.bin>=25)&(col.bin<75)],
                           col.mrcapercent.loc[(col.bin>=75)|(col.bin<25)])
sm, smp = stats.mannwhitneyu(Scol.mrcapercent.loc[(Scol.bin>=25)&(Scol.bin<75)],
                          Scol.mrcapercent.loc[(Scol.bin>=75)|(Scol.bin<25)])
# annotate the graph w/ quartile ages
first_quartile = col.mrcapercent.loc[(col.bin<25)].mean().round(3)
inter_quartile =  col.mrcapercent.loc[(col.bin>=25)&(col.bin<75)].mean().round(3)
last_quartile = col.mrcapercent.loc[(col.bin>=75)].mean().round(3)
ax.text(0.15, 0.05, first_quartile, size=24,\
        va="center", ha="center",transform=ax.transAxes)
ax.text(0.5, 0.05, inter_quartile, size=24,\
        va="center", ha="center",transform=ax.transAxes)
ax.text(0.85, 0.05, last_quartile, size=24,\
        va="center", ha="center",transform=ax.transAxes)

sfirst_quartile = Scol.mrcapercent.loc[(Scol.bin<25)].mean().round(3)
sinter_quartile =  Scol.mrcapercent.loc[(Scol.bin>=25)&(Scol.bin<75)].mean().round(3)
slast_quartile = Scol.mrcapercent.loc[(Scol.bin>=75)].mean().round(3)
ax.text(0.15, 0.12, sfirst_quartile, size=24,\
        va="center", ha="center",transform=ax.transAxes,  color = "grey")
ax.text(0.5, 0.12, sinter_quartile, size=24,\
        va="center", ha="center",transform=ax.transAxes,  color = "grey")
ax.text(0.85, 0.12, slast_quartile, size=24,\
        va="center", ha="center",transform=ax.transAxes, color = "grey")
n = len(col.enh_id.unique())
Sn = len(Scol.enh_id.unique())

ax.set(xlabel = "Normalized Enhancer Bin\nn = %s\nEnh mwu-p=%s\
\nShuf n= %s, Shuf mwu-p=%s"\
 % (n, round(mp, 3), Sn, round(smp, 3)), title = "Composite Complex  Enhancer")

ax.legend(bbox_to_anchor = (1,1))


plt.savefig("%sfig2b_complex_fantom_enh_shuf.pdf" % RE, bbox_inches = "tight")

#%%
plot_mrca(col2, Scol2)
