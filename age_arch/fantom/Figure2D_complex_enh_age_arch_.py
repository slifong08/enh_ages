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

RE ="/dors/capra_lab/projects/enhancer_ages/fantom/results/for_publication/arch/"

if os.path.exists(RE) == False:
    os.system("mkdir %s" % RE)


#%% In[3]:


syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t') # read the file
syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2", "mya", "mya2"]] # whittle down the df
syn_gen_bkgd


#%%
# # load the file


path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/"
f = "%sall_unique_fantom_erna_112_tissue.bed" % path


#%% In[6]: format the dataframe
 def format_df(df):

    enhdf = df.loc[df.chr_s != "chrX"] # autosomes only

    enhdf = enhdf.drop( ["mrca_2", "code"], axis = 1) # bad columns

    com = enhdf.loc[enhdf.core_remodeling == 1] # complex enhancer df
    sim = enhdf.loc[enhdf.core_remodeling == 0] # simple enhancer df


    core = com.groupby("enh_id")["mrca"].max().reset_index() # get core ages
    core.columns = ["enh_id", "mrcacore"]

    com = pd.merge(com, core, how = "left") # merge w/ copmlex enhancer dataframe

    return com, sim, core

# Function to bin enhancers into 100 bins

def bin_enh(df, enh_id, bins):
# split enhancer topology up into bins, return ages.

    test = df.loc[df.enh_id == enh_id].copy()

    test["len"] = test.end_e - test.start_e
    test["midpoint"] = (test.start_e + (test.len)/2).round(0).map(int)

    start, stop = test.start_e.unique().item(), test.end_e.unique().item()
    percent = (test.len/bins).unique()

    one_percent_series = np.arange(start, stop, step = percent, dtype = float)
    binpercent = np.arange(start = 0, stop = bins)
    one_percent_series = one_percent_series[:bins]
    binpercent = binpercent[:bins]

    # create a dataframe of the percent coordinates
    newdf = pd.DataFrame({"coorPercent": one_percent_series, "bin":binpercent})

    newdf["enh_id"] = test.enh_id.unique().item()

    newdf = pd.merge(newdf, test, how = "left")

    newdf["mrcapercent"] = "none"
    newdf.loc[(newdf.start_s <= newdf.coorPercent)\
            & (newdf.end_s >= newdf.coorPercent), "mrcapercent"] = newdf.mrca


    return_df = newdf[["enh_id", "coorPercent","bin",\
     "mrcapercent"]].loc[(newdf.mrcapercent!= "none")]

    return_df.mrcapercent = return_df.mrcapercent.astype(float)
    return_df.sort_values(by = "bin")

    if len(return_df) ==100:
        return_df["bin2"] = np.arange(-50,50,1) # for plotting

    else:
        print(len(return_df))
        return_df["bin2"] = np.arange(-51,50,1) # for plotting
    return return_df
    #return return_df

#%%
# function to plot mrca landscape

def plot_mrca(col2, Scol2):

# initialize N
    N = 2
    # All possible N combination tuples
    # Using list comprehesion + product()
    res = [ele for ele in product(range(0, N + 1), repeat = 2)]

    fig, axes = plt.subplots(nrows=3, ncols=3)

    for n, mrca in enumerate(col2.mrca_2.unique()):

        # WORK HERE TO  assign axis FOR FUNCTION
        plot = col2.loc[col2.mrca_2 == mrca]
        Splot = Scol2.loc[Scol2.mrca_2 == mrca] # compare with the shuffle df
        sns.lineplot(x = "bin2", y = "mrcapercent", data = plot.loc[plot.bin2>-51],
        color = "green", label = "enh")
        sns.lineplot(x = "bin2", y = "mrcapercent", data = Splot.loc[Splot.bin2>-51],
        color = "grey", label = "shuf")
        ax.axvline(-25, linewidth=2, color='k')
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
        ax.text(0.15, 0.05, first_quartile, size=30,\
                va="center", ha="center",transform=ax.transAxes)
        ax.text(0.5, 0.05, inter_quartile, size=30,\
                va="center", ha="center",transform=ax.transAxes)
        ax.text(0.85, 0.05, last_quartile, size=30,\
                va="center", ha="center",transform=ax.transAxes)
        n = len(plot.enh_id.unique())
        Sn = len(Splot.enh_id.unique())

        ax.set(xlabel = "Normalized Enhancer Bin\nn = %s\nEnh mwu-p=%s\
        \nShuf n= %s, Shuf mwu-p=%s"\
         % (n, round(mp, 3), Sn, round(smp, 3)), title = plot.taxon2.unique().item())

        ax.legend(bbox_to_anchor = (1,1))

    plt.savefig("%sfigS2.2_complex_mrca_2_%s_fantom_enh.pdf" % (RE, mrca), bbox_inches = "tight")

    return fig

#%% evaluate enhancers
df = pd.read_csv(f, sep ='\t', header = None) # read file

# name columns
df.columns = ["chr_s", "start_s", "end_s", "enh_id", "chr_e", "start_e", "end_e",
             "seg_index", "core_remodeling", "core", "mrca", "mrca_2", "code", "syn_id"]

#%%
com, sim, core = format_df(df)
#%%
collection_dict = {}

for enh_id in com.enh_id.unique():

    newdf = bin_enh(com, enh_id, 100)
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
sns.lineplot(x = "bin", y = "mrcapercent", data = col2, hue = "mrca_2")
ax.set(ylabel= "mean mrca", xlabel = "Normalized Enhancer Bin",
 title = "Composite Complex Enhancer per MRCA")

ax.legend(bbox_to_anchor=(1.3, 1))
plt.savefig("%sfig2_complex_ALLmrca_2_fantom_enh.pdf" % RE, bbox_inches = "tight")

#%% evaluate shuffled landscape
# "S" for Shuffle
Spath = "/dors/capra_lab/projects/enhancer_ages/fantom/data/shuffle/breaks/"
Sfs = glob.glob("%s*_age_breaks.bed" % Spath)
Sf = Sfs[0]
Sdf = pd.read_csv(Sf, sep ='\t', header = None) # read file

# name columns
Sdf.columns = ["chr_s", "start_s", "end_s", "enh_id", "chr_e", "start_e", "end_e",
             "seg_index", "core_remodeling", "core", "mrca"]
Sdf.mrca = Sdf.mrca.round(3)
Sdf = pd.merge(Sdf, syn_gen_bkgd[["mrca", "mrca_2", "taxon2"]], how = "left")

Sdf["code"] = "complex_core"
Sdf.loc[(Sdf.core_remodeling ==0), "code"] = "simple"
Sdf.loc[(Sdf.core == 0), "code"] = "derived"

Scom, Ssim, Score = format_df(Sdf)

#%%
Scollection_dict = {}

for enh_id in Scom.enh_id.unique():

    newdf = bin_enh(Scom, enh_id, 100)
    Scollection_dict[enh_id] = newdf

Scol = pd.concat(Scollection_dict.values())
#%%
Scol.mrcapercent = Scol.mrcapercent.round(3)
Score.mrcacore = Score.mrcacore.round(3)

Scol = pd.merge(Scol, Score, how = "left")
Scol2 = pd.merge(Scol, syn_gen_bkgd, how = "left", right_on = "mrca", left_on = "mrcacore")
Scol = Scol.sort_values(by = "mrcacore")
Scol2 = Scol2.sort_values(by = "mrcacore")

#%% In[13]:

fig, ax = plt.subplots(figsize = (8,8))

sns.lineplot(x = "bin2", y = "mrcapercent", data = col.loc[col.bin2>-51],\
color = "green", label = "FANTOM")
sns.lineplot(x = "bin2", y = "mrcapercent", data = Scol.loc[Scol.bin2>-51], \
color = "grey", label = "Shuffle")

sns.set(font_scale=2)
sns.set_style("white")

ax.axvline(-25, linewidth=2, color='k')
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
ax.text(0.15, 0.05, first_quartile, size=30,\
        va="center", ha="center",transform=ax.transAxes)
ax.text(0.5, 0.05, inter_quartile, size=30,\
        va="center", ha="center",transform=ax.transAxes)
ax.text(0.85, 0.05, last_quartile, size=30,\
        va="center", ha="center",transform=ax.transAxes)
n = len(col.enh_id.unique())
Sn = len(Scol.enh_id.unique())

ax.set(xlabel = "Normalized Enhancer Bin\nn = %s\nEnh mwu-p=%s\
\nShuf n= %s, Shuf mwu-p=%s"\
 % (n, round(mp, 3), Sn, round(smp, 3)), title = "Composite Complex  Enhancer")

ax.legend(bbox_to_anchor = (1,1))


plt.savefig("%sfig2b_complex_fantom_enh.pdf" % RE, bbox_inches = "tight")

#%%
# plot per age



#%% simple enhancer analysis

sim_dict = {}
for enh_id in sim.enh_id.unique():

    newdf = bin_enh(sim, enh_id, 100)
    sim_dict[enh_id] = newdf


#%% In[20]:


col_sim = pd.concat(sim_dict.values())
col_sim.mrcapercent = col_sim.mrcapercent.round(3)
core.mrcacore = core.mrcacore.round(3)

col_sim = pd.merge(col_sim, core, how = "left", on = "enh_id")
col_sim.head()


#%% In[23]:


fig, ax = plt.subplots(figsize = (8,8))
sns.lineplot(x = "bin", y = "mrcapercent", data = col_sim)
ax.set_ylabel("mean mrca")
ax.set_xlabel("Normalized Enhancer Bin")
ax.set_title("Composite Simple Enhancer")
ax.legend(bbox_to_anchor = (1.4, 1))
plt.savefig("%ssimple_fantom_enh.pdf" % RE, bbox_inches = "tight")


#%% In[24]:


fig, ax = plt.subplots(figsize = (8,8))
sns.lineplot(x = "bin", y = "mrcapercent", data = col_sim, hue = "mrcacore")
ax.set_ylabel("mean mrca")
ax.set_xlabel("Normalized Enhancer Bin")
ax.set_title("Composite Simple Enhancer Per MRCA")
ax.legend(bbox_to_anchor=(1, 1))
plt.savefig("%ssimple_mrcacore_fantom_enh.pdf" % RE, bbox_inches = "tight")


#%% In[38]:
# # TE analysis


te_path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/te/"
te_f = "%sall_unique_fantom_erna_112_tissue_x_te.bed" % te_path

tedf = pd.read_csv(te_f, sep = "\t", header = None) # read df

tedf.columns = ["chr_s", "start_s", "end_s", "enh_id", "chr_e", "start_e", "end_e",
             "seg_index", "core_remodeling", "core", "mrca", "mrca_2", "code", "syn_id",
             "chr_t", "start_t", "end_t", "te", "te_overlap", "mystery"]

tedf = tedf.loc[tedf.chr_s != "chrX"] # remove chrX
tedf = tedf.drop( ["mrca_2", "code", "mystery"], axis = 1) # bad columns

# create TE binary
tedf["te_bin"] = 0
tedf["te_bin"].loc[tedf.start_t>0] = 1

# get enh_id + TE binary
te_enh = tedf.groupby("enh_id")["te_bin"].max().reset_index()
te_enh.head()


#%% In[39]:


# merge te binary w/ collection dataframes
colte =pd.merge(col, te_enh, how = "left")
col2te =pd.merge(col2, te_enh, how = "left")

col2te.head()


#%% In[287]:


for mrca in col2.mrca_2.unique():
    print(mrca)
    plot = col2te.loc[col2te.mrca_2 == mrca]

    fig, ax = plt.subplots(figsize = (8,8))

    sns.lineplot(x = "bin", y = "mrcapercent", data = plot, hue = "te_bin")
    ax.axvline(25, linewidth=2, color='k')
    ax.axvline(75, linewidth=2, color='k')
    ax.set_ylabel("mrca core age = %s" % mrca)
    ax.set_xlabel("wo/te n = %s, w/te n = %s" % (len(plot.loc[plot.te_bin==0].enh_id.unique()),
                                                 len(plot.loc[plot.te_bin==1].enh_id.unique())))
    ax.legend(bbox_to_anchor=(1, 1))
    plt.savefig("%scomplex_mrca_2_TE_%s_fantom_enh.pdf" % (RE, mrca), bbox_inches = "tight")
