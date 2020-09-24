#!/usr/bin/env python
# coding: utf-8
# 20200921
# sarahfong

# evaluate syntenic blocks from simple and complex enhancers.

import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os, sys
from scipy import stats
import seaborn as sns

import datetime
last_run = datetime.datetime.now()
today = (datetime.date.today())
print("last run", datetime.datetime.now())

RE = "/dors/capra_lab/projects/enhancer_ages/results/landscape/fantom/age_arch/"
#%%
# sns colors
# for syntenic architectures
arch_colors = [ "amber", "dusty purple", "windows blue","greyish"]
archpal = sns.xkcd_palette(arch_colors)
sns.palplot(archpal)
archOrder =  ["simple_core", "complex_core", "complex_derived"]

# for shuffled syntenic comparisons
shuf_pal = [ "mustard",  "greyish purple","greyish blue", "greyish"]
spal = sns.xkcd_palette(shuf_pal)
sns.palplot(spal)
spalOrder = ["Shuffle", "FANTOM"]

# for enhancers v. shuffled
colors = ["greyish", "faded green"]
epal = sns.xkcd_palette(colors)
enhOrder = ["shuffle", "FANTOM"]

# sns graphing preferences
sns.set(color_codes=True)
sns.set(font_scale=1.5)
sns.set_style("white")
sns.despine(bottom=True, left=True)

#%% load the genomic background

syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t')
syn_gen_bkgd["mrca"]=syn_gen_bkgd["mrca"].round(3)
syn_gen_bkgd["mrca_2"]=syn_gen_bkgd["mrca_2"].round(3)
syn_gen_bkgd.head()
syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2"]]

# get enhancer file tissue/cell line descriptions

desc_file = "/dors/capra_lab/data/fantom/fantom5/facet_expressed_enhancers/sample_id_descriptions.txt"
desc_df= pd.read_csv(desc_file, sep = '\t', header = None)

#%% define
def format_df(df, syn_gen_bkgd, datasetName):

    # format the dataset columns
    df["datatype"] = datasetName

    # SYNTENIC BLOCKS ONLY
    if "syn_id" in list(df):

        df = df.loc[~df.syn_id.str.contains("chrX")] # remove chrX

        df["code"] = ""
        df.loc[df.core.astype(int) ==0, "code"] = "complex_derived"
        df.loc[df.core_remodeling.astype(int)==0, "code"] = "simple_core"
        df.loc[(df.core.astype(int) ==1)\
         & (df["core_remodeling"].astype(int) ==1), "code"] = "complex_core"

        df = pd.merge(df, syn_gen_bkgd,\
        how = "left", on = "mrca" ).sort_values(by="mrca_2")# add species annotation

        df["seg_index"] =  (df["seg_index"].map(int) + 1)
        # seg_index count reflects # age segments instead of breaks.
        df["seg_rep"] = df.syn_len.astype(int).divide(df.enh_len.astype(int))
        return df.drop_duplicates()

    # ENHANCERS ONLY - max age, max breaks.
    elif "enh_tf_u_count" in list(df):

        df = df.loc[~df.enh_id.str.contains("chrX")] # remove chrX

        df["arch"] = 0
        df.loc[(df["core_remodeling"] ==0), "arch"] = "simple_core"
        df.loc[(df["core_remodeling"] ==1), "arch"] = "complex_core"
        df = pd.merge(df, syn_gen_bkgd,\
         how = "left", on = "mrca" ).sort_values(by="mrca_2")
         # add species annotation

        return df.drop_duplicates()

#%% # run samples

path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/"

samples = "%sall_unique_fantom_erna_112_tissue.bed" % path
samples
#%%

# concatenate the TFBS dataframe
df = pd.read_csv(samples, sep='\t', header = None, \
error_bad_lines=False, low_memory=False)

#%% name the columns

df.columns = ["chr_s", "start_s", "end_s", "enh_id", "chr_e", "start_e",
"end_e", "seg_index", "core_remodeling", "core", "mrca", "core_mrca", "code",\
"syn_id"]
df.head()
#%%
if "mrca" in list(df):

    df = df.loc[df.mrca != "mrca"].copy() # rogue column name

    df["mrca"] = df["mrca"].astype(float).round(3) # round ages

#%% In[11]:

path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/download/tfbs/"

samples = glob.glob("%s*tfbs_density.bed" % path)
samples[0]

# concatenate the TFBS dataframe

shufEnhDict = {}
shufSynDict = {} # {sample_id: enh_tfbs_density df}

for sample in samples:

    df = pd.read_csv(sample, sep='\t', low_memory=False)

    if "mrca" in list(df):

        df = df.loc[df.mrca != "mrca"] # rogue column name

        df["mrca"] = df["mrca"].astype(float).round(3) # round ages

    dataName = "".join((sample.split("/")[-1]).split(".")[0])

    s_desc = (desc_df[1].loc[desc_df[0] == dataName]).to_string(index = False)# get the sample_description

    if "syn_tfbs" in sample:

        sdf = format_df(df, syn_gen_bkgd, dataName)
        sdf[['mrca','seg_index', 'core', 'core_remodeling', 'syn_len',
 'syn_tf_count', 'syn_tf_den', 'enh_len', 'syn_tf_u_count', 'syn_tf_u_den',
 'mrca_2',]] = sdf[[ 'mrca','seg_index', 'core', 'core_remodeling',
 'syn_len', 'syn_tf_count', 'syn_tf_den',
 'enh_len', 'syn_tf_u_count', 'syn_tf_u_den', 'mrca_2',]].astype(float)
        shufSynDict[dataName]= sdf

    elif "enh_tfbs" in sample: # deal with the enhancers and
        edf = format_df(df, syn_gen_bkgd, dataName)
        edf[[ 'core_remodeling', 'enh_len', 'enh_tf_count',
 'enh_tf_den', 'mrca', 'enh_tf_u_count', 'enh_tf_u_den',]]\
  = edf[[ 'core_remodeling', 'enh_len', 'enh_tf_count',
 'enh_tf_den', 'mrca', 'enh_tf_u_count', 'enh_tf_u_den',]].astype(float)

        shufEnhDict[dataName]= edf
#%%
samples
#%%

enhShufAll, synShufAll = pd.concat(shufEnhDict.values()), pd.concat(shufSynDict.values())
enhShuf = enhShufAll.drop(["datatype"], axis = 1).drop_duplicates()
synShuf = synShufAll.drop(["datatype"], axis = 1).drop_duplicates()

enhShuf["data"]= "SHUFFLE"
synShuf["data"]= "SHUFFLE"

#%%
dfEnhDict
#%%
enhAll, synAll,  = pd.concat(dfEnhDict.values()), pd.concat(dfSynDict.values())
enh = enhAll.drop(["datatype"], axis = 1).drop_duplicates()
syn = synAll.drop(["datatype"], axis = 1).drop_duplicates()

enh["data"]= "FANTOM"
syn["data"]= "FANTOM"

# # Simple v. Complex Enhancer lengths

# In[14]:


enh_len = enh[["enh_id", "core_remodeling","enh_len"]].drop_duplicates() # extract fantom enhancer lengths and ages
enh_len["datatype"] = "FANTOM" # mark datatype as fantom
shuf_len = enhShuf[["enh_id", "core_remodeling","enh_len"]].drop_duplicates()# extract shuffle lengths and ages
shuf_len["datatype"] = "shuffle" # mark datatype as shuffle
lens = pd.concat([enh_len, shuf_len]) # concatenate fantom + shuffle

#%%

colors = ["amber", "dusty purple"]
palette = sns.xkcd_palette(colors)
plt.rcParams.update({'font.size': 18})
fig, ax = plt.subplots(figsize = (10, 10))
sns.despine(bottom=True, left=True)
sns.boxplot(x = "core_remodeling", y = "enh_len", data = enh_len,\
             showfliers = False, palette = palette )
median = enh_len.groupby("core_remodeling")["enh_len"].median()

ax.set_title("Fantom Enhancer Lengths\nSimple v Complex")
ax.set_xlabel("Enhancer Architecture\n\n Median %s" % median)
ax.set_xticklabels(["Simple Enhancer", "Complex Enhancer"])
ax.set_ylabel("enhancer length")

plt.tight_layout()
plt.savefig("%s%sFANTOM_eRNA_enh_arch_lengths.pdf"% (RE,today))


# In[16]:


plt.rcParams.update({'font.size': 18})

fig, ax = plt.subplots(figsize = (10, 10))

sns.boxplot(x = "core_remodeling", y = "enh_len", data = lens,             hue = "datatype", showfliers = False, palette =epal, hue_order = enhOrder)
median = lens.groupby(["core_remodeling", "datatype"])["enh_len"].median()

ax.set_title("Fantom Enhancer Lengths\nSimple v Complex")
ax.set_xlabel("Enhancer Architecture\n\n Median %s" % median)
ax.set_xticklabels(["Simple Enhancer", "Complex Enhancer"])
ax.set_ylabel("enhancer length")
ax.legend(bbox_to_anchor = (1,1))
plt.tight_layout()
#plt.savefig("%s%sFANTOM_eRNA_enh_arch_lengths.pdf"% (ANALYSIS_PATH,today))


# # Simple v. Complex Enhancer lengths

# In[17]:


enh_len = enh[["enh_id", "core_remodeling","enh_len", "mrca_2", "taxon2"]].drop_duplicates() # extract fantom enhancer lengths and ages
enh_len["datatype"] = "FANTOM" # mark datatype as fantom
shuf_len = enhShuf[["enh_id", "core_remodeling","enh_len",  "mrca_2", "taxon2"]].drop_duplicates()# extract shuffle lengths and ages
shuf_len["datatype"] = "shuffle" # mark datatype as shuffle
lens = pd.concat([enh_len, shuf_len]) # concatenate fantom + shuffle


# In[18]:


breaks = syn.groupby(["enh_id", "core_remodeling", "data"])["mrca", "seg_index"].max().reset_index()
breaks = pd.merge(breaks, syn_gen_bkgd, how = "left", on = "mrca")
breaksShuf = synShuf.groupby(["enh_id", "core_remodeling", "data"])["mrca", "seg_index"].max().reset_index()
breaksShuf = pd.merge(breaksShuf, syn_gen_bkgd, how = "left", on = "mrca")
plotOrder = ["SHUFFLE", "FANTOM"]
breaks = pd.concat([breaks, breaksShuf])
fig, ax = plt.subplots(figsize = (8, 8))
sns.pointplot(x = "taxon2", y = "seg_index", data = breaks.loc[breaks.core_remodeling ==1].sort_values(by = "mrca_2"),
            hue = "data", hue_order = plotOrder, palette = epal, join = False, size = 12)
ax.set_xticklabels(ax.get_xticklabels(),rotation = 270)


# In[19]:



# What frequency of enhancers are remodeled in this dataset?

shuffle_remodeling = enhShuf[["enh_id", "core_remodeling"]].drop_duplicates()
shuffle_remodeling["core_remodeling"] = shuffle_remodeling["core_remodeling"].astype(int)
shuffle_remodeling["dataset"] = "Shuffle"

shuffle_remodeling["freq"] = np.where((shuffle_remodeling["core_remodeling"] ==1),
                                         (len(shuffle_remodeling.loc[shuffle_remodeling["core_remodeling"] ==1])/len(shuffle_remodeling)),\
                                         (len(shuffle_remodeling.loc[shuffle_remodeling["core_remodeling"] ==0])/len(shuffle_remodeling)))

remodeling = enh[["enh_id", "core_remodeling"]].drop_duplicates()
remodeling["core_remodeling"] = remodeling["core_remodeling"].astype(int)
remodeling["dataset"] = "FANTOM"

remodeling["freq"] = np.where((remodeling["core_remodeling"] ==1),
                                         (len(remodeling.loc[remodeling["core_remodeling"] ==1])/len(remodeling)),\
                                         (len(remodeling.loc[remodeling["core_remodeling"] ==0])/len(remodeling)))
remodeling_plot = pd.concat([remodeling, shuffle_remodeling])

remodeling_plot["arch"] = ""
remodeling_plot["arch"].loc[remodeling_plot["core_remodeling"] ==0] = "simple"
remodeling_plot["arch"].loc[remodeling_plot["core_remodeling"] ==1] = "complex"

remodeling_plot.head()


# In[20]:



arch_colors = [ "amber", "dusty purple", "windows blue","greyish"]
archpal = sns.xkcd_palette(arch_colors)

# What is the length of core v derived regions
arch_colors = [ "amber", "dusty purple", "windows blue","greyish"]
archpal = sns.xkcd_palette(arch_colors)
fig,(ax1, ax2) = plt.subplots(2,figsize = (8, 16))
p1 = sns.boxplot(y = "syn_len", x = "code", data = syn, order=archOrder,                 ax = ax1,  palette = archpal, showfliers=False)
median1 = syn.groupby("code")["syn_len"].median()
ax1.set_xlabel("Enhancer Architecture\n\nMedian Segment Length %s"% median1)
ax1.set_ylabel("syntenic length")
ax1.set_title("Enhancer Architecture Syntenic Lengths\n FANTOM")

# what is the frequency of core/derived across the len of the core
sns.boxplot(y = "seg_rep", x = "code", data=syn, ax = ax2,order=archOrder,  palette = archpal)

ax2.set_title("Enhancer Architecture Segment Occupancy\nFANTOM", fontsize = 16)
ax2.set_ylabel("% of enhancer")
median2 = syn.groupby("code")["seg_rep"].median()
ax2.set_xlabel("Enhancer Architecture\n\nMedian Syntenic Occupancy %s"% median2)

plt.savefig("%sFANTOM_enh_breaks_v_no_breaks_len_coverage.pdf"% (RE))
plt.show()


# In[21]:



# What is the length of core v derived regions
arch_colors = [ "amber", "dusty purple", "windows blue","greyish"]
archpal = sns.xkcd_palette(arch_colors)
fig,(ax1) = plt.subplots(1,figsize = (8, 8))
p1 = sns.boxplot(y = "syn_len", x = "code", data = syn, order=archOrder,                 ax = ax1,  palette = archpal, showfliers=False)
median1 = syn.groupby("code")["syn_len"].median()
ax1.set_xlabel("Enhancer Architecture\n\nMedian Segment Length %s"% median1)
ax1.set_xticklabels(["Simple Enhancer", "Complex Core", "Complex Derived"])
ax1.set_ylabel("syntenic length bp")
ax1.set_title("Enhancer Architecture Syntenic Lengths\n FANTOM")

plt.savefig("%sFANTOM_enh_breaks_v_no_breaks_len.pdf"% (RE))
plt.show()


# In[22]:


# What is the length of core v derived regions AGE STRATIFIED
fig,(ax1, ax2) = plt.subplots(2,figsize = (16, 16))
sns.pointplot(y = "syn_len", x = "taxon2", data=synShuf.sort_values(by="mrca_2"), ax = ax1            , hue = "code", hue_order =archOrder, palette = spal)
ax1.set_xticklabels(ax1.get_xticklabels(), rotation = 330, horizontalalignment = "left")
ax1.legend(bbox_to_anchor = (1,1))
ax2.set_title("what is the length of cores v derived segments?\n%speak" % (today))
ax2.set_xticklabels(ax2.get_xticklabels(), rotation = 330, horizontalalignment = "left")
# what is the frequency of core/derived across the len of the core
sns.pointplot(y = "seg_rep", x = "taxon2", data=syn.sort_values(by="mrca_2"),            ax = ax2, hue = "code",hue_order =archOrder,  palette = spal)

ax2.set_title("what percent of enhancer is cores v derived?", fontsize = 16)
ax2.set_ylabel("% of enhancer")
ax2.legend(bbox_to_anchor = (1,1))
plt.tight_layout()


# In[23]:


# What is the length of core v derived regions AGE STRATIFIED
arch_colors = [ "amber",  "dusty purple","windows blue","greyish"]
archpal = sns.xkcd_palette(arch_colors)
fig,(ax1, ax2) = plt.subplots(ncols = 2,figsize = (16, 8))
sns.pointplot(y = "syn_len", x = "taxon2", data=syn.sort_values(by="mrca_2"), ax = ax1            , hue = "code", hue_order =archOrder, palette = archpal, join = False, size = 15)
ax1.set_xticklabels(ax1.get_xticklabels(), rotation = 270, horizontalalignment = "left")
ax1.legend().remove()
ax2.set_title("what is the length of cores v derived segments?\n%speak" % (today))
ax2.set_xticklabels(ax2.get_xticklabels(), rotation = 270, horizontalalignment = "left")

sns.pointplot(y = "seg_rep", x = "taxon2", data=syn.sort_values(by="mrca_2"),            ax = ax2, hue = "code", hue_order =archOrder, palette = archpal, join = False, size = 15)

ax2.set_title("what percent of enhancer is cores v derived?", fontsize = 16)
ax2.set_ylabel("% of enhancer")
ax2.legend(bbox_to_anchor = (1,1))


# In[24]:


syntf = pd.concat([synShuf[["code", "data", "syn_tf_den"]], syn[["code", "data", "syn_tf_den"]]])


# In[25]:



fig, ax= plt.subplots(figsize = (6,6))
sns.despine(bottom=True, left=True)
sns.boxplot(x = "code", y = "syn_tf_den", data=syn,
            notch = True, order = archOrder, palette = archpal, showfliers=False)
medians = syn.groupby(["code"])['syn_tf_den'].median().reset_index()
means = syn.groupby(["code"])['syn_tf_den'].mean().reset_index()
print("medians", medians)
print("means", means)
plt.ylim(-0.05, 0.2)
ax.set_title("FANTOM eRNA")

# within remodeled enhancers is the difference between cores and derived regions
s_core, p_core = stats.mannwhitneyu(syn["syn_tf_den"].loc[(syn["code"]=="simple_core")],                          syn["syn_tf_den"].loc[(syn["code"]=="complex_core")])

s_complex, p_complex= stats.mannwhitneyu(syn["syn_tf_den"].loc[(syn["code"]=="complex_derived")],                          syn["syn_tf_den"].loc[(syn["code"]=="complex_core")])

kwstat, kw_p = stats.kruskal(syn["syn_tf_den"].loc[(syn["code"]=="simple_core")],                             syn["syn_tf_den"].loc[(syn["code"]=="complex_core")],                             syn["syn_tf_den"].loc[(syn["code"]=="complex_derived")])

ax.set_xticklabels(["Simple", "Complex\nCore", "Complex\nDerived"])
ax.set_xlabel("")
ax.set_ylim(0, 0.2)
ax.set_ylabel("TF ChIP-Seq sites/bp")
plt.tight_layout()
plt.savefig("%sFANTOM_enh_arch_tfbsden.pdf"% (RE))
print("MWU of complex core v. derived", s_complex, p_complex)
print("MWU of complex core v. simple_core",s_core, p_core)
plt.show()


# In[26]:


synShuf.head()


# In[27]:


plot = synShuf[["code", "data", "syn_tf_den", "syn_id"]].drop_duplicates()
plot.head()


# In[28]:


order


# In[ ]:



fig, ax= plt.subplots(figsize = (8,8))
sns.despine(bottom=True, left=True)
sns.barplot(x = "code", y = "syn_tf_den", data=plot, palette= spal, order =archOrder )
plt.savefig("%sFANTOM_SHUF_enh_arch_tfbsden.pdf"% (RE))


# In[ ]:


shuf_count = enhShuf.groupby(["enh_id"])["mrca"].max().reset_index()
shuf_count.head()


# In[ ]:


# count of shuffle enh_id
shuf_count = pd.merge(shuf_count, syn_gen_bkgd[["mrca", "taxon"]], how = "left", on = "mrca")
shuf_count = shuf_count.groupby(["mrca", "taxon"])["enh_id"].count().reset_index()
shuf_count["freq"] = shuf_count["enh_id"].divide( shuf_count["enh_id"].sum())
shuf_count["dataset"] = "Shuffle"
shuf_count.head()


# In[ ]:


# count of FANTOM enh_id
enh_count = enh.groupby(["enh_id"])["mrca"].max().reset_index()
enh_count = pd.merge(enh_count, syn_gen_bkgd[["mrca", "taxon"]], how = "left", on = "mrca")

enh_count = enh_count.groupby(["mrca", "taxon"])["enh_id"].count().reset_index()

enh_count["freq"] = enh_count["enh_id"].divide(enh_count["enh_id"].sum())
enh_count["dataset"] = "FANTOM"
enh_count.head()


# In[ ]:


enh_plot = pd.concat([enh_count, shuf_count])
enh_plot.head()


# In[ ]:


colors = ["greyish", "faded green"]
palette = sns.xkcd_palette(colors)

# plot count of syn_id
fig, (ax, ax2) = plt.subplots(2,figsize = (16,16))

sns.barplot(x= "taxon", y = "freq", data = enh_plot, hue = "dataset",
            palette = epal, hue_order = enhOrder, ax=ax)
ax.set_title("FANTOM Enhancer Distribution\nv. Shuffle background")
ax.set_ylabel("Frequency")
ax.set_xticklabels(ax.get_xticklabels(), rotation = 330, horizontalalignment = "left")
for p in ax.patches:
    ax.annotate('{:.0f}'.format(p.get_height()*100), (p.get_x()+0.005, p.get_height()+0.01), color = "grey",               fontsize = 12)

sns.pointplot(x= "taxon", y = "freq", data = enh_plot, hue = "dataset", palette = epal, hue_order = enhOrder,
              ax = ax2)
ax2.set_title("FANTOM Enhancer Distribution\nv. Shuffle background")
ax2.set_ylabel("Frequency")
ax2.set_xticklabels(ax2.get_xticklabels(), rotation = 330, horizontalalignment = "left")
mwu_stat, pval = stats.mannwhitneyu(syn["mrca"],synShuf["mrca"])
ax2.set_xlabel("taxon\n MWU = %s, pval = %s" % (mwu_stat, pval))
plt.tight_layout()
#plt.savefig("%sFANTOM_shuf_syn_mrca.pdf"% (RE))


# In[ ]:


#####FANTOM#####
# group by age & architecture
syn_counts= syn.groupby(["mrca", "code"])["syn_id"].count().reset_index()

# group by architecture
totals = syn_counts.groupby(["code"])["syn_id"].sum().reset_index()
totals.columns = ["code", "total_count"]

# calculate age v. architecture
concat = pd.merge(syn_counts, totals, how = "left", on = "code")
concat["frequency"] = concat["syn_id"].divide(concat["total_count"])
concat["dataset"] = "FANTOM"
concat = pd.merge(concat, syn_gen_bkgd, how = "left", on = "mrca")


#####FANTOM#####
# group by age & architecture
synShuf_counts= synShuf.groupby(["mrca", "code"])["syn_id"].count().reset_index()

# group by architecture
totals = synShuf_counts.groupby(["code"])["syn_id"].sum().reset_index()
totals.columns = ["code", "total_count"]

# calculate age v. architecture
concatShuf = pd.merge(synShuf_counts, totals, how = "left", on = "code")
concatShuf["frequency"] = concatShuf["syn_id"].divide(concatShuf["total_count"])
concatShuf["dataset"] = "Shuffle"
concatShuf = pd.merge(concatShuf, syn_gen_bkgd, how = "left", on = "mrca")

concat = pd.concat([concat, concatShuf])
concat["code2"] = concat.code + "-" + concat.dataset


# In[ ]:


concat = pd.concat([concat, concatShuf])
concat["code2"] = concat.code + "-" + concat.dataset


# In[ ]:


arch_colors = [ "windows blue", "amber", "dusty purple", "purplish", "lavender"]
archpal = sns.xkcd_palette(arch_colors)
fig, ax = plt.subplots(figsize=(16,8))
sns.barplot(x = "taxon", y = "syn_id", data = concat.sort_values(by=["mrca", "code"]),            hue = "code", alpha = 0.0,  palette = archpal, hue_order = archOrder)
#other plotating tools
ax.set_xticklabels(ax.get_xticklabels(), rotation=330, horizontalalignment = "left")
ax.set_xlabel("mrca")
ax.set_ylabel("enh syntenic block count")

ax.set_title("Age distribution of simple, complex enhancer segments")
ax2 = ax.twinx()
sns.barplot(x = "taxon", y = "frequency", data = concat.sort_values(by=["mrca", "code"]), hue = "code", ax = ax2,  palette = archpal)
ax2.set_ylabel("freq of complex|simple core dataset")
plt.tight_layout()
#plt.savefig("%sfantom_mrca_arch_dist.pdf"% (RE))


# In[ ]:


mwu_stat, pval = stats.mannwhitneyu(syn["mrca"].loc[(syn["code"].str.contains("simple_core"))],                   shuffle["mrca"].loc[(shuffle["code"].str.contains("simple_core"))])
print(mwu_stat, pval)
#ax2.set_xlabel("Simple core fantom v shuffle age dist\nMWU = %s, pval = %s" %(mwu_stat, pval))


# In[ ]:


pivot = pd.pivot_table(concat, values='frequency',                           index=['mrca', 'code', 'taxon'], columns=['dataset']).reset_index()
pivot["fold_change"] = pivot["FANTOM"].divide(pivot["Shuffle"]) -1
pivot.head()


# In[ ]:


arch_colors = [ "windows blue","greyish", "amber", "yellowish", "dusty purple",  "lavender"]
archpal = sns.xkcd_palette(arch_colors)
fig, (ax1, ax2, ax3) = plt.subplots(ncols = 3, figsize=(20,6))

#####DERIVED#####
sns.pointplot(x = "taxon", y = "frequency", data = concat.loc[concat["code2"].str.contains("derived")].sort_values(by=["mrca", "code2"]),              hue = "code2", ax = ax1,  palette = archpal)
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=330, horizontalalignment = "left", fontsize = 12)
mwu_stat, pval = stats.mannwhitneyu(syn["mrca"].loc[(syn["code"].str.contains("derived"))],                   synShuf["mrca"].loc[(synShuf["code"].str.contains("derived"))])
ax1.set_xlabel("mrca\nMWU = %s, pval = %s" %(mwu_stat, pval))
ax1.set_title("Complex Derived Age Distribution")
ax1.set_ylabel("Freq")
ax1.legend().remove()

#####COMPLEX_CORE#####
arch_colors = ["dusty purple",  "greyish"]
archpal = sns.xkcd_palette(arch_colors)
sns.pointplot(x = "taxon", y = "frequency", data = concat.loc[concat["code2"].str.contains("complex_core")].sort_values(by=["mrca", "code2"]),              hue = "code2", ax = ax2,  palette = archpal)
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=330, horizontalalignment = "left", fontsize = 12)
mwu_stat, pval = stats.mannwhitneyu(syn["mrca"].loc[(syn["code"].str.contains("complex_core"))],                   synShuf["mrca"].loc[(synShuf["code"].str.contains("complex_core"))])
ax2.set_xlabel("mrca\nMWU = %s, pval = %s" %(mwu_stat, pval))
ax2.set_title("Complex Core Age Distribution")
ax2.set_ylabel("Freq")
ax2.legend().remove()

#####SIMPLE_CORE#####
arch_colors = ["amber", "greyish"]
archpal = sns.xkcd_palette(arch_colors)
sns.pointplot(x = "taxon", y = "frequency", data = concat.loc[concat["code2"].str.contains("simple")].sort_values(by=["mrca", "code2"]),              hue = "code2", ax = ax3,  palette = archpal)
ax3.set_xticklabels(ax3.get_xticklabels(), rotation=300, horizontalalignment = "left", fontsize = 12)
mwu_stat, pval = stats.mannwhitneyu(syn["mrca"].loc[(syn["code"].str.contains("simple"))],                   synShuf["mrca"].loc[(synShuf["code"].str.contains("simple"))])
ax3.set_xlabel("mrca\nMWU = %s, pval = %s" %(mwu_stat, pval))
ax3.set_title("Simple Core Age Distribution")
ax3.set_ylabel("Freq")
ax3.legend().remove()

plt.tight_layout()
plt.savefig("%sMRCA_fantom_v_synShuf_core_derived_line.pdf"% (RE))


# In[ ]:


arch_colors = [ "windows blue","greyish", "amber", "yellowish", "dusty purple",  "lavender"]
archpal = sns.xkcd_palette(arch_colors)
fig, (ax1, ax2, ax3) = plt.subplots(ncols = 3, figsize=(20,6))

#####DERIVED#####
sns.barplot(x = "taxon", y = "fold_change",            data = pivot.loc[pivot["code"].str.contains("derived")].sort_values(by=["mrca", "code"]),              hue = "code", ax = ax1,  palette = archpal)

mwu_stat, pval = stats.mannwhitneyu(syn["mrca"].loc[(syn["code"].str.contains("derived"))],                   synShuf["mrca"].loc[(synShuf["code"].str.contains("derived"))])
ax1.set_xlabel("mrca\nMWU = %s, pval = %s" %(mwu_stat, pval))
ax1.set_title("Complex Derived Age Distribution")
ax1.set_ylabel("Fold Change")
ax1.legend().remove()

#####COMPLEX_CORE#####
arch_colors = ["dusty purple",  "greyish"]
archpal = sns.xkcd_palette(arch_colors)
sns.barplot(x = "taxon", y = "fold_change",            data = pivot.loc[pivot["code"].str.contains("complex_core")].sort_values(by=["mrca", "code"]),              hue = "code", ax = ax2,  palette = archpal)
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=270,  fontsize = 12)

mwu_stat, pval = stats.mannwhitneyu(syn["mrca"].loc[(syn["code"].str.contains("complex_core"))],                   synShuf["mrca"].loc[(synShuf["code"].str.contains("complex_core"))])
ax2.set_xlabel("mrca\nMWU = %s, pval = %s" %(mwu_stat, pval))
ax2.set_title("Complex Core Age Distribution")
ax2.set_ylabel("Fold Change")
ax2.legend().remove()

#####SIMPLE_CORE#####
arch_colors = ["amber", "greyish"]
archpal = sns.xkcd_palette(arch_colors)
sns.barplot(x = "taxon", y = "fold_change",            data = pivot.loc[pivot["code"].str.contains("simple")].sort_values(by=["mrca", "code"]),              hue = "code", ax = ax3,  palette = archpal)
ax3.set_xticklabels(ax3.get_xticklabels(), rotation=270,fontsize = 12)
ax1.set_xticklabels(ax3.get_xticklabels(), rotation=270, fontsize = 12)
mwu_stat, pval = stats.mannwhitneyu(syn["mrca"].loc[(syn["code"].str.contains("simple"))],                   synShuf["mrca"].loc[(synShuf["code"].str.contains("simple"))])
ax3.set_xlabel("mrca\nMWU = %s, pval = %s" %(mwu_stat, pval))
ax3.set_title("Simple Age Distribution")
ax3.set_ylabel("Fold Change")
ax3.legend().remove()

plt.tight_layout()
plt.savefig("%sFANTOM_v_synShuf_mrca_arch_fold_change.pdf"% (RE), bbox_inches = "tight")


# In[ ]:


arch_colors = ["amber", "dusty purple", "windows blue", ]
archpal = sns.xkcd_palette(arch_colors)
fig, ax1 = plt.subplots(figsize= (8,8))

sns.pointplot(x = "taxon", y = "fold_change",            data = pivot.sort_values(by=["mrca", "code"]),              hue = "code", ax = ax1,  palette = archpal, hue_order= archOrder)
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=270, horizontalalignment = "left", fontsize = 12)

mwu_stat, pval = stats.mannwhitneyu(syn["mrca"].loc[(syn["code"].str.contains("derived"))],                   synShuf["mrca"].loc[(synShuf["code"].str.contains("derived"))])
ax1.set_xlabel("mrca\nMWU = %s, pval = %s" %(mwu_stat, pval))
ax1.set_title("Age Distribution")
ax1.set_ylabel("Freq")
#ax1.legend().remove()
plt.savefig("%sFANTOM_v_shuffle_mrca_arch_fold_change_overlay.pdf"% (RE), bbox_inches = "tight")


# In[ ]:


concat.head()


# In[ ]:


#####FANTOM#####
# group by age & architecture
syn_counts= syn.groupby(["mrca_2", "code"])["syn_id"].count().reset_index()


# In[ ]:


# group by architecture
totals = syn_counts.groupby(["code"])["syn_id"].sum().reset_index()
totals.columns = ["code", "total_count"]
totals


# In[ ]:


# calculate age v. architecture
concat = pd.merge(syn_counts, totals, how = "left", on = "code")
concat["frequency"] = concat["syn_id"].divide(concat["total_count"])
concat["dataset"] = "FANTOM"
concat = pd.merge(concat, syn_gen_bkgd[["mrca_2", "taxon2"]], how = "left", on = "mrca_2").drop_duplicates()
concat


# In[ ]:


#####FANTOM#####
# group by age & architecture
synShuf_counts= synShuf.groupby(["mrca_2", "code"])["syn_id"].count().reset_index()

# group by architecture
totals = synShuf_counts.groupby(["code"])["syn_id"].sum().reset_index()
totals.columns = ["code", "total_count"]

# calculate age v. architecture
concatShuf = pd.merge(synShuf_counts, totals, how = "left", on = "code")
concatShuf["frequency"] = concatShuf["syn_id"].divide(concatShuf["total_count"])
concatShuf["dataset"] = "Shuffle"
concatShuf = pd.merge(concatShuf,  syn_gen_bkgd[["mrca_2", "taxon2"]], how = "left", on = "mrca_2").drop_duplicates()

concat = pd.concat([concat, concatShuf])
concat["code2"] = concat.code + "-" + concat.dataset


# In[ ]:


tx2 = concat.groupby(['mrca_2', 'code', 'taxon2', 'dataset' , "code2"])["frequency"].sum().reset_index()

tx2


# In[ ]:


arch_colors = [ "windows blue","greyish", "amber", "yellowish", "dusty purple",  "lavender"]
archpal = sns.xkcd_palette(arch_colors)
fig, (ax3, ax2, ax1) = plt.subplots(ncols = 3, figsize=(20,6))

#####DERIVED#####
sns.pointplot(x = "taxon2", y = "frequency",
              data = tx2.loc[tx2["code"].str.contains("derived")].sort_values(by=["mrca_2", "code"]),\
              hue = "code2", ax = ax1,  palette = archpal)
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=270, horizontalalignment = "left", fontsize = 12)
mwu_stat, pval = stats.mannwhitneyu(syn["mrca_2"].loc[(syn["code"].str.contains("derived"))],                   synShuf["mrca_2"].loc[(synShuf["code"].str.contains("derived"))])
ax1.set_xlabel("mrca\nMWU = %s, pval = %s" %(mwu_stat, pval))
ax1.set_title("Complex Derived Age Distribution")
ax1.set_ylabel("Freq")
ax1.legend().remove()

#####COMPLEX_CORE#####
arch_colors = ["dusty purple",  "greyish"]
archpal = sns.xkcd_palette(arch_colors)
sns.pointplot(x = "taxon2", y = "frequency",
              data = tx2.loc[tx2["code"].str.contains("complex_core")].sort_values(by=["mrca_2", "code"]),\
              hue = "code2", ax = ax2,  palette = archpal)
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=270, horizontalalignment = "left", fontsize = 12)
mwu_stat, pval = stats.mannwhitneyu(syn["mrca_2"].loc[(syn["code"].str.contains("complex_core"))],                   synShuf["mrca_2"].loc[(synShuf["code"].str.contains("complex_core"))])
ax2.set_xlabel("mrca\nMWU = %s, pval = %s" %(mwu_stat, pval))
ax2.set_title("Complex Core Age Distribution")
ax2.set_ylabel("Freq")
ax2.legend().remove()

#####SIMPLE_CORE#####
arch_colors = ["amber", "greyish"]
archpal = sns.xkcd_palette(arch_colors)
sns.pointplot(x = "taxon2", y = "frequency",
              data = tx2.loc[tx2["code"].str.contains("simple")].sort_values(by=["mrca_2", "code"]),\
              hue = "code2", ax = ax3,  palette = archpal)
ax3.set_xticklabels(ax3.get_xticklabels(), rotation=270, horizontalalignment = "left", fontsize = 12)
mwu_stat, pval = stats.mannwhitneyu(syn["mrca_2"].loc[(syn["code"].str.contains("simple"))],                   synShuf["mrca_2"].loc[(synShuf["code"].str.contains("simple"))])
ax3.set_xlabel("mrca\nMWU = %s, pval = %s" %(mwu_stat, pval))
ax3.set_title("Simple Core Age Distribution")
ax3.set_ylabel("Freq")
ax3.legend().remove()


plt.savefig("%sMRCA_fantom_v_synShuf_core_derived_grouped2.pdf"% (RE))


# In[ ]:


pivot2 = pd.pivot_table(tx2, values='frequency',                           index=['mrca_2', 'code', 'taxon2'], columns=['dataset']).reset_index()
pivot2["fold_change"] = pivot2["FANTOM"].divide(pivot2["Shuffle"]) -1
pivot2=pivot2.fillna(0)

arch_colors = [ "windows blue","greyish", "amber", "yellowish", "dusty purple",  "lavender"]
archpal = sns.xkcd_palette(arch_colors)
fig, (ax3, ax2, ax1) = plt.subplots(ncols = 3, figsize=(24,8))

#####DERIVED#####
sns.barplot(x = "taxon2", y = "fold_change",            data = pivot2.loc[pivot2["code"].str.contains("derived")].sort_values(by=["mrca_2", "code"]),              hue = "code", ax = ax1,  palette = archpal)
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=270, horizontalalignment = "left", fontsize = 12)
mwu_stat, pval = stats.mannwhitneyu(syn["mrca"].loc[(syn["code"].str.contains("derived"))],                   synShuf["mrca"].loc[(synShuf["code"].str.contains("derived"))])
ax1.set_xlabel("mrca\nMWU = %s, pval = %s" %(mwu_stat, pval))
ax1.set_title("Complex Derived Age Distribution")
ax1.set_ylabel("Fold Change")
ax1.legend().remove()

#####COMPLEX_CORE#####
arch_colors = ["dusty purple",  "greyish"]
archpal = sns.xkcd_palette(arch_colors)
sns.barplot(x = "taxon2", y = "fold_change",            data = pivot2.loc[pivot2["code"].str.contains("complex_core")].sort_values(by=["mrca_2", "code"]),              hue = "code", ax = ax2,  palette = archpal)
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=270, horizontalalignment = "left", fontsize = 12)
mwu_stat, pval = stats.mannwhitneyu(syn["mrca"].loc[(syn["code"].str.contains("complex_core"))],                   synShuf["mrca"].loc[(synShuf["code"].str.contains("complex_core"))])
ax2.set_xlabel("mrca\nMWU = %s, pval = %s" %(mwu_stat, pval))
ax2.set_title("Complex Core Age Distribution")
ax2.set_ylabel("Fold Change")
ax2.legend().remove()

#####SIMPLE_CORE#####
arch_colors = ["amber", "greyish"]
archpal = sns.xkcd_palette(arch_colors)
sns.barplot(x = "taxon2", y = "fold_change",            data = pivot2.loc[pivot2["code"].str.contains("simple")].sort_values(by=["mrca_2", "code"]),              hue = "code", ax = ax3,  palette = archpal)
ax3.set_xticklabels(ax3.get_xticklabels(), rotation=270, horizontalalignment = "left", fontsize = 12)
mwu_stat, pval = stats.mannwhitneyu(syn["mrca"].loc[(syn["code"].str.contains("simple"))],                   synShuf["mrca"].loc[(synShuf["code"].str.contains("simple"))])
ax3.set_xlabel("mrca\nMWU = %s, pval = %s" %(mwu_stat, pval))
ax3.set_title("Simple Age Distribution")
ax3.set_ylabel("Fold Change")
ax3.legend().remove()


plt.savefig("%sFANTOM_v_synShuf_mrca_arch_fold_change_tx2.pdf"% (RE), bbox_inches = "tight")


# In[ ]:


arch_colors = ["amber","dusty purple", "windows blue", ]
archpal = sns.xkcd_palette(arch_colors)
fig, ax1 = plt.subplots(figsize= (8,8))

sns.pointplot(x = "taxon2", y = "fold_change",            data = pivot2.sort_values(by=["mrca_2", "code"]),              hue = "code", ax = ax1,  palette = archpal, hue_order = archOrder)
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=270, horizontalalignment = "left", fontsize = 12)

mwu_stat, pval = stats.kruskal(syn["mrca_2"].loc[(syn["code"].str.contains("derived"))],                   syn["mrca_2"].loc[(syn["code"].str.contains("complex_core"))],
                              syn["mrca_2"].loc[(syn["code"].str.contains("simple_core"))])
ax1.set_xlabel("mrca\nKruskal = %s, pval = %s" %(mwu_stat, pval))
ax1.set_title("Enhancer Architecture Age Distribution")
ax1.set_ylabel("Fold-Change")
#ax1.legend().remove()
plt.savefig("%sFANTOM_v_shuffle_mrca_arch_fold_change_tx2_overlay.pdf"% (RE), bbox_inches = "tight")


# In[ ]:


# abs enh tf density v. break density per segment
arch_colors = ["amber", "dusty purple", "windows blue","greyish"]
archpal = sns.xkcd_palette(arch_colors)
fig, (ax, ax2) = plt.subplots(2, figsize = (15, 20))

sns.barplot(y = "syn_tf_den", x = "taxon", data = syn.sort_values(by = "mrca"), ax = ax,            hue = "code", palette = archpal, hue_order = archOrder)
ax.set_ylabel("tfbs den/syntenic block")
ax.set_xticklabels(ax.get_xticklabels(), rotation=270, horizontalalignment = "left", fontsize =10)
ax.set_title("Derived v. Core TFBS density at each age")
sns.pointplot(y = "syn_tf_den", x = "taxon", data = syn.sort_values(by = "mrca"), ax= ax2,              hue = "code", palette = archpal, hue_order = archOrder)
ax2.set_ylabel("tfbs den/syntenic block")
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=270, horizontalalignment = "left", fontsize =10)
ax2.set_title("Derived v. Core TFBS density at each age")


plt.savefig("%stfbs_density_per_mrca_plus_arch_bar.pdf"% (RE))
plt.show()


# In[ ]:


syn.head()


# In[ ]:


# abs enh tf density v. break density per segment
arch_colors = ["amber", "dusty purple", "windows blue","greyish"]
archpal = sns.xkcd_palette(arch_colors)
fig, ax = plt.subplots( figsize = (8, 8))

sns.barplot(y = "syn_tf_den", x = "taxon2", data = syn.sort_values(by = "mrca_2"), ax = ax,            hue = "code", palette = archpal, hue_order = archOrder)
ax.set_ylabel("tfbs den/syntenic block")
ax.set_xticklabels(ax.get_xticklabels(), rotation=270,  fontsize =10)
ax.set_title("Derived v. Core TFBS density at each age")


plt.savefig("%stfbs_density_per_mrca_plus_arch_bar.pdf"% (RE))
