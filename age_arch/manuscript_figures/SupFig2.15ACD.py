#!/usr/bin/env python
# coding: utf-8

# In[32]:


import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os, sys

from scipy import stats
import seaborn as sns
import statsmodels
import statsmodels.api as sm

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

faded_green = "#7bb274"
slate_grey = "#59656d"
amber = "#feb308"
dusty_lavender = "#ac86a8"
dusty_purple = "#825f87"
windows_blue = "#3778bf"
greyish = "#a8a495"
greyish_blue = "#5e819d"
greyish_purple = "#887191"
dull_yellow = "#eedc5b"

# sns graphing preferences
sns.set(color_codes=True)
sns.set(font_scale=1.5)
sns.set_style("white")
sns.despine(bottom=True, left=True)

import datetime
last_run = datetime.datetime.now()

today = (datetime.date.today())
print("last run", datetime.

datetime.now())

RE ="/dors/capra_lab/projects/enhancer_ages/roadmap_encode/results/for_publication/non-genic/trimmed/"


# In[33]:


es_colors = [ "slate grey","greyish"]
es_pal = sns.xkcd_palette(es_colors)
sns.palplot(es_pal)


# In[34]:


path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19/download/h3k27ac_plus_h3k4me3_minus_peaks/all_roadmap_enh/trimmed/trimmed-310-all_roadmap_enh/"


# # load dataframes

# In[36]:
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


final_merge = pd.read_csv(enh, sep = '\t', header = None) # open enhancer file

final_merge.columns = columns # rename columns

print("enhancer relative simple cut off", final_merge.seg_index.median())

final_merge.loc[final_merge.seg_index == 1, "core_remodeling"] = 0

final_merge.loc[final_merge.seg_index == 1, "arch"] = "simple"

final_merge.mrca_2 = final_merge.mrca_2.round(3) # round values

final_merge["enh_len"] = final_merge.end_enh - final_merge.start_enh
final_merge["datatype"] = "roadmap_310"


#%% count of simple and complex enhancers per shuffled dataset.
shuf_arch = shuffle[["enh_id", "core_remodeling", "shuf_id"]].drop_duplicates()
shuf_arch_freq = shuf_arch.groupby(["core_remodeling", "shuf_id"])["enh_id"].count().reset_index()


totals = shuf_arch.groupby(["shuf_id"])["enh_id"].count().reset_index()
totals.columns = ["shuf_id", "totals"]
shuf_arch_freq = pd.merge(shuf_arch_freq, totals, how = "left")
shuf_arch_freq["freq"] = shuf_arch_freq["enh_id"].divide(shuf_arch_freq.totals)
shuf_arch_freq["dataset"] = "SHUFFLE"

# count of simple and complex enhancers in roadmap_310
arch = final_merge[["enh_id", "core_remodeling"]].drop_duplicates()
arch_freq = arch.groupby("core_remodeling")["enh_id"].count().reset_index()
totals = len(arch)
arch_freq["freq"] = arch_freq["enh_id"].divide(totals)
arch_freq["dataset"] = "roadmap_310"

archs = pd.concat([shuf_arch_freq, arch_freq]) # combine datasets for plotting


# In[37]:


final_merge.seg_index.min()


# In[ ]:


shuf_colors = [ "amber", "greyish",]
shuf_pal = sns.xkcd_palette(shuf_colors)
hue_order = ["roadmap_310", "SHUFFLE"]

fig, ax = plt.subplots(figsize = (8, 8))
sns.barplot(x = "core_remodeling", y="freq", data = archs.loc[archs.core_remodeling ==0],
            hue = "dataset", hue_order = hue_order, palette = shuf_pal)

ax.set_ylabel("Frequency of Dataset")
for p in ax.patches:
    x=p.get_bbox().get_points()[:,0]
    y=p.get_bbox().get_points()[1,1]
    ax.annotate('{:.0f}%'.format(100.*y), (x.mean(), y),
            ha='left', va='bottom', color = "k", alpha = 0.4, fontsize = 20) # set the alignment of the text
#ax.legend(bbox_to_anchor=(1.5, 1.0))

ax.set_xticklabels("")
ax.set_xlabel("")
ax.set_title("")
ax.get_legend().remove()
sns.set(font_scale=2)
sns.set_style("white")
#plt.savefig("%sroadmap_310_enh_shuf_freq.pdf" %RE, bbox_inches = "tight")


# In[7]:


fig, ax = plt.subplots(figsize = (8, 8))
sns.barplot(x = "core_remodeling", y="freq", data = arch_freq, palette = palette)
ax.set_xticklabels(["Simple", "Complex\nEnhancer"])
ax.set_xlabel("")
ax.set_ylabel("Frequency of Dataset")
ax.set_title("roadmap_310 Enhancer Architectures")
for p in ax.patches:
    x=p.get_bbox().get_points()[:,0]
    y=p.get_bbox().get_points()[1,1]
    ax.annotate('{:.0f}%'.format(100.*y), (x.mean(), y),
            ha='left', va='bottom', color = "k", alpha = 0.4, fontsize = 20) # set the alignment of the text
ax.set_xticklabels("")
ax.set_xlabel("")
ax.set_title("")
#ax.get_legend().remove()
plt.savefig("%sroadmap_310_enh_arch_freq.pdf" %RE, bbox_inches = "tight")


# # breaks

# In[38]:


breaks = final_merge[["enh_id", "core_remodeling", "seg_index"]].drop_duplicates()
breaks.head()
#%%
sbreaks = shuffle[["enh_id", "core_remodeling", "shuf_id", "seg_index"]]
sbreaks.head()


# In[39]:


sbreaks.min()


# In[40]:


# get shuffled breaks distribution

shuf_break_freq = sbreaks.groupby(["seg_index", "shuf_id"])["enh_id"].count().reset_index()

shuf_totals = shuf_break_freq.groupby("shuf_id")["enh_id"].sum().reset_index()
shuf_totals.columns = ["shuf_id", "shuf_id_totals"]

shuf_break_freq = pd.merge(shuf_break_freq, shuf_totals, how = "left", on = "shuf_id")

shuf_break_freq["freq"] = shuf_break_freq["enh_id"].divide(shuf_break_freq.shuf_id_totals)
shuf_break_freq["dataset"] = "Shuffle"


# In[41]:


sbreaks.seg_index.mean()


# In[42]:


breaks.seg_index.mean()


#%% # plot cumulative dist

enh_seg = final_merge[["enh_id", "seg_index"]].drop_duplicates()
enh_seg["cdf"] = enh_seg.seg_index.rank(pct =True)
enh_seg["dataset"] = "roadmap_310"

shuf_seg = shuffle[["enh_id", "seg_index"]]
shuf_seg["cdf"] = shuf_seg.seg_index.rank(pct =True)
shuf_seg["dataset"] = "shuffle"

plot_cdf = pd.concat([enh_seg, shuf_seg])
plot_cdf.head()



#%%

fig, ax = plt.subplots(figsize = (8,8))
x = "seg_index"
y = "cdf"


sns.lineplot(x, y, data = plot_cdf, ax = ax, hue = "dataset", palette = shuf_pal)
ax.set(xticks = (np.arange(0, plot_cdf.seg_index.max(), step = 5)),
    xlabel = "number of segments",
    ylabel = "cumulative distribution",
    xlim = (0,10))

#plt.savefig("%sroadmap_310_CDF.pdf" %RE, bbox_inches = 'tight')



# In[48]:



OR_dict = {}

for seg_index in plot_cdf.seg_index.unique():

    enh_totals = len(plot_cdf.loc[plot_cdf.dataset == "roadmap_310"])
    shuf_totals = len(plot_cdf.loc[plot_cdf.dataset == "shuffle"])

    a = len(plot_cdf.loc[(plot_cdf.seg_index ==seg_index)
    & (plot_cdf.dataset == "roadmap_310")]) # num simple enhancers
    b = enh_totals - a # num complex enhancers
    c = len(plot_cdf.loc[(plot_cdf.seg_index ==seg_index)
    & (plot_cdf.dataset == "shuffle")]) # num simple shuffle
    d = shuf_totals - c # num complex shuffle

    obs = [[a,b], [c,d]]
    OR, P = stats.fisher_exact(obs)
    table = sm.stats.Table2x2(obs) # get confidence interval
    odds_ci = table.oddsratio_confint()
    newdf = pd.DataFrame({"seg_index":[seg_index], "a":[a], "b":[b], "c":[c], "d":[d],
                         "OR":[OR], "P":[P], "ci_lower" :[odds_ci[0]],
                        "ci_upper" :[odds_ci[1]]})
    OR_dict[seg_index] = newdf
    print(seg_index, obs, OR, P)


#%%


ORdf = pd.concat(OR_dict.values())
ORdf["yerr"] = ORdf.ci_upper-ORdf.ci_lower
ORdf['log'] = np.log2(ORdf.OR)
ORdf.head()


# In[53]:


fig, ax = plt.subplots(figsize =(8,8))
sns.set("poster")

firstfive = ORdf.loc[ORdf.seg_index <6]
sns.barplot(x = "seg_index", y = "log", data = firstfive.loc[firstfive.seg_index <6],
            linewidth=2.5, facecolor=(1, 1, 1, 0),
             edgecolor=".2",  yerr=(firstfive["ci_upper"] - firstfive["ci_lower"]))
ax.set_ylabel("Fold Change v. Bkgd\n(log2-scaled)")
ax.set_xlabel("Number of Age Segments")
plt.axhline(0, color = "grey", linewidth = 2.5)
from matplotlib.ticker import MultipleLocator
import matplotlib.ticker as ticker
ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(round(2**x, 2)))
ax.yaxis.set_major_formatter(ticks)
ax.yaxis.set_major_locator(MultipleLocator(0.05))
#ax.set_ylim(-1.2,0.5)

plt.savefig("%sSupFig2.15A-roadmap_310_age_seg_fold_change.pdf" %RE, bbox_inches = "tight")


# In[54]:




shuf_colors = [ "amber", "greyish","dusty purple", "windows blue","greyish"]
shuf_pal = sns.xkcd_palette(shuf_colors)

hue_order = ["roadmap_310", "Shuffle"]
fig, ax = plt.subplots(figsize = (8,8))
sns.set("poster")
sns.barplot(x = "seg_index", y = "cdf", data = plot_cdf, hue = "dataset", hue_order = hue_order, palette = es_pal)
ax.set_xlabel("Number of Age Segments per Enhancer")
ax.set_ylabel("Cumulative Frequency of Dataset")
ax.set_title("roadmap_310 Enhancer Age Segment Count\n Cumulative distribution")
plt.legend()
ax.set_xlim(-1, 5.5)
#ax.set_xticklabels("")

ax.set_title("")
ax.get_legend().remove()
sns.set("poster")
plt.savefig("%sSupFig2.15B-roadmap_310_age_seg_cdf.pdf" %RE, bbox_inches = "tight")


# In[55]:

#%%
print("hello")
#%%

enh_lens = final_merge[["enh_id", "core_remodeling","enh_len", "arch", "mrca_2", "taxon2", "mya2"]].drop_duplicates()
print(enh_lens.shape)

enh_lens["datatype"]="roadmap_310"
enh_lens.head()


# In[56]:

shuf_len = shuffle[["enh_id", "core_remodeling", "enh_len", "arch", "mrca_2", "taxon2", "mya2"]].drop_duplicates()
print(shuf_len.shape)

shuf_len["datatype"]="SHUFFLE"

print(shuf_len.shape)

sample_shuf = shuf_len.sample(frac = 0.1)

lens = pd.concat([enh_lens, sample_shuf])
lens.head()


# In[57]:


def get_age_arch_freq(df, datatype, dataset):

    a = df.loc[df[datatype] == dataset].groupby(["mrca_2", "core_remodeling"])["enh_len"].count().reset_index()

    b = a.pivot(index = "mrca_2", columns = "core_remodeling", values = "enh_len").reset_index()

    b["totals"] = b[0] + b[1].fillna(0)
    b["simple_freq"] = b[0]/ b["totals"]
    b["complex_freq"] = b[1]/ b["totals"]
    b["total_freq"]= 1
    b["dataset"]= dataset

    b = b.fillna(0)

    return b


# In[58]:


lens_dict = {}

# enhancer frequencies

lens_dict["roadmap_310"] = get_age_arch_freq(lens, "datatype", "roadmap_310")
lens_dict["roadmap_310"].dataset.unique()


#%%
lens.head()
#%%

for i in lens.datatype.unique():
    print(i)

    d = get_age_arch_freq(lens, "datatype", str(i))
    lens_dict[i] = d
#%%
for_stats = pd.concat(lens_dict.values())
roadmap_310_stats = for_stats.loc[for_stats.dataset == "roadmap_310"]
for_stats = for_stats.loc[for_stats.dataset != "roadmap_310"]

m, p = stats.mannwhitneyu(roadmap_310_stats.simple_freq, for_stats.simple_freq)
print(m, p)

m, p = stats.kruskal(roadmap_310_stats.simple_freq.loc[roadmap_310_stats["mrca_2"] == 0.000],
                     roadmap_310_stats.simple_freq.loc[roadmap_310_stats["mrca_2"] == 0.126],
                     roadmap_310_stats.simple_freq.loc[roadmap_310_stats["mrca_2"] == 0.131],
                     roadmap_310_stats.simple_freq.loc[roadmap_310_stats["mrca_2"] == 0.152],
                     roadmap_310_stats.simple_freq.loc[roadmap_310_stats["mrca_2"] == 0.175],
                     roadmap_310_stats.simple_freq.loc[roadmap_310_stats["mrca_2"] == 0.308],
                     roadmap_310_stats.simple_freq.loc[roadmap_310_stats["mrca_2"] == 0.38],
                     roadmap_310_stats.simple_freq.loc[roadmap_310_stats["mrca_2"] == 0.49],
                     roadmap_310_stats.simple_freq.loc[roadmap_310_stats["mrca_2"] == 0.656],
                     roadmap_310_stats.simple_freq.loc[roadmap_310_stats["mrca_2"] == 0.957],
                         )
print(m, p)

roadmap_310_stats.simple_freq.mean()

for_stats.simple_freq.mean()



#%%
hue_order = ["roadmap_310", "SHUFFLE"]
fig, ax = plt.subplots(figsize = (8, 8))
sns.barplot(x = "taxon2", y = "enh_len", hue = "datatype",
            data = lens.sort_values(by = 'mrca_2'), palette = es_pal, hue_order = hue_order)

ax.set_xlabel("Architecture")
ax.set_ylabel("Enhancer Length")
ax.legend(bbox_to_anchor=(1.0, 1.0))
ax.set_xticklabels(ax.get_xticklabels(), rotation = 90, horizontalalignment = "left")


# Get medians

lens_data = lens.groupby(["core_remodeling", "datatype"])["enh_len"].median().reset_index()
lens_data_mean = lens.groupby(["core_remodeling", "datatype"])["enh_len"].mean().reset_index()
print("medians", lens_data,)
print("means", lens_data_mean,)
ax.set_ylim(225,375)
ax.set_xlabel("")
ax.set_title("")
ax.get_legend().remove()
#plt.savefig("%sfig1d-roadmap_310_enh_mrca_lens.pdf" %RE, bbox_inches = "tight")
#%%


# In[71]:


x = len(enh_lens)
enh_len_freq = enh_lens.groupby(["mrca_2", "core_remodeling"])["enh_id"].count().divide(x).round(3).reset_index()
enh_len_freq.columns = ["mrca_2", "core_remodeling", "enh_freq"]

y = len(shuf_len)
shuf_len_freq = shuf_len.groupby(["mrca_2", "core_remodeling"])["enh_id"].count().divide(y).round(3).reset_index()
shuf_len_freq.columns = ["mrca_2", "core_remodeling", "shuf_freq"]

len_freq = pd.merge(enh_len_freq, shuf_len_freq, how = 'left' )
len_freq["dif"] = len_freq.enh_freq / len_freq.shuf_freq
len_freq

shuf_len = shuf_len.loc[~(shuf_len.arch.str.contains("complex")&(shuf_len.mrca_2 == 0.0))]
enh_lens = pd.merge(enh_lens, len_freq, how = "left")
shuf_len = pd.merge(shuf_len, len_freq, how = "left")
shuf_len.head()


# In[72]:


enh_lens.drop_duplicates().shape


# In[73]:



e_colors = [ "amber", "faded green"]
e_pal = sns.xkcd_palette(e_colors)
s_colors = [ "greyish", "slate grey"]
s_pal = sns.xkcd_palette(s_colors)

hue_order = ["roadmap_310", "Shuffle"]
fig,(ax1) = plt.subplots(figsize = (8, 8))
order =["Simple", "Complexenh"]
sns.barplot(y = "enh_len", x = "taxon2",
data = enh_lens.sort_values(by = "mrca_2"), ax = ax1,
            hue = "arch",  palette = e_pal, estimator = np.median)#showfliers=False)


ax1.set(ylabel="Enhancer Length (bp)",
ylim= (190,400), xlabel = "", title = "")
ax1.set_xticklabels(ax1.get_xticklabels(), rotation = 90, horizontalalignment = "left")
ax1.get_legend().remove()
plt.savefig("%sSupFig2.15C-roadmap_310_ENH_MRCA_x_LEN_ENH.pdf" % RE, bbox_inches = "tight")


# In[74]:


len_concat = pd.concat([enh_lens, shuf_len])

len_concat.head()

len_concat["datatype2"] = len_concat.datatype + "-"+ len_concat.arch
len_concat["datatype2"].unique()


e_colors = [ "amber","greyish", "faded green", "slate grey"]
e_pal = sns.xkcd_palette(e_colors)
hue_order = ['roadmap_310-simple','SHUFFLE-simple','roadmap_310-complexenh','SHUFFLE-complexenh']


fig,(ax1) = plt.subplots(figsize = (8, 8))
order =["Simple", "Complexenh"]
sns.barplot(y = "enh_len", x = "taxon2",
            data = len_concat.sort_values(by = "mrca_2"), ax = ax1,
            hue = "datatype2",  palette = e_pal,
            hue_order = hue_order,
            estimator = np.median)#showfliers=False)


ax1.set_ylabel("Enhancer Length (bp)")
ax1.set_xticklabels(ax1.get_xticklabels(), rotation = 90, horizontalalignment = "left")

ax1.set_xlabel("")
ax1.set_title("")
ax1.legend().remove()
plt.savefig("%sSupFig2.15.7-roadmap_310_ENH_MRCA_x_LEN_ENH.pdf" % RE, bbox_inches = "tight")

#%%

print(len_concat.core_remodeling.unique())
print(len_concat["datatype2"].unique())

#%%


enh_lens.groupby("core_remodeling")["mrca_2"].mean()


# In[78]:


shuf_len.groupby("core_remodeling")["mrca_2"].mean()
#%%
list(shuffle)
#%%

# # enhancer age segments

core_breaks = final_merge[["enh_id", "core_remodeling", "mrca_2", "seg_index", "taxon2"]]


core_breaks["datatype"] = "roadmap_310"
core_breaks=core_breaks.drop_duplicates()

smrca = shuffle[["enh_id", "core_remodeling", "mrca_2", "seg_index", "taxon2"]]

# DATA CLEAN UP - remove the enh_id errors in the shuffle db where the oldest complexenh age is homosapiens
smrca_list = list(smrca.loc[(smrca.core_remodeling == 1)
                 & (smrca.mrca_2 == 0.00)]["enh_id"])

smrca = smrca[~smrca.enh_id.isin(smrca_list)]

sseg = shuffle[["enh_id", "core_remodeling", "datatype", "seg_index"]]

shuf_core_breaks = pd.merge(smrca, sseg, how ="left", on = ("enh_id", "core_remodeling", "seg_index"))



shuf_core_breaks=shuf_core_breaks.drop_duplicates()
#%%
shuf_core_breaks.head()
#%%

plot_breaks=pd.concat([core_breaks,shuf_core_breaks])
plot_breaks.datatype.unique()
#%%
plot_breaks.groupby(["core_remodeling", "taxon2", "datatype"])["seg_index"].mean()
#%%
f, ax = plt.subplots(figsize = (8,8))

hue_order = ["roadmap_310", "SHUFFLE"]
sns.barplot(x = "taxon2", y = "seg_index",
            data =plot_breaks.loc[plot_breaks.core_remodeling ==1].sort_values(by = "mrca_2"),
            palette = cs_pal, hue = "datatype", hue_order = hue_order)

ax.set_xticklabels(ax.get_xticklabels(), rotation = 90, horizontalalignment = "left")
ax.set_ylabel("# age segments")
ax.set_title("roadmap_310 age segments v MRCA")
ax.set_xlabel("")
ax.set_title("")
ax.get_legend().remove()
ax.ylim = (0,3.2)
plt.savefig("%sFigS2.4B-roadmap_310_age_segments_v_MRCA_bar.pdf" %RE, bbox_inches = "tight")


# In[83]:


f, ax = plt.subplots(figsize = (8,8))

hue_order = ["roadmap_310", "SHUFFLE"]
sns.barplot(x = "datatype", y = "seg_index",
            data =plot_breaks[plot_breaks.core_remodeling ==1],
            palette = cs_pal)

ax.set_xticklabels(ax.get_xticklabels(), rotation = 90, horizontalalignment = "left")
ax.set_ylabel("# age segments")
ax.set_title("roadmap_310 v Shuffle age segments")
ax.set_xlabel("")
ax.set_title("")
ax.legend().remove()
ax.ylim = (0,3.2)
plt.savefig("%sFigS2.4A-roadmap_310_age_segments_bar.pdf" %RE, bbox_inches = "tight")


# In[84]:


from scipy.stats import mstats
kw_list = []
roadmap_310_list = []
shuffle_list = []
for i in plot_breaks.mrca_2.unique():
    for d in plot_breaks.datatype.unique():
        mrca_list = plot_breaks.loc[(plot_breaks.mrca_2 == i)
                                    & (plot_breaks.datatype == d)
                                    & (plot_breaks.core_remodeling ==1) , "seg_index"].to_list() # collect the # segments per age
        kw_list.append(mrca_list)
        if d == "roadmap_310":
            roadmap_310_list.append(mrca_list) # collect the # segments per age roadmap_310
        elif d == "SHUFFLE":
            shuffle_list.append(mrca_list)# collect the # segments per age shuffle

args=[l for l in kw_list]
argsF=[l for l in roadmap_310_list]
argsS=[l for l in shuffle_list]
stats.mstats.kruskalwallis(*args) # KW on age segment differences for all age and datatype complex enhancers

stats.mstats.kruskalwallis(*argsF) # KW on age segment differences among roadmap_310 complex enhancers

stats.mstats.kruskalwallis(*argsS) # KW on age segment differences among Shuffle complex enhancers


# In[85]:


# MWU on age segments in roadmap_310 v. Shuffle
stats.mannwhitneyu(plot_breaks.loc[(plot_breaks.datatype == "roadmap_310") & (plot_breaks.core_remodeling ==1), "seg_index"],
                  plot_breaks.loc[(plot_breaks.datatype == "SHUFFLE") & (plot_breaks.core_remodeling ==1), "seg_index"],)


# In[86]:


plot_breaks.groupby(["datatype", "core_remodeling"])["seg_index"].mean()


# In[87]:


plot_breaks.groupby(["datatype"])["seg_index"].mean()


# # Enhancer age fold-change

# In[88]:


enh_ages = final_merge.groupby(['enh_id'])["mrca_2"].max().reset_index()
enh_ages["data_type"] = "roadmap_310"
shuf_ages = shuffle.groupby(['enh_id', 'shuf_id',])["mrca_2"].max().reset_index()
shuf_ages["data_type"] = "shuffle"

shuf_ages.head()


# In[89]:


ANNOTATIONS = 0

shuffle_dist = shuffle[["enh_id", "datatype", "mrca_2", "taxon2"]]


shuffle_dist2 = shuffle_dist.groupby(["mrca_2", "datatype", "taxon2"])["enh_id"].count().reset_index()
shuffle_dist2.columns = ["mrca_2", "datatype", "taxon2", "mrca_count"]

totals = shuffle_dist2.groupby(["datatype"])["mrca_count"].sum().reset_index()
totals.columns = ["datatype", "id_totals"]

shuffle_dist2 = pd.merge(shuffle_dist2, totals, how = "left")
shuffle_dist2['freq'] = shuffle_dist2.mrca_count.divide(shuffle_dist2.id_totals)
shuffle_dist2["dataset"] = "shuffle"


# In[90]:


plot_dist = final_merge.groupby("enh_id")["mrca_2"].max().reset_index() # calculate age distribution of roadmap_310 enhancers by enh_id

plot_dist["mrca_2"] = plot_dist["mrca_2"].round(3) # round ages

totals = plot_dist.groupby("mrca_2")["enh_id"].count().reset_index()
totals.columns = ["mrca_2", "mrca_count"]
totals['freq'] = totals.mrca_count.divide(totals.mrca_count.sum())

plot_dist = pd.merge(plot_dist, totals, how = "left", on = "mrca_2")

plot_dist["dataset"] = "roadmap_310"
dist = pd.concat([shuffle_dist2, plot_dist]) # concat shuffle and roadmap_310 age distributions
dist.head()


# In[91]:





# calculate MWU shuffle v. roadmap_310 age distributions (from each enhancer)

mwustat, mwupval = stats.mannwhitneyu(shuffle_dist2["mrca_2"].loc[shuffle_dist2["dataset"] =="shuffle"],
                                     plot_dist["mrca_2"].loc[plot_dist["dataset"] =="roadmap_310"] )
dist = dist[["taxon2", "mrca_2", "dataset", "freq"]].drop_duplicates()

dist.sort_values(by = "mrca_2").head()


# In[92]:


# Fold change of roadmap_310/shuffle age frequency
fold_change_dict = {}
for i in dist.dataset.unique():

    fold_change = dist.pivot(index="mrca_2", columns='dataset')['freq'].reset_index()

    fold_change["fold_change"] = np.log2(fold_change["roadmap_310"].divide(fold_change["shuffle"]))
    fold_change = pd.merge(fold_change, dist[["mrca_2", "taxon2"]], how = "left", on = "mrca_2").drop_duplicates()
    fold_change_dict[i] = fold_change

fc = pd.concat(fold_change_dict.values())


dist

#%%
colors = ["greyish", "slate grey",]
palette = sns.xkcd_palette(colors)
sns.set("poster")
#
mwustat, mwupval = stats.mannwhitneyu(shuffle_dist2["mrca_2"].loc[shuffle_dist2["dataset"] =="shuffle"],
                                     plot_dist["mrca_2"].loc[plot_dist["dataset"] =="roadmap_310"] )
fig, ( ax2, ax1) = plt.subplots(ncols = 2, figsize=(16,8))

# line plot

sns.barplot(x = "mrca_2", y= "freq",  data = dist.sort_values(by="mrca_2"),
            hue = "dataset", palette = palette, ax = ax2)

ax2.set_xticklabels(ax2.get_xticklabels(), rotation = 90, horizontalalignment="left")
ax2.set_ylabel("Frequency")
ax2.set_xlabel("Taxon \n MWU = %s, pval = %s" % (mwustat, mwupval))
ax2.set_ylim(0, 0.6)
ax2.set_title("roadmap_310 Enhancer Age Distribution")
ax2.set_xlabel("")
ax2.set_title("")
ax2.legend(frameon = False)

# fold change plot
colors = ["slate grey"]
palette_fc = sns.xkcd_palette(colors)
sns.barplot(x = "mrca_2", y= "fold_change",  data = fc.sort_values(by="mrca_2"),
          palette = palette_fc, ax = ax1)
if ANNOTATIONS == 1:
    for p in ax1.patches:
        x = p.get_x().round(2)
        y = (p.get_height().round(2))
        ax1.annotate('{:.1f}'.format(p.get_height()), (x, y), color = "black", fontsize = 20)
ax1.set_xticklabels(ax1.get_xticklabels(), rotation = 90, horizontalalignment="left")

ax1.set_ylabel("log2(Fold-Change)")
ax1.set_xlabel("Taxon \n MWU = %s, pval = %s" % (mwustat, mwupval))
ax1.set_ylim(-2.7, 2)
ax1.set_title("roadmap_310 Enhancer Age Fold-change")
ax1.set_xlabel("")
ax1.set_title("")

plt.savefig("%sfig1b_1c-roadmap_310_WHOLE_ENHANCER_MRCA_DIST_TAXON2.pdf" % (RE), bbox_inches='tight')


# In[95]:


plot_ages = pd.concat([enh_ages, shuf_ages])

print(plot_ages.groupby("data_type")["mrca_2"].mean())


m, mp = stats.mannwhitneyu(plot_ages.loc[plot_ages.data_type == "shuffle", "mrca_2"],
                        plot_ages.loc[plot_ages.data_type == "roadmap_310", "mrca_2"])
print(m, mp)


# # landscape

#%%
enh_ = final_merge[["enh_id", "mrca_2", "core_remodeling"]]
#shuf_ = shuffle[["enh_id", "mrca_2"]]

enh_mrca = enh_.groupby(["mrca_2", "core_remodeling"])["enh_id"].count().reset_index()
enh_mrca.columns = ["mrca_2", "core_remodeling", "mrca_counts"]
enh_totals = enh_mrca.groupby(["core_remodeling"])["mrca_counts"].sum().reset_index()
enh_totals.columns = ["core_remodeling", "total_arch"]
enh_mrca = pd.merge(enh_mrca, enh_totals, how = "left", on = "core_remodeling")
enh_mrca["freq"] = enh_mrca.mrca_counts.divide(enh_mrca.total_arch)
#%%
enh_mrca
#%%
fig, ax = plt.subplots(figsize = (6,6))

sns.barplot(x = "mrca_2",
y = "freq", data = enh_mrca,
hue = "core_remodeling", palette = e_pal)

xlabs = [ "Prim", "Euar", "Bore", "Euth", "Ther", "Mam", "Amni", "Tetr", "Vert"]
ax.set(ylabel = "% total arch", xlabel = "")
ax.legend().remove()
ax.set_xticklabels(xlabs, rotation = 90)
plt.savefig("%sSupFig2.15D-roadmap_310_mrca_x_arch.pdf" % (RE), bbox_inches='tight')
