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

RE ="/dors/capra_lab/projects/enhancer_ages/fantom/results/for_publication/age_breaks/no_exon/"


# In[33]:


es_colors = [ "slate grey","greyish"]
es_pal = sns.xkcd_palette(es_colors)
sns.palplot(es_pal)


# In[34]:


path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/non-genic/"


# # load the genomic background

# In[4]:


syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t') # read the file
syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2", "mya", "mya2"]] # whittle down the df
syn_gen_bkgd


# # get enhancer file tissue/cell line descriptions

# In[5]:


desc_file = "/dors/capra_lab/data/fantom/fantom5/facet_expressed_enhancers/sample_id_descriptions.txt"
desc_df= pd.read_csv(desc_file, sep = '\t', header = None)


# # load dataframes

# In[36]:


enh = "%sFANTOM_NOEXON_enh_age_arch_full_matrix.tsv" % path
summaryEnh = "%sFANTOM_NOEXON_enh_age_arch_summary_matrix.tsv" % path

shuf = "%sSHUFFLE_NOEXON_FANTOM_enh_age_arch_full_matrix.tsv" % path
summaryShuf = "%sSHUFFLE_NOEXON_FANTOM_enh_age_arch_summary_matrix.tsv" % path

shuffle = pd.read_csv(shuf, sep = '\t')
shuffle.mrca_2 = shuffle.mrca_2.round(3)

final_merge = pd.read_csv(enh, sep = '\t')
final_merge.mrca_2 = final_merge.mrca_2.round(3)

# count of simple and complex enhancers per shuffled dataset.
shuf_arch = shuffle[["enh_id", "core_remodeling", "shuf_id"]].drop_duplicates()
shuf_arch_freq = shuf_arch.groupby(["core_remodeling", "shuf_id"])["enh_id"].count().reset_index()


totals = shuf_arch.groupby(["shuf_id"])["enh_id"].count().reset_index()
totals.columns = ["shuf_id", "totals"]
shuf_arch_freq = pd.merge(shuf_arch_freq, totals, how = "left")
shuf_arch_freq["freq"] = shuf_arch_freq["enh_id"].divide(shuf_arch_freq.totals)
shuf_arch_freq["dataset"] = "SHUFFLE"

# count of simple and complex enhancers in fantom
arch = final_merge[["enh_id", "core_remodeling"]].drop_duplicates()
arch_freq = arch.groupby("core_remodeling")["enh_id"].count().reset_index()
totals = len(arch)
arch_freq["freq"] = arch_freq["enh_id"].divide(totals)
arch_freq["dataset"] = "FANTOM"

archs = pd.concat([shuf_arch_freq, arch_freq]) # combine datasets for plotting


# In[37]:


final_merge.seg_index.min()


# In[ ]:


shuf_colors = [ "amber", "greyish",]
shuf_pal = sns.xkcd_palette(shuf_colors)
hue_order = ["FANTOM", "SHUFFLE"]

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
plt.savefig("%sfantom_enh_shuf_freq.pdf" %RE, bbox_inches = "tight")


# In[7]:


fig, ax = plt.subplots(figsize = (8, 8))
sns.barplot(x = "core_remodeling", y="freq", data = arch_freq, palette = palette)
ax.set_xticklabels(["Simple", "Complex\nEnhancer"])
ax.set_xlabel("")
ax.set_ylabel("Frequency of Dataset")
ax.set_title("FANTOM Enhancer Architectures")
for p in ax.patches:
    x=p.get_bbox().get_points()[:,0]
    y=p.get_bbox().get_points()[1,1]
    ax.annotate('{:.0f}%'.format(100.*y), (x.mean(), y),
            ha='left', va='bottom', color = "k", alpha = 0.4, fontsize = 20) # set the alignment of the text
ax.set_xticklabels("")
ax.set_xlabel("")
ax.set_title("")
#ax.get_legend().remove()
plt.savefig("%sfantom_enh_arch_freq.pdf" %RE, bbox_inches = "tight")


# # breaks

# In[38]:


breaks = final_merge.groupby(["enh_id", "core_remodeling"])["seg_index"].max().reset_index()
breaks.head()

sbreaks = shuffle.groupby(["enh_id", "core_remodeling", "shuf_id"])["seg_index"].max().reset_index()
sbreaks.head()


# In[39]:


sbreaks.min()


# In[40]:


# get shuffled breaks distribution
shuf_break = shuffle.groupby(["enh_id", "shuf_id"])["seg_index"].max().reset_index()

shuf_break_freq = shuf_break.groupby(["seg_index", "shuf_id"])["enh_id"].count().reset_index()

shuf_totals = shuf_break_freq.groupby("shuf_id")["enh_id"].sum().reset_index()
shuf_totals.columns = ["shuf_id", "shuf_id_totals"]

shuf_break_freq = pd.merge(shuf_break_freq, shuf_totals, how = "left", on = "shuf_id")

shuf_break_freq["freq"] = shuf_break_freq["enh_id"].divide(shuf_break_freq.shuf_id_totals)
shuf_break_freq["dataset"] = "Shuffle"


# In[41]:


sbreaks.seg_index.mean()


# In[42]:


breaks.seg_index.mean()


# # plot cumulative summary

# In[43]:


def get_cumsum(df):

    cdf= np.cumsum(df.freq)/1
    newdf = pd.DataFrame({"cdf": cdf})
    testdf = pd.merge(df, newdf, left_index = True, right_index = True)
    return testdf


# In[44]:


cumsum_dict = {}
for sid in shuf_break_freq.shuf_id.unique():
    testdf = shuf_break_freq.loc[shuf_break_freq.shuf_id == sid]
    newdf = get_cumsum(testdf)
    cumsum_dict[sid] = newdf
plot = pd.concat(cumsum_dict.values())
plot.head()


# In[45]:


breaks = final_merge.groupby(["enh_id"])["seg_index"].max().reset_index()
breaks_freq = breaks.groupby("seg_index")["enh_id"].count().reset_index()
enh_totals = len(breaks)

breaks_freq["freq"] = breaks_freq["enh_id"].divide(enh_totals)

breaks_freq["dataset"] = "FANTOM"
breaks_freq["shuf_id"] = "FANTOM"
breaks_freq["cdf"]= np.cumsum(breaks_freq.freq)/1

breaks_freq.head()


# In[46]:


shuf_break_freq.head()


# In[47]:


shuf_cdfplot = plot[["seg_index", "freq", "dataset", "shuf_id", "cdf"]]
shuf_cdfplot.columns = ['seg_index', 'freq', 'dataset', 'shuf_id', 'cdf']
#%%
concat = [breaks_freq, shuf_cdfplot]
plot_cdf = pd.concat(concat)
fig, ax = plt.subplots(figsize = (8,8))
x = "seg_index"
y = "cdf"


sns.lineplot(x, y, data = plot_cdf, ax = ax, hue = "dataset", palette = shuf_pal)
ax.set(xticks = (np.arange(0, plot_cdf.seg_index.max(), step = 5)), xlabel = "number of segments", ylabel = "cumulative distribution",
xlim = (0,10))
plt.savefig("%sFantom_CDF.pdf" %RE, bbox_inches = 'tight')


# In[48]:


plot.columns=["seg_index", "shuf_id", "shuf_count", "shuf_totals", "shuf_freq", "shuf_dataset", "shuf_cdf"]
plot.head()


# In[49]:


archs = pd.merge(plot, breaks_freq, how = "left", on = "seg_index" )
archs.head()


# # FET - More simple enhancers than expected?

# In[50]:


total_shuf_breaks = shuf_break_freq.groupby(["seg_index"])["enh_id"].sum().reset_index()

total_shuf_breaks.loc[total_shuf_breaks.seg_index ==1, "enh_id"][0]


# In[51]:


OR_dict = {}

for seg_index in breaks_freq.seg_index.unique():
    a = breaks_freq.loc[breaks_freq.seg_index ==seg_index, "enh_id"].item() # num simple enhancers
    b = enh_totals - a # num complex enhancers
    c = total_shuf_breaks.loc[total_shuf_breaks.seg_index ==seg_index, "enh_id"].item() # num simple shuffle
    d = total_shuf_breaks.enh_id.sum() - c # num complex shuffle

    obs = [[a,b], [c,d]]
    OR, P = stats.fisher_exact(obs)
    table = sm.stats.Table2x2(obs) # get confidence interval
    odds_ci = table.oddsratio_confint()
    newdf = pd.DataFrame({"seg_index":[seg_index], "a":[a], "b":[b], "c":[c], "d":[d],
                         "OR":[OR], "P":[P], "ci_lower" :[odds_ci[0]],
                        "ci_upper" :[odds_ci[1]]})
    OR_dict[seg_index] = newdf
    print(seg_index, obs, OR, P)


# In[52]:


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
ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(round(2**x, 1)))
ax.yaxis.set_major_formatter(ticks)
ax.yaxis.set_major_locator(MultipleLocator(0.5))
ax.set_ylim(-1.2,0.5)

plt.savefig("%sfig2b-fantom_age_seg_fold_change_matplotlib.pdf" %RE, bbox_inches = "tight")


# In[54]:


archs = pd.concat([plot, breaks_freq])

shuf_colors = [ "amber", "greyish","dusty purple", "windows blue","greyish"]
shuf_pal = sns.xkcd_palette(shuf_colors)

hue_order = ["FANTOM", "Shuffle"]
fig, ax = plt.subplots(figsize = (8,8))
sns.set("poster")
sns.barplot(x = "seg_index", y = "cdf", data = archs, hue = "dataset", hue_order = hue_order, palette = es_pal)
ax.set_xlabel("Number of Age Segments per Enhancer")
ax.set_ylabel("Cumulative Frequency of Dataset")
ax.set_title("FANTOM Enhancer Age Segment Count\n Cumulative distribution")
plt.legend()
ax.set_xlim(-1, 5.5)
#ax.set_xticklabels("")

ax.set_title("")
ax.get_legend().remove()
sns.set("poster")
plt.savefig("%sfig2b-fantom_age_seg_cum_dist_matplotlib.pdf" %RE, bbox_inches = "tight")


# In[55]:


enh_lens = final_merge.groupby(["enh_id", "core_remodeling","enh_len", "arch"])["mrca_2"].max().reset_index()

enh_lens = pd.merge(enh_lens, syn_gen_bkgd[[ "mrca_2","taxon2", "mya2"]], how = "left", on = "mrca_2")
enh_lens = enh_lens.drop_duplicates()
print(enh_lens.shape)

enh_lens["datatype"]="FANTOM"
enh_lens.head()


# In[56]:


shuf_len = shuffle.groupby(["enh_id", "core_remodeling", "enh_len", "shuf_id", "arch"])["mrca_2"].max().reset_index()
print(shuf_len.shape)

shuf_len["datatype"]="SHUFFLE"
shuf_len.mrca_2 =shuf_len.mrca_2.round(3)
shuf_len = pd.merge(shuf_len, syn_gen_bkgd[["mrca_2", "taxon2", "mya2"]], how = "left", on = "mrca_2")

shuf_len = shuf_len.drop_duplicates()
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

    b = pd.merge(b, syn_gen_bkgd[["mrca_2", "taxon2"]], how = "left", on = "mrca_2").drop_duplicates()
    b = b.fillna(0)

    return b


# In[58]:


lens_dict = {}

# enhancer frequencies

lens_dict["FANTOM"] = get_age_arch_freq(lens, "datatype", "FANTOM")
lens_dict["FANTOM"].dataset.unique()


# In[61]:


shuf_list = lens.shuf_id.unique()[1:]
for i in shuf_list:
    print(i)
    if len(i)> 8:
        new_i = (i.split("-")[1]).split("_")[-3]
    d = get_age_arch_freq(lens, "shuf_id", str(i))
    lens_dict[i] = d

for_stats = pd.concat(lens_dict.values())
fantom_stats = for_stats.loc[for_stats.dataset == "FANTOM"]
for_stats = for_stats.loc[for_stats.dataset != "FANTOM"]

m, p = stats.mannwhitneyu(fantom_stats.simple_freq, for_stats.simple_freq)
print(m, p)

m, p = stats.kruskal(fantom_stats.simple_freq.loc[fantom_stats["mrca_2"] == 0.000],
                     fantom_stats.simple_freq.loc[fantom_stats["mrca_2"] == 0.126],
                     fantom_stats.simple_freq.loc[fantom_stats["mrca_2"] == 0.131],
                     fantom_stats.simple_freq.loc[fantom_stats["mrca_2"] == 0.152],
                     fantom_stats.simple_freq.loc[fantom_stats["mrca_2"] == 0.175],
                     fantom_stats.simple_freq.loc[fantom_stats["mrca_2"] == 0.308],
                     fantom_stats.simple_freq.loc[fantom_stats["mrca_2"] == 0.38],
                     fantom_stats.simple_freq.loc[fantom_stats["mrca_2"] == 0.49],
                     fantom_stats.simple_freq.loc[fantom_stats["mrca_2"] == 0.656],
                     fantom_stats.simple_freq.loc[fantom_stats["mrca_2"] == 0.957],
                         )
print(m, p)

fantom_stats.simple_freq.mean()

for_stats.simple_freq.mean()


green = ["faded green"]
g_p = sns.xkcd_palette(green)
sns.palplot(g_p)
amber = ["amber"]
a_p = sns.xkcd_palette(amber)
sns.palplot(a_p)
grey = ["greyish"]
sg_p = sns.xkcd_palette(grey)
sns.palplot(sg_p)
slate = ["slate grey"]
sa_p = sns.xkcd_palette(slate)
sns.palplot(sa_p)

shuf = (pd.concat(lens_dict.values()))
shuf = shuf.loc[shuf.dataset.str.contains("shuf")]

fig, (ax, ax2) = plt.subplots( ncols = 2, figsize = (18,10))

sns.barplot(x = "taxon2", y  ="total_freq", data = lens_dict["FANTOM"].sort_values(by="mrca_2"), palette = g_p, ax = ax)
sns.barplot(x = "taxon2", y  = "simple_freq", data = lens_dict["FANTOM"].sort_values(by="mrca_2"), palette = a_p, ax = ax)

sns.barplot(x = "taxon2", y  ="total_freq", data = lens_dict["shuf-all_fantom_enh_112_tissues-96_age_breaks"].sort_values(by="mrca_2"), palette = sg_p, ax=ax2)
sns.barplot(x = "taxon2", y  = "simple_freq", data = shuf.sort_values(by="mrca_2"), palette = sa_p, ax = ax2)
ax2.set_xticklabels(ax2.get_xticklabels(), rotation = 90)
ax2.set_ylabel("Architecture Frequency ")
ax2.set_xlabel("MRCA")


ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
ax.set_ylabel("Architecture Frequency ")
ax.set_xlabel("MRCA")


ax.set_xlabel("")
ax.set_title("")
#ax.get_legend().remove()

ax2.set_xlabel("")
ax2.set_title("")
#ax2.get_legend().remove()
plt.savefig("%sfig2F-age_x_arch_freq.pdf" % RE, bbox_inches = "tight")


# # enhancer length

# In[62]:


##### fig, ax = plt.subplots(figsize = (10,10))
fig, ax = plt.subplots(figsize = (8, 8))
sns.boxplot(x = "core_remodeling", y = "enh_len",
            data = enh_lens, palette = palette,
           showfliers = False)
ax.set_xticklabels(["Simple", "Complex"])
ax.set_xlabel("Architecture")
ax.set_ylabel("Enhancer Length (bp)")
ax.legend(bbox_to_anchor=(1.0, 1.0))


# Get medians

lens_data = lens.groupby(["core_remodeling", "datatype"])["enh_len"].median().reset_index()
lens_data_mean = lens.groupby(["core_remodeling", "datatype"])["enh_len"].mean().round().reset_index()
print("medians", lens_data,)
print("means", lens_data_mean,)
ms, msp = stats.mannwhitneyu(enh_lens.enh_len.loc[enh_lens.core_remodeling==1],
                            enh_lens.enh_len.loc[enh_lens.core_remodeling==0])
print("simple", ms, msp)
ax.set_xlabel("")
ax.set_title("")
#ax.get_legend().remove()
plt.savefig("%sfantom_enh_lens.pdf" %RE, bbox_inches = "tight")


# In[63]:


##### fig, ax = plt.subplots(figsize = (10,10))
fig, ax = plt.subplots(figsize = (8, 8))
sns.boxplot(x = "core_remodeling", y = "enh_len", hue = "datatype",
            data = lens, palette = es_pal,
           showfliers = False, notch = True)
ax.set_xticklabels(["Simple", "Complex"])
ax.set_xlabel("Architecture")
ax.set_ylabel("Enhancer Length (bp)")
ax.legend(bbox_to_anchor=(1.0, 1.0))

# Get medians

lens_data = lens.groupby(["core_remodeling", "datatype"])["enh_len"].median().reset_index()
lens_data_mean = lens.groupby(["core_remodeling", "datatype"])["enh_len"].mean().reset_index()
print("medians", lens_data,)
print("means", lens_data_mean,)

ax.set_xlabel("")
ax.set_title("")
ax.get_legend().remove()
plt.savefig("%sfantom_enh_lens.pdf" %RE, bbox_inches = "tight")


# In[64]:


stats.mannwhitneyu(lens.loc[(lens.mrca_2>0.175) & (lens.datatype == "FANTOM"),"enh_len"],
                  lens.loc[(lens.mrca_2<=0.175) & (lens.datatype == "FANTOM"),"enh_len"],)

lens.loc[lens.mrca_2>0.175].groupby("datatype")["enh_len"].median()

lens.loc[lens.mrca_2<=0.175].groupby("datatype")["enh_len"].median()

stats.mannwhitneyu(lens.loc[(lens.mrca_2<=0.175) & (lens.datatype == "FANTOM"),"enh_len"],
                  lens.loc[(lens.mrca_2<=0.175) & (lens.datatype == "SHUFFLE"),"enh_len"],)


# In[65]:


hue_order = ["FANTOM", "SHUFFLE"]
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
plt.savefig("%sfig1d-fantom_enh_mrca_lens.pdf" %RE, bbox_inches = "tight")


# In[66]:


stats.mannwhitneyu(lens.loc[(lens.arch == "simple") & (lens.datatype == "FANTOM"), "enh_len"],
                  lens.loc[(lens.arch == "simple") & (lens.datatype == "SHUFFLE"), "enh_len"], )


# In[67]:


enh_lens.head()


# In[68]:


#MRCA LINEAR REGRESSION
xs = enh_lens.mrca_2.loc[enh_lens.arch.str.contains("simple")]
ys = enh_lens.enh_len.loc[enh_lens.arch.str.contains("simple")]

slope_s, intercept_s, r_value_s, p_value_s, std_err_s = stats.linregress(xs,ys)

xc = enh_lens.mrca_2.loc[enh_lens.arch.str.contains("complex")]
yc = enh_lens.enh_len.loc[enh_lens.arch.str.contains("complex")]

slope_c, intercept_c, r_value_c, p_value_c, std_err_c  = stats.linregress(xc,yc)
slope_c, intercept_c, r_value_c, p_value_c, std_err_c

xsshuf = shuf_len.mrca_2.loc[shuf_len.arch.str.contains("simple")]
ysshuf = shuf_len.enh_len.loc[shuf_len.arch.str.contains("simple")]

slope_ss, intercept_ss, r_value_ss, p_value_ss, std_err_ss = stats.linregress(xsshuf,ysshuf)

xcs = shuf_len.mrca_2.loc[shuf_len.arch.str.contains("complex")]
ycs = shuf_len.enh_len.loc[shuf_len.arch.str.contains("complex")]

slope_cs, intercept_cs, r_value_cs, p_value_cs, std_err_cs = stats.linregress(xcs,ycs)

print("Simple enhancer slope = %s, r2 = %s, pval = %s\nComplex enhancer slope = %s, r2 = %s, pval = %s"              % (round(slope_s, 0), round(r_value_s,2), p_value_s, round(slope_c,0),  round(r_value_c,2), p_value_c))
print("SHUFFLE Simple enhancer slope = %s, r2 = %s, pval = %s\nSHUFFLE Complex enhancer slope = %s, r2 = %s, pval = %s"              % (round(slope_ss, 0), round(r_value_ss,2), p_value_ss, round(slope_cs,0),  round(r_value_cs,2), p_value_cs))


# In[69]:


enh_lens.head()


# In[70]:


# MYA LINEAR REGRESSION
xs = enh_lens.mya2.loc[enh_lens.arch.str.contains("simple")]
ys = enh_lens.enh_len.loc[enh_lens.arch.str.contains("simple")]

slope_s, intercept_s, r_value_s, p_value_s, std_err_s = stats.linregress(xs,ys)

xc = enh_lens.mya2.loc[enh_lens.arch.str.contains("complex")]
yc = enh_lens.enh_len.loc[enh_lens.arch.str.contains("complex")]

slope_c, intercept_c, r_value_c, p_value_c, std_err_c  = stats.linregress(xc,yc)
slope_c, intercept_c, r_value_c, p_value_c, std_err_c

xsshuf = shuf_len.mya2.loc[shuf_len.arch.str.contains("simple")]
ysshuf = shuf_len.enh_len.loc[shuf_len.arch.str.contains("simple")]

slope_ss, intercept_ss, r_value_ss, p_value_ss, std_err_ss = stats.linregress(xsshuf,ysshuf)

xcs = shuf_len.mya2.loc[shuf_len.arch.str.contains("complex")]
ycs = shuf_len.enh_len.loc[shuf_len.arch.str.contains("complex")]

slope_cs, intercept_cs, r_value_cs, p_value_cs, std_err_cs = stats.linregress(xcs,ycs)

print("Simple enhancer slope = %s, r2 = %s, pval = %s\nComplex enhancer slope = %s, r2 = %s, pval = %s"              % (round(slope_s, 3), round(r_value_s,2), p_value_s, round(slope_c,3),  round(r_value_c,2), p_value_c))
print("SHUFFLE Simple enhancer slope = %s, r2 = %s, pval = %s\nSHUFFLE Complex enhancer slope = %s, r2 = %s, pval = %s"              % (round(slope_ss, 3), round(r_value_ss,2), p_value_ss, round(slope_cs,3),  round(r_value_cs,2), p_value_cs))


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

hue_order = ["FANTOM", "Shuffle"]
fig,(ax1) = plt.subplots(figsize = (8, 8))
order =["Simple", "Complexenh"]
sns.barplot(y = "enh_len", x = "taxon2", data = enh_lens.sort_values(by = "mrca_2"), ax = ax1,            hue = "arch",  palette = e_pal, estimator = np.median)#showfliers=False)


ms, msp = stats.mannwhitneyu(enh_lens.enh_len.loc[enh_lens.arch.str.contains("imple")],
                            shuf_len.enh_len.loc[shuf_len.arch.str.contains("imple")])
print("simple", ms, msp)

mc, mcp = stats.mannwhitneyu(enh_lens.enh_len.loc[enh_lens.arch.str.contains("omplex")],
                            shuf_len.enh_len.loc[shuf_len.arch.str.contains("omplex")])
print("complex", mc, mcp)
ax1.set_ylabel("Enhancer Length (bp)")
ax1.set_xticklabels(ax1.get_xticklabels(), rotation = 90, horizontalalignment = "left")
ax1.set_ylim(190,400)
ax1.set_xlabel("")
ax1.set_title("")
ax1.get_legend().remove()
plt.savefig("%sfig2c-Fantom_ENH_MRCA_x_LEN_ENH.pdf" % RE, bbox_inches = "tight")


# In[74]:


len_concat = pd.concat([enh_lens, shuf_len])

len_concat.head()

len_concat["datatype2"] = len_concat.datatype + "-"+ len_concat.arch
len_concat["datatype2"].unique()


e_colors = [ "amber","greyish", "faded green", "slate grey"]
e_pal = sns.xkcd_palette(e_colors)
hue_order = ['FANTOM-simple','SHUFFLE-simple','FANTOM-complexenh','SHUFFLE-complexenh']


fig,(ax1) = plt.subplots(figsize = (8, 8))
order =["Simple", "Complexenh"]
sns.barplot(y = "enh_len", x = "taxon2", data = len_concat.sort_values(by = "mrca_2"), ax = ax1,            hue = "datatype2",  palette = e_pal,
            hue_order = hue_order,
            estimator = np.median)#showfliers=False)


ax1.set_ylabel("Enhancer Length (bp)")
ax1.set_xticklabels(ax1.get_xticklabels(), rotation = 90, horizontalalignment = "left")

ax1.set_xlabel("")
ax1.set_title("")
ax1.legend().remove()
plt.savefig("%sfigS2.7-Fantom_ENH_MRCA_x_LEN_ENH.pdf" % RE, bbox_inches = "tight")


# In[75]:


e_colors = [ "slate grey","greyish",]
e_pal = sns.xkcd_palette(e_colors)
hue_order = ['FANTOM', 'SHUFFLE']

order =["Simple", "Complexenh"]
sns.lmplot(y = "enh_len", x = "mya2", data = len_concat,            hue = "datatype",  palette = e_pal,
            hue_order = hue_order, scatter = False, x_estimator = np.median)
plt.savefig("%sfigS1.2-LM_enh_len.pdf" % RE, bbox_inches = "tight")


# In[76]:


e_colors = [ "amber","greyish", "faded green", "slate grey"]
e_pal = sns.xkcd_palette(e_colors)
hue_order = ['FANTOM-simple','SHUFFLE-simple','FANTOM-complexenh','SHUFFLE-complexenh']
#fig,(ax1) = plt.subplots(figsize = (8, 8))
order =["Simple", "Complexenh"]
sns.lmplot(y = "enh_len", x = "mya2", data = len_concat,            hue = "datatype2",  palette = e_pal,
            hue_order = hue_order, scatter = False)
plt.savefig("%sfigS2.7-LM_enh_len.pdf" % RE, bbox_inches = "tight")


# In[77]:


enh_lens.groupby("core_remodeling")["mrca_2"].mean()


# In[78]:


shuf_len.groupby("core_remodeling")["mrca_2"].mean()


# # enhancer age segments

# In[79]:


plot_len= pd.concat([enh_lens, shuf_len])

plot_len.groupby(["arch"]).count()


# In[80]:


plot_len.head()


# In[81]:


ms, msp = stats.mannwhitneyu(enh_lens.enh_len.loc[enh_lens.arch.str.contains("imple")],
                            shuf_len.enh_len.loc[shuf_len.arch.str.contains("imple")])
print("simple", ms, msp)

mc, mcp = stats.mannwhitneyu(enh_lens.enh_len.loc[enh_lens.arch.str.contains("omplex")],
                            shuf_len.enh_len.loc[shuf_len.arch.str.contains("omplex")])
print("complex", mc, mcp)


print("Simple enhancer slope = %s, r2 = %s, pval = %s\nComplex enhancer slope = %s, r2 = %s, pval = %s"              % (round(slope_s, 0), round(r_value_s,2), p_value_s, round(slope_c,0),  round(r_value_c,2), p_value_c))
print("SHUFFLE Simple enhancer slope = %s, r2 = %s, pval = %s\nSHUFFLE Complex enhancer slope = %s, r2 = %s, pval = %s"              % (round(slope_ss, 0), round(r_value_ss,2), p_value_ss, round(slope_cs,0),  round(r_value_cs,2), p_value_cs))


# In[82]:


core_breaks = final_merge.groupby(["enh_id", "core_remodeling"])["mrca_2", "seg_index"].max().reset_index()

core_breaks = pd.merge(core_breaks, syn_gen_bkgd[["mrca_2", "taxon2"]])
core_breaks["datatype"] = "FANTOM"
core_breaks=core_breaks.drop_duplicates()

smrca = shuffle.groupby(["enh_id", "core_remodeling", "shuf_id"])["mrca_2"].max().reset_index()

# DATA CLEAN UP - remove the enh_id errors in the shuffle db where the oldest complexenh age is homosapiens
smrca_list = list(smrca.loc[(smrca.core_remodeling == 1)                             & (smrca.mrca_2 == 0.00)]["enh_id"])

smrca = smrca[~smrca.enh_id.isin(smrca_list)]

sseg = shuffle.groupby(["enh_id", "core_remodeling", "shuf_id"])["seg_index"].max().reset_index()

shuf_core_breaks = pd.merge(smrca, sseg, how ="left", on = ("enh_id", "core_remodeling", "shuf_id"))

shuf_core_breaks = pd.merge(shuf_core_breaks, syn_gen_bkgd[["mrca_2", "taxon2"]])
shuf_core_breaks["datatype"] = "SHUFFLE"
shuf_core_breaks=shuf_core_breaks.drop_duplicates()


plot_breaks=pd.concat([core_breaks,shuf_core_breaks])

f, ax = plt.subplots(figsize = (8,8))

hue_order = ["FANTOM", "SHUFFLE"]
sns.barplot(x = "taxon2", y = "seg_index",
            data =plot_breaks[plot_breaks.core_remodeling ==1].sort_values(by = "mrca_2"),
            palette = cs_pal, hue = "datatype", hue_order = hue_order)

ax.set_xticklabels(ax.get_xticklabels(), rotation = 90, horizontalalignment = "left")
ax.set_ylabel("# age segments")
ax.set_title("FANTOM age segments v MRCA")
ax.set_xlabel("")
ax.set_title("")
ax.get_legend().remove()
ax.ylim = (0,3.2)
plt.savefig("%sFigS2.4B-FANTOM_age_segments_v_MRCA_bar.pdf" %RE, bbox_inches = "tight")


# In[83]:


f, ax = plt.subplots(figsize = (8,8))

hue_order = ["FANTOM", "SHUFFLE"]
sns.barplot(x = "datatype", y = "seg_index",
            data =plot_breaks[plot_breaks.core_remodeling ==1],
            palette = cs_pal)

ax.set_xticklabels(ax.get_xticklabels(), rotation = 90, horizontalalignment = "left")
ax.set_ylabel("# age segments")
ax.set_title("FANTOM v Shuffle age segments")
ax.set_xlabel("")
ax.set_title("")
ax.legend().remove()
ax.ylim = (0,3.2)
plt.savefig("%sFigS2.4A-FANTOM_age_segments_bar.pdf" %RE, bbox_inches = "tight")


# In[84]:


from scipy.stats import mstats
kw_list = []
fantom_list = []
shuffle_list = []
for i in plot_breaks.mrca_2.unique():
    for d in plot_breaks.datatype.unique():
        mrca_list = plot_breaks.loc[(plot_breaks.mrca_2 == i)
                                    & (plot_breaks.datatype == d)
                                    & (plot_breaks.core_remodeling ==1) , "seg_index"].to_list() # collect the # segments per age
        kw_list.append(mrca_list)
        if d == "FANTOM":
            fantom_list.append(mrca_list) # collect the # segments per age fantom
        elif d == "SHUFFLE":
            shuffle_list.append(mrca_list)# collect the # segments per age shuffle

args=[l for l in kw_list]
argsF=[l for l in fantom_list]
argsS=[l for l in shuffle_list]
stats.mstats.kruskalwallis(*args) # KW on age segment differences for all age and datatype complex enhancers

stats.mstats.kruskalwallis(*argsF) # KW on age segment differences among FANTOM complex enhancers

stats.mstats.kruskalwallis(*argsS) # KW on age segment differences among Shuffle complex enhancers


# In[85]:


# MWU on age segments in FANTOM v. Shuffle
stats.mannwhitneyu(plot_breaks.loc[(plot_breaks.datatype == "FANTOM") & (plot_breaks.core_remodeling ==1), "seg_index"],
                  plot_breaks.loc[(plot_breaks.datatype == "SHUFFLE") & (plot_breaks.core_remodeling ==1), "seg_index"],)


# In[86]:


plot_breaks.groupby(["datatype", "core_remodeling"])["seg_index"].mean()


# In[87]:


plot_breaks.groupby(["datatype"])["seg_index"].mean()


# # Enhancer age fold-change

# In[88]:


enh_ages = final_merge.groupby(['enh_id'])["mrca_2"].max().reset_index()
enh_ages["data_type"] = "fantom"
shuf_ages = shuffle.groupby(['enh_id', 'shuf_id',])["mrca_2"].max().reset_index()
shuf_ages["data_type"] = "shuffle"

shuf_ages.head()


# In[89]:


ANNOTATIONS = 0

shuffle_dist = shuffle.groupby(["enh_id", "shuf_id"])["mrca_2"].max().reset_index() # calculate age distribution shuffled enhancers by enh_id

shuffle_dist["mrca_2"] = shuffle_dist["mrca_2"].round(3) # round ages


shuffle_dist2 = shuffle_dist.groupby(["mrca_2", "shuf_id"])["enh_id"].count().reset_index()
shuffle_dist2.columns = ["mrca_2", "shuf_id", "mrca_count"]

totals = shuffle_dist2.groupby(["shuf_id"])["mrca_count"].sum().reset_index()
totals.columns = ["shuf_id", "id_totals"]

shuffle_dist2 = pd.merge(shuffle_dist2, totals, how = "left")
shuffle_dist2['freq'] = shuffle_dist2.mrca_count.divide(shuffle_dist2.id_totals)
shuffle_dist2["dataset"] = "shuffle"


# In[90]:


plot_dist = final_merge.groupby("enh_id")["mrca_2"].max().reset_index() # calculate age distribution of FANTOM enhancers by enh_id

plot_dist["mrca_2"] = plot_dist["mrca_2"].round(3) # round ages

totals = plot_dist.groupby("mrca_2")["enh_id"].count().reset_index()
totals.columns = ["mrca_2", "mrca_count"]
totals['freq'] = totals.mrca_count.divide(totals.mrca_count.sum())

plot_dist = pd.merge(plot_dist, totals, how = "left", on = "mrca_2")
plot_dist["shuf_id"] = "FANTOM"
plot_dist["dataset"] = "FANTOM"
dist = pd.concat([shuffle_dist2, plot_dist]) # concat shuffle and fantom age distributions
dist.head()


# In[91]:


dist = pd.merge(dist, syn_gen_bkgd, how = "left", on = "mrca_2")


# calculate MWU shuffle v. FANTOM age distributions (from each enhancer)

mwustat, mwupval = stats.mannwhitneyu(shuffle_dist2["mrca_2"].loc[shuffle_dist2["dataset"] =="shuffle"],
                                     plot_dist["mrca_2"].loc[plot_dist["dataset"] =="FANTOM"] )
dist = dist[["taxon2", "mrca_2", "shuf_id", "dataset", "freq"]].drop_duplicates()

dist.sort_values(by = "mrca_2").head()


# In[92]:


# Fold change of fantom/shuffle age frequency
fold_change_dict = {}
for i in dist.shuf_id.unique():
    if i != "FANTOM":
        fold_change = dist.loc[(dist.shuf_id.str.contains(i))|(dist.shuf_id == "FANTOM")].pivot(index="mrca_2", columns='dataset')['freq'].reset_index()
        fold_change = pd.merge(fold_change, syn_gen_bkgd[["mrca_2", "taxon2"]], how = "left").drop_duplicates()
        fold_change["fold_change"] = np.log2(fold_change["FANTOM"].divide(fold_change["shuffle"]))
        fold_change_dict[i] = fold_change

fc = pd.concat(fold_change_dict.values())


# In[93]:


dist.head()


# In[94]:


colors = ["greyish", "slate grey",]
palette = sns.xkcd_palette(colors)
sns.set("poster")
#
mwustat, mwupval = stats.mannwhitneyu(shuffle_dist2["mrca_2"].loc[shuffle_dist2["dataset"] =="shuffle"],
                                     plot_dist["mrca_2"].loc[plot_dist["dataset"] =="FANTOM"] )
fig, ( ax2, ax1) = plt.subplots(ncols = 2, figsize=(16,8))

# line plot

sns.barplot(x = "taxon2", y= "freq",  data = dist.sort_values(by="mrca_2"),              hue = "dataset", palette = palette, ax = ax2)

ax2.set_xticklabels(ax2.get_xticklabels(), rotation = 90, horizontalalignment="left")
ax2.set_ylabel("Frequency")
ax2.set_xlabel("Taxon \n MWU = %s, pval = %s" % (mwustat, mwupval))
ax2.set_ylim(0, 0.6)
ax2.set_title("FANTOM Enhancer Age Distribution")
ax2.set_xlabel("")
ax2.set_title("")
ax2.legend(frameon = False)

# fold change plot
colors = ["slate grey"]
palette_fc = sns.xkcd_palette(colors)
sns.barplot(x = "taxon2", y= "fold_change",  data = fc.sort_values(by="mrca_2"),              palette = palette_fc, ax = ax1)
if ANNOTATIONS == 1:
    for p in ax1.patches:
        x = p.get_x().round(2)
        y = (p.get_height().round(2))
        ax1.annotate('{:.1f}'.format(p.get_height()), (x, y), color = "black", fontsize = 20)
ax1.set_xticklabels(ax1.get_xticklabels(), rotation = 90, horizontalalignment="left")

ax1.set_ylabel("log2(Fold-Change)")
ax1.set_xlabel("Taxon \n MWU = %s, pval = %s" % (mwustat, mwupval))
ax1.set_ylim(-2.7, 2)
ax1.set_title("FANTOM Enhancer Age Fold-change")
ax1.set_xlabel("")
ax1.set_title("")

plt.savefig("%sfig1b_1c-FANTOM_WHOLE_ENHANCER_MRCA_DIST_TAXON2.pdf" % (RE), bbox_inches='tight')


# In[95]:


plot_ages = pd.concat([enh_ages, shuf_ages])

print(plot_ages.groupby("data_type")["mrca_2"].mean())


m, mp = stats.mannwhitneyu(plot_ages.loc[plot_ages.data_type == "shuffle", "mrca_2"],
                        plot_ages.loc[plot_ages.data_type == "fantom", "mrca_2"])
print(m, mp)


# # landscape

# In[96]:


ssyn_lens = shuffle[["syn_id", "code","syn_len", "mrca_2"]].drop_duplicates()
ssyn_lens.head()


# In[97]:


ssyn_lens = shuffle[["syn_id", "code","syn_len", "mrca_2"]].drop_duplicates()
ssyn_lens["datatype"] = "SHUFFLE"
order = ["simple", "complex_core", "derived"]
fig, ax = plt.subplots(figsize = (8,8))

sns.boxplot(x = "code", y = "syn_len",
            data = ssyn_lens, palette = es_pal, order = order,
           showfliers = False)
#ax.set_xticklabels(["Simple", "Complexenh"])
ax.set_xlabel("Architecture")
ax.set_ylabel("Enhancer Length (bp)")
ax.set_xticklabels(["Simple", "Complex\nCore", "Derived"])
ax.set_xlabel("")
ax.set_title("")
plt.savefig("%sfantom_shuffle_syn_lens.pdf" %RE, bbox_inches = "tight")

ssyn_lens_df = ssyn_lens.groupby(["code"])["syn_len"].median().reset_index()
print(ssyn_lens_df)


# In[98]:


syn_lens = final_merge[["syn_id", "code","syn_len", "mrca_2"]].drop_duplicates()

order = ["simple", "complex_core", "derived"]
fig, ax = plt.subplots(figsize = (8,8))

sns.boxplot(x = "code", y = "syn_len",
            data = syn_lens, palette = arch_palette, order = order,
           showfliers = False)
#ax.set_xticklabels(["Simple", "Complexenh"])
ax.set_xlabel("Architecture")
ax.set_ylabel("Enhancer Length (bp)")
ax.set_xlabel("")
ax.set_title("")
ax.set_xticklabels(["Simple", "Complex\nCore", "Derived"])
plt.savefig("%sfantom_syn_lens.pdf" %RE, bbox_inches = "tight")

syn_lens_df = syn_lens.groupby(["code"])["syn_len"].median().reset_index()
print(syn_lens_df)


# In[99]:


# manipulate enh_lens_df to concat with syn_lens_df
enh_lens["code"] = ""
enh_lens["code"].loc[enh_lens["core_remodeling"] == 0] = "simple"
enh_lens["code"].loc[enh_lens["core_remodeling"] == 1] = "complexenh"
enh_lens_df = enh_lens[["enh_id", "code", "mrca_2", "enh_len"]]
enh_lens_df.columns = ["shuf_id", "code", "mrca_2",  "len"]

enh_lens_df.head()


# In[100]:


# rename syn_lens columns to concat w/ enh_lens_df
syn_lens.columns = ["shuf_id", "code",  "len", "mrca_2"]

syn_lens.head()


# In[101]:


# concat enh_lens and syn_lens together, plot on single plot
lens = pd.concat([enh_lens_df, syn_lens]).drop_duplicates()
lens = pd.merge(lens, syn_gen_bkgd[["mrca_2", "taxon2"]], how = "left", on = "mrca_2")
lens.sort_values(by="mrca_2").head(20)


# In[102]:


order = ["simple", "complexenh", "complex_core", "derived"]
fig, ax = plt.subplots(figsize = (8,8))

sns.boxplot(x = "code", y = "len",
            data = lens, palette = palette,
            showfliers = False, order = order)
ax.set_title("All Fantom Architecture Lengths")
ax.set_xlabel("Architecture")
ax.set_ylabel("Enhancer Length (bp)")
ax.set_xlabel("")
ax.set_title("")
ax.set_xticklabels(["Simple", "Complex\nEnhancer", "Complex\nCore", "Derived"])
#plt.savefig("%sfantom_enh_and_syn_lens.pdf" %RE, bbox_inches = "tight")


# In[103]:


final_merge["seg_rep"] = final_merge.syn_len.divide(final_merge.enh_len)
shuffle["seg_rep"] = shuffle.syn_len.divide(shuffle.enh_len)
sns.distplot(final_merge["seg_rep"])


# In[104]:


percent_enh = final_merge.groupby(["enh_id", "code"])["seg_rep"].sum().reset_index()
percent_shuf = shuffle.groupby(["enh_id", "code"])["seg_rep"].sum().reset_index()
percent_enh.head()


# In[126]:


#####FANTOM#####
# group by age & architecture
syn_counts= final_merge.groupby(["mrca", "code"])["syn_id"].count().reset_index()

# group by architecture
totals = syn_counts.groupby(["code"])["syn_id"].sum().reset_index()
totals.columns = ["code", "total_count"]

# calculate age v. architecture
concat2 = pd.merge(syn_counts, totals, how = "left", on = "code")
concat2["frequency"] = concat2["syn_id"].divide(concat2["total_count"])
concat2["dataset"] = "FANTOM"


#####SHUFFLE#####
# group by age & architecture
shuf_arch= shuffle.groupby(["mrca", "code"])["syn_id"].count().reset_index()

# group by architecture
shuf_totals = shuf_arch.groupby(["code"])["syn_id"].sum().reset_index()
shuf_totals.columns = ["code", "total_count"]

# calculate age v. architecture
shuf_concat = pd.merge(shuf_arch, shuf_totals, how = "left", on = "code")
shuf_concat["frequency"] = shuf_concat["syn_id"].divide(shuf_concat["total_count"])
shuf_concat["dataset"] = "shuffle"
concat2 = pd.concat([concat2, shuf_concat])
concat2.mrca = concat2.mrca.round(3)
concat_bkgd2=pd.merge(concat2, syn_gen_bkgd, how = "left", on = "mrca")

#####CONCAT2#####
concat_bkgd2["code2"] = concat_bkgd2["code"] + "-"+ concat_bkgd2["dataset"]

#####df_empty##### add in empty data
df_empty = pd.DataFrame({'mrca' :[0.957, 0.957, 0.000, 0.000], 'code':["derived", "derived", "complex_core", "complex_core"],
 'syn_id':[0, 0, 0, 0],
 'total_count':[0, 0, 0, 0],
 'frequency':[0, 0, 0, 0],
 'dataset':["FANTOM", "shuffle", "FANTOM", "shuffle"],
 'taxon':["Vertebrata (615)", "Vertebrata (615)", "Homo sapien (0)", "Homo sapien (0)"],
 'mrca_2':[0.957, 0.957, 0.000, 0.000],
 'taxon2':["Vertebrata (615)", "Vertebrata (615)", "Homo sapien (0)", "Homo sapien (0)"],
 'code2':["derived-FANTOM", "derived-shuffle", "complex_core-FANTOM", "complex_core-shuffle"]})
concat_bkgd2 = concat_bkgd2.append(df_empty, ignore_index=True)
concat_bkgd2.sort_values(by = "mrca").head()

tx2 = concat_bkgd2.groupby(['mrca_2', 'code', 'taxon2', 'dataset' , "code2"])["frequency"].sum().reset_index()

tx2.head(20)


# In[127]:


# what are the frequencies of architectures in each enhancer type?
matharch = tx2.loc[tx2.dataset =="FANTOM"].drop_duplicates()
matharch


# In[128]:


# frequency of simple enhancers younger than eutherian
matharch.loc[(matharch.code.str.contains("simple")) & (matharch.mrca_2<0.175)]["frequency"].sum()


# In[129]:


# frequency of simple enhancers of eutherian age
matharch.loc[(matharch.code.str.contains("simple")) & (matharch.mrca_2==0.175)]["frequency"].sum()


# In[130]:


# frequency of simple enhancers older than eutherian
matharch.loc[(matharch.code.str.contains("simple")) & (matharch.mrca_2>0.175)]["frequency"].sum()


# In[131]:


arch_colors = [ "windows blue","greyish", "amber", "yellowish", "dusty purple",  "lavender"]
arch_palette = sns.xkcd_palette(arch_colors)
fig, (ax3, ax2, ax1) = plt.subplots(ncols = 3, figsize=(27,8))

#####DERIVED#####
sns.pointplot(x = "taxon2", y = "frequency",
              data = tx2.loc[tx2["code"].str.contains("derived")].sort_values(by=["mrca_2", "code"]),\
              hue = "code2", ax = ax1,  palette = arch_palette)
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90, horizontalalignment = "left", fontsize = 20)
mwu_stat, pval = stats.mannwhitneyu(final_merge["mrca_2"].loc[(final_merge["code"].str.contains("derived"))],                   shuffle["mrca_2"].loc[(shuffle["code"].str.contains("derived"))])
ax1.set_xlabel("mrca\nMWU = %s, pval = %s" %(mwu_stat, pval))
ax1.set_title("Age Distribution")
ax1.set_ylabel("Frequency")
ax1.legend().remove()
ax1.set_xlabel("")
ax1.set_title("")
#####COMPLEX_CORE#####
arch_colors = ["dusty purple",  "greyish"]
arch_palette = sns.xkcd_palette(arch_colors)
sns.pointplot(x = "taxon2", y = "frequency",
              data = tx2.loc[tx2["code"].str.contains("complex_core")].sort_values(by=["mrca_2", "code"]),\
              hue = "code2", ax = ax2,  palette = arch_palette)
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=90, horizontalalignment = "left", fontsize = 20)
mwu_stat, pval = stats.mannwhitneyu(final_merge["mrca_2"].loc[(final_merge["code"].str.contains("complex_core"))],                   shuffle["mrca_2"].loc[(shuffle["code"].str.contains("complex_core"))])
ax2.set_xlabel("mrca\nMWU = %s, pval = %s" %(mwu_stat, pval))
ax2.set_title("Age Distribution")
ax2.set_ylabel("Frequency")
ax2.legend().remove()
ax2.set_xlabel("")
ax2.set_title("")
#####SIMPLE_CORE#####
arch_colors = ["amber", "greyish"]
arch_palette = sns.xkcd_palette(arch_colors)
sns.pointplot(x = "taxon2", y = "frequency",
              data = tx2.loc[tx2["code"].str.contains("simple")].sort_values(by=["mrca_2", "code"]),\
              hue = "code2", ax = ax3,  palette = arch_palette)
ax3.set_xticklabels(ax3.get_xticklabels(), rotation=90, horizontalalignment = "left", fontsize = 20)
mwu_stat, pval = stats.mannwhitneyu(final_merge["mrca_2"].loc[(final_merge["code"].str.contains("simple"))],                   shuffle["mrca_2"].loc[(shuffle["code"].str.contains("simple"))])
ax3.set_xlabel("mrca\nMWU = %s, pval = %s" %(mwu_stat, pval))
ax3.set_title("Age Distribution")
ax3.set_ylabel("Frequency")
ax3.legend().remove()
ax3.set_xlabel("")
ax3.set_title("")

plt.savefig("%sMRCA_fantom_v_shuffle_core_derived_grouped2.pdf"% (RE))


# In[132]:


pivot2 = pd.pivot_table(tx2, values='frequency',                           index=['mrca_2', 'code', 'taxon2'], columns=['dataset']).reset_index()
pivot2["fold_change"] = pivot2["FANTOM"].divide(pivot2["shuffle"]) -1
pivot2=pivot2.fillna(0)

derived = [ "windows blue","greyish",]
derivedp = sns.xkcd_palette(derived)
pivot2.head()

fig, (ax3, ax2, ax1) = plt.subplots(ncols = 3, figsize=(27,8))

#####DERIVED#####
sns.barplot(x = "taxon2", y = "fold_change",            data = pivot2.loc[pivot2["code"].str.contains("derived")].sort_values(by=["mrca_2", "code"]),              hue = "code", ax = ax1,  palette = derivedp)

ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90, horizontalalignment = "left", fontsize = 20)

mwu_stat, pval = stats.mannwhitneyu(final_merge["mrca"].loc[(final_merge["code"].str.contains("derived"))],                   shuffle["mrca"].loc[(shuffle["code"].str.contains("derived"))])

ax1.set_xlabel("mrca\nMWU = %s, pval = %s" %(mwu_stat, pval))
ax1.set_ylabel("Fold Change")
ax1.legend().remove()

#####COMPLEX_CORE#####
complexcore = ["dusty purple",  "greyish"]
complexcorep = sns.xkcd_palette(complexcore)
sns.barplot(x = "taxon2", y = "fold_change",            data = pivot2.loc[pivot2["code"].str.contains("complex_core")].sort_values(by=["mrca_2", "code"]),              hue = "code", ax = ax2,  palette = complexcorep)
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=90, horizontalalignment = "left", fontsize = 20)
mwu_stat, pval = stats.mannwhitneyu(final_merge["mrca"].loc[(final_merge["code"].str.contains("complex_core"))],                   shuffle["mrca"].loc[(shuffle["code"].str.contains("complex_core"))])
ax2.set_xlabel("mrca\nMWU = %s, pval = %s" %(mwu_stat, pval))

ax2.set_ylabel("Fold Change")

#####SIMPLE_CORE#####
simple = ["amber", "greyish"]
simplep = sns.xkcd_palette(simple)
sns.barplot(x = "taxon2", y = "fold_change",            data = pivot2.loc[pivot2["code"].str.contains("simple")].sort_values(by=["mrca_2", "code"]),              hue = "code", ax = ax3,  palette = simplep)
ax3.set_xticklabels(ax3.get_xticklabels(), rotation=90, horizontalalignment = "left", fontsize = 20)
mwu_stat, pval = stats.mannwhitneyu(final_merge["mrca"].loc[(final_merge["code"].str.contains("simple"))],                   shuffle["mrca"].loc[(shuffle["code"].str.contains("simple"))])
ax3.set_xlabel("mrca\nMWU = %s, pval = %s" %(mwu_stat, pval))

ax1.legend().remove()
ax1.set_xlabel("")
ax1.set_title("")
ax2.legend().remove()
ax2.set_xlabel("")
ax2.set_title("")
ax3.set_ylabel("Fold Change")
ax3.legend().remove()
ax3.set_xlabel("")
ax3.set_title("")

plt.savefig("%sFANTOM_v_shuffle_mrca_arch_fold_change_tx2.pdf"% (RE), bbox_inches = "tight")


# In[133]:


tx2.loc[(tx2["dataset"].str.contains("FANTOM")) &( tx2.code.str.contains("complex")| tx2.code.str.contains("simple"))].sort_values(by=["mrca_2", "code"] )


# In[134]:


fig, ax = plt.subplots(figsize=(8,8))

fig2e_plot = tx2.loc[(tx2["dataset"]                      .str.contains("FANTOM"))                      &( tx2.code.str.contains("complex")                       | tx2.code.str.contains("simple"))]                    .sort_values(by=["mrca_2", "code"] )
fig2e_plot= fig2e_plot.loc[fig2e_plot.taxon2 != "Homo"]
#####COMPLEX ENH#####
order = ["simple", "complex_core"]
arch_colors = [ "amber", "faded green"]
arch_palette = sns.xkcd_palette(arch_colors)
sns.set("talk")
sns.barplot(x = "taxon2", y = "frequency",
              data = fig2e_plot, hue = "code",palette = arch_palette, hue_order = order)
ax.set_xticklabels(ax.get_xticklabels(),
                    rotation=90, horizontalalignment = "left", fontsize = 20)

mwu_stat, pval = stats.mannwhitneyu(final_merge["mrca_2"].loc[(final_merge["code"].str.contains("complex_core"))],                   final_merge["mrca_2"].loc[(final_merge["code"].str.contains("simple"))])
print("complexenh v simple", mwu_stat, pval)

#ax.set_xlabel("mrca\nMWU = %s, pval = %s" %(mwu_stat, pval))
ax.set_title("Age Distribution")
ax.set_ylabel("Frequency")
ax.legend().remove()
ax.set_ylim(-0.01, 0.75)

ax.legend().remove()
ax.set_xlabel("")
ax.set_title("")


plt.savefig("%sfig2e-MRCA_fantom_complex_v_simple.pdf"% (RE), bbox_inches = "tight")


# In[135]:


arch_colors = [ "windows blue","greyish", "amber", "yellowish", "dusty purple",  "lavender"]
arch_palette = sns.xkcd_palette(arch_colors)
fig, (ax3, ax2) = plt.subplots(ncols = 2, figsize=(16,8))



#####COMPLEX ENH#####
arch_colors = ["faded green",  "greyish"]
arch_palette = sns.xkcd_palette(arch_colors)
sns.barplot(x = "taxon2", y = "frequency",
              data = tx2.loc[tx2["code"].str.contains("complex_core")].sort_values(by=["mrca_2", "code"]),\
              hue = "code2", ax = ax2,  palette = arch_palette)
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=90, horizontalalignment = "left", fontsize = 20)
mwu_stat, pval = stats.mannwhitneyu(final_merge["mrca_2"].loc[(final_merge["code"].str.contains("complex_core"))],                   shuffle["mrca_2"].loc[(shuffle["code"].str.contains("complex_core"))])
print("complexenh", mwu_stat, pval)

#ax2.set_xlabel("mrca\nMWU = %s, pval = %s" %(mwu_stat, pval))
ax2.set_title("Age Distribution")
ax2.set_ylabel("Frequency")
ax2.legend().remove()
ax2.set_ylim(-0.01, 0.75)

#####SIMPLE_CORE#####
arch_colors = ["amber", "greyish"]
arch_palette = sns.xkcd_palette(arch_colors)
sns.barplot(x = "taxon2", y = "frequency",
              data = tx2.loc[tx2["code"].str.contains("simple")].sort_values(by=["mrca_2", "code"]),\
              hue = "code2", ax = ax3,  palette = arch_palette)

ax3.set_xticklabels(ax3.get_xticklabels(), rotation=90, horizontalalignment = "left", fontsize = 20)

mwu_stat, pval = stats.mannwhitneyu(final_merge["mrca_2"].loc[(final_merge["code"].str.contains("simple"))],                   shuffle["mrca_2"].loc[(shuffle["code"].str.contains("simple"))])
print("simple", mwu_stat, pval)

#ax3.set_xlabel("mrca\nMWU = %s, pval = %s" %(mwu_stat, pval))
ax3.set_title("Age Distribution")
ax3.set_ylabel("Frequency")

"""
ax1.legend().remove()
ax1.set_xlabel("")
ax1.set_title("")
"""
ax2.legend().remove()
ax2.set_xlabel("")
ax2.set_title("")

ax3.legend().remove()
ax3.set_xlabel("")
ax3.set_title("")

plt.savefig("%sfigS2-MRCA_fantom_v_shuffle_core_derived_grouped_enh.pdf"% (RE), bbox_inches = "tight")


# In[136]:


pivot2 = pd.pivot_table(tx2, values='frequency',                           index=['mrca_2', 'code', 'taxon2'], columns=['dataset']).reset_index()
pivot2["fold_change"] = np.log2(pivot2["FANTOM"].divide(pivot2["shuffle"]))
pivot2=pivot2.fillna(0)

fig, (ax3, ax2) = plt.subplots(ncols = 2, figsize=(16,8))
ANNOTATIONS = 0

#####COMPLEX_CORE#####
complexcore = ["faded green",  "greyish"]
complexcorep = sns.xkcd_palette(complexcore)

sns.barplot(x = "taxon2", y = "fold_change",            data = pivot2.loc[pivot2["code"].str.contains("complex_core")].sort_values(by=["mrca_2", "code"]),              hue = "code", ax = ax2,  palette = complexcorep)
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=90, horizontalalignment = "left", fontsize = 20)
ax2.set_ylabel("log2(Fold Change)")
ax2.legend().remove()
ax2.set_ylim(-2.8, 2.1)
ax2.set_xlabel("")
ax2.set_title("")

mwu_stat, pval = stats.mannwhitneyu(final_merge["mrca"].loc[(final_merge["code"].str.contains("complex_core"))],                   shuffle["mrca"].loc[(shuffle["code"].str.contains("complex_core"))])

print("complexenh", mwu_stat, pval)

#####ANNOTATIONS#####
if ANNOTATIONS ==1:
    for p in ax2.patches:
        x = p.get_x().round(2)
        y = (p.get_height().round(2))
        ax2.annotate('{:.1f}'.format(p.get_height()), (x, y), color = "black", fontsize = 20)

#####SIMPLE_CORE#####
simple = ["amber", "greyish"]
simplep = sns.xkcd_palette(simple)
sns.barplot(x = "taxon2", y = "fold_change",            data = pivot2.loc[pivot2["code"].str.contains("simple")].sort_values(by=["mrca_2", "code"]),              hue = "code", ax = ax3,  palette = simplep)


ax3.set_ylabel("log2(Fold Change)")
ax3.set_ylim(-2.8, 2.1)
ax3.legend().remove()
ax3.set_xlabel("")
ax3.set_title("")
ax3.set_xticklabels(ax3.get_xticklabels(), rotation=90, horizontalalignment = "left", fontsize = 20)

mwu_stat, pval = stats.mannwhitneyu(final_merge["mrca"].loc[(final_merge["code"].str.contains("simple"))],                   shuffle["mrca"].loc[(shuffle["code"].str.contains("simple"))])
print("simple", mwu_stat, pval)

#####ANNOTATIONS#####
if ANNOTATIONS ==1:
    for p in ax3.patches:
        x = p.get_x().round(2)
        y = (p.get_height().round(2))
        ax3.annotate('{:.1f}'.format(p.get_height()), (x, y), color = "black", fontsize = 20)

plt.savefig("%sfigS2e-FANTOM_v_shuffle_mrca_arch_fold_change_tx2_enh.pdf"% (RE), bbox_inches = "tight")
