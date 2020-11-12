#!/usr/bin/env python
# coding: utf-8

# In[ ]:


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
RE ="/dors/capra_lab/projects/enhancer_ages/roadmap_encode/results/for_publication/non-genic/age_breaks/"


arch_colors = ["amber", "dusty purple", "windows blue","greyish"]
arch_palette = sns.xkcd_palette(arch_colors)
sns.palplot(arch_palette)

colors = [ "amber", "faded green",]
palette = sns.xkcd_palette(colors)
sns.palplot(palette)

cs = ["faded green", "greyish"]
cs_pal = sns.xkcd_palette(cs)
sns.palplot(cs_pal)

#%% Files


path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/"

fs = glob.glob("%sbreaks/ROADMAP_E*_enh_and_shuf_age_arch_summary_matrix.bed" %path)

sid_list = []

for f in fs:
    sid = ((f.split("/")[-1]).split("_")[1]).split(".")[0]
    if sid not in sid_list:
        sid_list.append(sid)



# In[2]:



#%% other summary files

# age and taxon file
syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t') # read the file
syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2", "mya", "mya2"]] # whittle down the df

#%% sample description file
desc_df = pd.read_csv("/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/roadmap_hg19_sample_id_desc.csv", header =None)
desc_df.columns = ["sid", "desc"]
desc_df.head()

#%% relative simple architecture map file
pdf = pd.read_csv("/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/breaks/all_ROADMAP_breaks_percentiles.bed", sep = '\t')


# In[10]:




#FUNCTIONS#

def get_relative_simple(sid, pdf):


    test = pdf.loc[(pdf.sid2 == sid)]

    relative_simple = test.loc[test.pct < 0.5, "seg_index"].max()

    percentile = test.loc[test.pct < 0.5, "pct"].max()

    print(sid, relative_simple, percentile)

    return relative_simple, percentile


def get_balanced_shuffles(df):

    chr_count = df[["shuf_id", "dataset", "chr_enh"]].drop_duplicates()
    chr_count = chr_count.groupby(["shuf_id", "dataset"])["chr_enh"].count().reset_index()
    drop_counts = chr_count.loc[(chr_count.chr_enh<21) & (chr_count.dataset =="shuf"), "shuf_id"].to_list()

    # keep only complete shuffle datasets by keeping datasets w/ line counts within stdev of mean.

    count = df.groupby(["shuf_id", "dataset"])["enh_id"].count().reset_index() # count lines per shuffle

    count["z"] = (count["enh_id"].mean() - count["enh_id"]).divide(count["enh_id"].std()) # get z-score

    drop_datasets = count.loc[(count.z.abs() > 1) &(count.dataset =="shuf"), "shuf_id"].to_list() # identify datasets within 1 sd

    drops = drop_counts + drop_datasets
    print(drops)

    df = df.loc[~df.shuf_id.isin(drops)]# keep onlly those datasets

    return df


def get_combined_df(sid, path, relative_simple):

    f = "%sbreaks/ROADMAP_%s_enh_and_shuf_age_arch_summary_matrix.bed" % (path, sid)

    df = pd.read_csv(f, sep = '\t', header = None, error_bad_lines = False)

    # name columns
    if len(list(df)) == 18:
        df.columns = ["chr_enh", "start_enh", "end_enh", "dataset_id", "enh_id",
                       "shuf_id", "core_remodeling", "arch", "seg_index", "mrca",
                       "enh_len", "taxon", "mrca_2", "taxon2","mya", "mya2",
                       "seg_density", "dataset" ]

    elif len(list(df)) == 17:
        df.columns = ["chr_enh", "start_enh", "end_enh",
                       "shuf_id","enh_id", "core_remodeling", "arch", "seg_index", "mrca",
                       "enh_len", "taxon", "mrca_2", "taxon2","mya", "mya2",
                       "seg_density", "dataset" ]
    sex_chrom = ["chrX", "chrY", "chrM"]

    df = df.loc[~df.chr_enh.isin(sex_chrom)]
    drop_enh_id = df.loc[~(df.core_remodeling ==1) &(df.mrca_2==0.0), "enh_id"].to_list()
    print(len(drop_enh_id))
    df.loc[~df.enh_id.isin(drop_enh_id)]


    df.loc[df.dataset.str.contains("shuf"), "dataset"] = "shuf"
    df.loc[~df.dataset.str.contains("shuf"), "dataset"] = "roadmap"

    # clean up
    df.enh_id = df.chr_enh + ":" + df.start_enh.map(str) + "-" + df.end_enh.map(str)
    df.loc[df.shuf_id.astype(str).str.contains("chr"), "shuf_id"] = sid + "age_breaks"

    df_bal = get_balanced_shuffles(df)

    df_bal.mrca_2 = df_bal.mrca_2.round(3)

    df_bal["relative_arch"] = "rel_simple"
    df_bal.loc[df_bal.seg_index.astype(float) >= float(relative_simple), "relative_arch"] = "rel_complex"
    df_bal.loc[df_bal.relative_arch == "rel_simple", "core_remodeling"] = 0

    shuffle = df_bal.loc[df_bal.dataset.str.contains("shuf")]


    df = df_bal.loc[~df_bal.dataset.str.contains("shuf")]


    return shuffle, df


def get_simple_OR(df, shuffle, sid):

    enh_totals = len(df) # total enhancers

    shuf_totals = len(shuffle) # total shuffles

    a = df.loc[df.core_remodeling ==0, "enh_id"].count() # num simple enhancers
    b = enh_totals - a # num complex enhancers

    c = shuffle.loc[shuffle.core_remodeling ==0, "enh_id"].count() # num simple shuffle
    d = shuf_totals - c # num complex shuffle

    obs = [[a,b], [c,d]]
    OR, P = stats.fisher_exact(obs)
    table = sm.stats.Table2x2(obs) # get confidence interval
    odds_ci = table.oddsratio_confint()
    ORdf = pd.DataFrame({ "a":[a], "b":[b], "c":[c], "d":[d],
                         "OR":[OR], "P":[P], "ci_lower" :[odds_ci[0]],
                        "ci_upper" :[odds_ci[1]]})

    print(sid, obs, OR, P, round((a/enh_totals), 2), round((c/shuf_totals), 2))

    ORdf["yerr"] = ORdf.ci_upper-ORdf.ci_lower
    ORdf['log'] = np.log2(ORdf.OR)
    ORdf['sid'] = sid

    return ORdf


def get_enh_len(df, sid):

    lens = df[["enh_id", "core_remodeling", "enh_len", "mrca_2", "taxon2"]].drop_duplicates()

    lensArch = lens.groupby(["core_remodeling","mrca_2", "taxon2"])["enh_len"].mean().reset_index()

    lensArch["sid"] = sid

    return lensArch


def get_freq_arch(df, sid):

    freq = df[["enh_id", "core_remodeling", "mrca_2", "taxon2"]].drop_duplicates()

    Narch = freq.groupby("core_remodeling")["enh_id"].count().reset_index()
    Narch.columns = ["core_remodeling", "total_arch"]

    freq = freq.groupby(["core_remodeling", "mrca_2", "taxon2"])["enh_id"].count().reset_index()
    freq.columns = ["core_remodeling", "mrca_2", "taxon2", "arch_mrca_count"]

    freq = pd.merge(freq, Narch, how = "left", on = "core_remodeling")

    freq["arch_freq"] = freq.arch_mrca_count.divide(freq.total_arch) # get the frequency of architecture per age

    freq["sid"] = sid

    return freq


def get_mean_seg(df, shuffle, sid):

    mean_enh, sd_enh = df.seg_index.mean(), df.seg_index.std()

    mean_shuf, sd_shuf = shuffle.seg_index.mean(), shuffle.seg_index.std()

    ratio = round(mean_enh/mean_shuf, 2) # the ratio of mean enh segs: mean shuf segs


    segdf = pd.DataFrame({"sid" : sid,
                         "mean_enh_seg":[mean_enh],
                         "sd_enh_seg": [sd_enh],
                         "mean_shuf_seg": [mean_shuf],
                         "sd_shuf_seg": [sd_shuf],
                         "fold_change":[ratio]})


    return segdf


def get_proportion_simple(enhdf, relative_simple, shuffle, sid):

    enh_percentile = enhdf.loc[enhdf.seg_index<= relative_simple].count()

    shuf_percentile = shuffle.loc[shuffle.seg_index<= relative_simple].count()

    proportiondf = pd.DataFrame({"sid": [sid],
                                "enh_proportion": [enh_percentile],
                                "shuffle_proportion":[shuf_percentile]})

    return proportiondf


def get_data(sid, path, pdf):


    relative_simple, percentile = get_relative_simple(sid, pdf) # get relative simple

    shuffle, enhdf = get_combined_df(sid, path, relative_simple) # get shuffles

    ORdf = get_simple_OR(enhdf, shuffle, sid)

    enhLens = get_enh_len(enhdf, sid)

    enhFreq = get_freq_arch(enhdf, sid)

    enhSeg = get_mean_seg(enhdf, shuffle, sid)

    enhProp = get_proportion_simple(enhdf, relative_simple, shuffle, sid)

    return ORdf, enhLens, enhFreq, enhSeg, enhProp


# In[ ]:


# combined

OR ={}
lens = {}
freq = {}
seg = {}
prop = {}


for sid in sid_list:
    ORdf, enhLens, enhFreq, enhSeg, enhProp = get_data(sid, path, pdf)

    OR[sid] = ORdf
    lens[sid] = enhLens
    freq[sid] = enhFreq
    seg[sid] = enhSeg
    prop[sid] = enhProp

ORdf = pd.concat(OR.values())
lensdf = pd.concat(lens.values())
freqdf = pd.concat(freq.values())
segdf = pd.concat(seg.values())
propdf = pd.concat(prop.values())


# format the OR dataframe

ORdf = pd.merge(ORdf, desc_df, how = "left", on = "sid") # add tissue descriptions

ORdf["sid2"] = ORdf.sid + "-" + ORdf.desc # make variable for dataset + tissue descriptions

ORdf["sig"] = ""
ORdf.loc[ORdf.P<0.05, "sig"] = "*"
ORdf["simple"]= "simple"
ORdf.loc[ORdf.OR<1, "simple"] = "complex"
ORdf.head()


# In[13]:


sig = ORdf.loc[ORdf.P<0.05]
print("n simple enriched", len(ORdf.loc[ORdf.log >0]), "\nn total", len(ORdf))
print("n sig simple enriched", len(sig.loc[sig.log >0]), "\nn sig total", len(sig))


# In[31]:


sns.set("poster")
sns.set_style("white")
fig, ax = plt.subplots(figsize=(8, 40))
sns.barplot( y= "desc", x ="log",
            data = ORdf.sort_values(by="log", ascending = False),
           linewidth=2.5, facecolor=(1, 1, 1, 0),
             edgecolor=".2", xerr=(ORdf["ci_upper"] - ORdf["ci_lower"]),
             ci="sd", errwidth = 5)

ax.set(xlabel = "Fold enrichment (log2-scaled)", ylabel = '')

# annotate graph with significant enrichment #

rects = ax.patches

labels = ORdf.sort_values(by="log", ascending = False)["sig"]

val = 0
for rect, label in zip(rects, labels):
    h = rect.get_height() + val
    w = rect.get_width()
    if w>0:
        w +=0.2
    else:

        w -=0.2

    ax.text(w,h, label,            ha='center', va='bottom')
    val +=1

ax.axvline(0, color = 'k') # add horizontal line
from matplotlib.ticker import MultipleLocator
import matplotlib.ticker as ticker
ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(round(2**x, 1)))
ax.xaxis.set_major_formatter(ticks)
ax.xaxis.set_major_locator(MultipleLocator(0.5))
plt.savefig("%sFigS2.5Arelative_simple_ROADMAP_fold_enrichment_all.pdf" % RE, bbox_inches = "tight")


# # per tissue sig enrichments box plot

# In[15]:


fig, ax = plt.subplots(figsize=(6, 6))
sns.boxplot( y= "OR", x ="simple",data = ORdf, linewidth=2.5, palette = palette)
sns.swarmplot( y= "OR", x ="simple",data = sig, linewidth=2.5, palette = palette)
ax.set(ylabel = "Fold enrichment\n(log2-scaled)", xlabel = "")
simplen = len(sig.loc[sig.simple == "simple"])
complen = len(sig.loc[sig.simple != "simple"])
labels = ["simple n=%s"%simplen, "complex n=%s" %complen]
ax.set_xticklabels(labels)
ax.axhline(1, color = "grey",)# linestyle = "--")

plt.savefig("%sFigS2.5Arelative_simple_ROADMAP_fold_enrichment_sig.pdf" % RE, bbox_inches = "tight")


# # enh length bar plot

# In[16]:


# drop complex enhancers of human age. Messy data.

indexNames  = lensdf[ (lensdf.mrca_2 ==0) & (lensdf.core_remodeling ==1) ].index
lensdf.drop(indexNames , inplace=True)


# In[17]:


fig, ax = plt.subplots(figsize=(8, 8))
sns.barplot( x= "taxon2", y ="enh_len",
            data = lensdf.sort_values(by="mrca_2"),
            hue = "core_remodeling",
           palette = palette)
ax.set(ylabel = "Length(bp)", xlabel = "")
ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
ax.legend().remove()
plt.savefig("%sFigS2.5Brelative_simple_ROADMAP_enh_lengths.pdf" % RE, bbox_inches = "tight")


# # age freq bar plot

# In[18]:


fig, ax = plt.subplots(figsize=(8, 8))
sns.barplot(y ="arch_freq", x= "taxon2",
            data = freqdf.sort_values(by="mrca_2"),
            hue = "core_remodeling",
           palette = palette)
ax.set(ylabel = "Total Architecture (%)", xlabel = "")
ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
ax.legend().remove()
plt.savefig("%sFigS2.5Drelative_simple_ROADMAP_age_freq.pdf" % RE, bbox_inches = "tight")


# # mean seg ratio bar plot

# In[34]:


len(segdf.loc[segdf.fold_change<=1]), len(segdf.loc[segdf.fold_change>1])


# In[32]:


fig, ax = plt.subplots(figsize=(8, 30))
sns.barplot(x ="fold_change", y= "sid",
            data = segdf.sort_values(by="fold_change"),
            linewidth=2.5, facecolor=(1, 1, 1, 0),
             edgecolor=".2",)

ax.set(ylabel = "", xlabel = "simple enh:shuf")
#ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
ax.axvline(1, color = 'k')
plt.savefig("%sFigS2.xrelative_simple_ROADMAP_seg.pdf" % RE, bbox_inches = "tight")


# In[41]:


def get_proportion_simple(enhdf, relative_simple, shuffle, sid):

    enh_percentile = len(enhdf.loc[enhdf.seg_index<= relative_simple])

    shuf_percentile = len(shuffle.loc[shuffle.seg_index<= relative_simple])

    proportiondf = pd.DataFrame({"sid": [sid],
                                "enh_proportion": [enh_percentile],
                                "shuffle_proportion":[shuf_percentile]})

    return proportiondf

def get_data(sid, path, pdf):


    relative_simple, percentile = get_relative_simple(sid, pdf) # get relative simple

    shuffle, enhdf = get_combined_df(sid, path, relative_simple) # get shuffles

    propdf = get_proportion_simple(enhdf, relative_simple, shuffle, sid)

    return propdf


# In[42]:


prop ={}
for sid in sid_list:
    enhProp = get_data(sid, path, pdf)
    prop[sid] = enhProp
propdf = pd.concat(prop.values())


# In[40]:


propdf.head()


# In[38]:


sns.distplot(propdf.enh_proportion)
sns.distplot(propdf.shuf_proportion)


# ### dev stuff

# In[ ]:
"""

sid = "E115"
relative_simple, percentile = get_relative_simple(sid, pdf)
f = "%sbreaks/to_delete/no-exon_%s.bed" % (path, sid)

df = pd.read_csv(f, sep = '\t', header = None)

# name columns
if len(list(df)) == 18:
    df.columns = ["chr_enh", "start_enh", "end_enh", "enh_id", "peak_counts",
                   "shuf_id", "core_remodeling", "arch", "seg_index", "mrca",
                   "enh_len", "taxon", "mrca_2", "taxon2","mya", "mya2",
                   "seg_density", "dataset" ]

elif len(list(df)) == 17:
    df.columns = ["chr_enh", "start_enh", "end_enh",
                   "shuf_id","enh_id", "core_remodeling", "arch", "seg_index", "mrca",
                   "enh_len", "taxon", "mrca_2", "taxon2","mya", "mya2",
                   "seg_density", "dataset" ]
sex_chrom = ["chrX", "chrY", "chrM"]

df = df.loc[~df.chr_enh.isin(sex_chrom)]
drop_enh_id = df.loc[~(df.core_remodeling ==1) &(df.mrca_2==0.0), "enh_id"].to_list()
print(len(drop_enh_id))
df.loc[~df.enh_id.isin(drop_enh_id)]
df.loc[df.dataset.str.contains("shuf"), "dataset"] = "shuf"
df.loc[~df.dataset.str.contains("shuf"), "dataset"] = "roadmap"

# clean up
df.enh_id = df.chr_enh + ":" + df.start_enh.map(str) + "-" + df.end_enh.map(str)
df.loc[df.shuf_id.astype(str).str.contains("chr"), "shuf_id"] = sid + "age_breaks"

df_bal = get_balanced_shuffles(df)

df_bal.mrca_2 = df_bal.mrca_2.round(3)

df_bal["relative_arch"] = "rel_simple"
df_bal.loc[df_bal.seg_index.astype(float) >= float(relative_simple), "relative_arch"] = "rel_complex"
df_bal.loc[df_bal.relative_arch == "rel_simple", "core_remodeling"] = 0

shuffle = df_bal.loc[df_bal.dataset.str.contains("shuf")]


df = df_bal.loc[~df_bal.dataset.str.contains("shuf")]



df_bal.shuf_id.unique()

enh_totals = len(df) # total enhancers

shuf_totals = len(shuffle) # total shuffles

a = df.loc[df.core_remodeling ==0, "enh_id"].count() # num simple enhancers
b = enh_totals - a # num complex enhancers

c = shuffle.loc[shuffle.core_remodeling ==0, "enh_id"].count() # num simple shuffle
d = shuf_totals - c # num complex shuffle

obs = [[a,b], [c,d]]
OR, P = stats.fisher_exact(obs)
table = sm.stats.Table2x2(obs) # get confidence interval
odds_ci = table.oddsratio_confint()
ORdf = pd.DataFrame({ "a":[a], "b":[b], "c":[c], "d":[d],
                     "OR":[OR], "P":[P], "ci_lower" :[odds_ci[0]],
                    "ci_upper" :[odds_ci[1]]})

print(sid, obs, OR, P, round((a/enh_totals), 2), round((c/shuf_totals), 2))

ORdf["yerr"] = ORdf.ci_upper-ORdf.ci_lower
ORdf['log'] = np.log2(ORdf.OR)
ORdf['sid'] = sid


# In[ ]:


sids = ["E044", "E012"]
inpath = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/non-genic/breaks/"
for sid in sids:

    f = "%sno-exon_%s_parallel_breaks_enh_age_arch_summary_matrix.bed" % (inpath, sid)
    df = pd.read_csv(f, sep ='\t', header = None, usecols = [0,1,2,5,3,6,7,8,9,10,11,12,13,14,15,16,17])
    outf = "%s/to_delete/%s.bed" % (inpath, sid)
    df.to_csv(outf, sep = '\t', header = False, index = False)
    temp = "%sto_delete/temp_%s.bed" % (inpath, sid)
    final_f = "%sto_delete/no-exon_ROADMAP_%s_enh_and_shuf_age_arch_summary_matrix.bed"% (inpath, sid)
    cmd = "cat %s %s > %s && mv %s %s" %(inpath, final_f, temp, temp, final_f)
    print(cmd)
"""
