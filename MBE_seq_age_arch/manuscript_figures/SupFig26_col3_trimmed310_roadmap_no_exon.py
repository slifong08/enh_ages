#!/usr/bin/env python
# coding: utf-8

# In[20]:


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

colors = ["amber", "faded green",]
palette = sns.xkcd_palette(colors)
sns.palplot(palette)

cs = ["faded green", "greyish"]
cs_pal = sns.xkcd_palette(cs)
sns.palplot(cs_pal)

#%% Files


path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19/download/h3k27ac_plus_h3k4me3_minus_peaks/"

fs = glob.glob("%s*/trimmed/no-exon*.bed" %path)


shufs = glob.glob("%s*/trimmed/shuffle/breaks/no-exon_*.bed" %path)


sid_list = []

for f in fs:
    sid = ((f.split("/")[-1]).split("_")[1]).split(".")[0]
    if sid not in sid_list:
        sid_list.append(sid)



# In[21]:


len(fs), len(shufs)


# In[24]:


outpath


# In[22]:


outpath = "%sstats/" % path


#%% other summary files

# age and taxon file
syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t') # read the file
syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2", "mya", "mya2"]] # whittle down the df

#%% sample description file
desc_df = pd.read_csv("/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19/roadmap_hg19_sample_id_desc.csv", header =None)
desc_df.columns = ["sid", "desc"]
desc_df.head()


# In[4]:


sid_fs = {}
missing = []
fss = fs + shufs
for f in fss:
    fid = (f.split("/")[-1]).split(".")[0]

    print(fid)
    if "shuf" in fid:
        sid = (fid.split("_")[1]).split("-")[1]
        print("shuf", sid)
    else:
        sid = (fid.split("-")[-1])
       #print(sid)

    if sid not in sid_fs:
        sid_fs[sid] = [f]


    else:
        f_list = sid_fs[sid]
        if f not in f_list:
            f_list.append(f)
            sid_fs[sid] = f_list
for sid, fs_list in sid_fs.items():
    if len(fs_list)>1:
        print(sid)
    else:
        print("where are the shuffles?", sid)
        missing.append(sid)


# In[5]:


print(missing)


# In[9]:


#FUNCTIONS#


def combine_and_format_sid(sid, fs_list): # put the enhancer and shuffle datasets together.

    sample_df = {}

    for f in fs_list:
        fid = (f.split("/")[-1]).split(".")[0]

        cols = pd.read_csv(f, sep = '\t', header = None, nrows = 1)
        print(fid, len(list(cols)))
        if len(list(cols)) ==18:
            df = pd.read_csv(f, sep = '\t', low_memory=False,
             error_bad_lines = False, usecols=[0,1,2,3,6,8,9])
            df.columns = ['chr_enh','start_enh','end_enh', 'enh_id','core_remodeling', 'seg_index',
                   'mrca']
            df["id"] = fid

        elif "shuf" in fid:
            df = pd.read_csv(f, sep = '\t', header = None,low_memory=False,
             error_bad_lines = False, usecols=[0,1,2,3,4,5,6,8])
            df.columns = ['chr_enh','start_enh','end_enh', "id", 'enh_id', 'seg_index',
                  'core_remodeling', 'mrca']

            df = df.loc[df.mrca != "max_age"]
        else:
            df = pd.read_csv(f, sep = '\t', low_memory=False,
             error_bad_lines = False, usecols=[0,1,2,3,4,6,7])
            df.columns = ['chr_enh','start_enh','end_enh', 'enh_id','core_remodeling', 'seg_index',
                   'mrca']
            df["id"] = fid


        df.mrca = df.mrca.astype(float).round(3)
        df = pd.merge(df, syn_gen_bkgd[["mrca", "mrca_2", "taxon2"]], how = "left", on = "mrca")

        df["enh_len"] = df["end_enh"] - df["start_enh"]
        df["datatype"] = fid

        df = df[['chr_enh','start_enh','end_enh', 'enh_id', "enh_len", 'core_remodeling', 'seg_index',
               'mrca_2', "taxon2", "id", "datatype"]].drop_duplicates()


        if sid not in sample_df.keys():
            sample_df[sid] = df # add the dataframe
        else:
            df_ = sample_df[sid]
            newdf = pd.concat([df_, df]) # concatenate the old and new datasets together.
            sample_df[sid] = newdf

    return sample_df[sid]


def cleanup_df(df, sid,  path):


    sex_chrom = ["chrX", "chrY", "chrM"]

    df = df.loc[~df.chr_enh.isin(sex_chrom)]
    drop_enh_id = df.loc[~(df.core_remodeling ==1) &(df.mrca_2==0.0), "enh_id"].to_list()
    print(len(drop_enh_id))
    df.loc[~df.enh_id.isin(drop_enh_id)]


    df.loc[df.datatype.str.contains("shuf"), "dataset"] = "shuf"
    df.loc[~df.datatype.str.contains("shuf"), "dataset"] = "roadmap"

    # clean up
    df.enh_id = df.chr_enh + ":" + df.start_enh.map(str) + "-" + df.end_enh.map(str)

    shuffle = df.loc[df.dataset.str.contains("shuf")]

    enhdf = df.loc[~df.dataset.str.contains("shuf")]


    return shuffle, enhdf


def get_simple_OR(enhdf, shuffle, sid):

    enh_totals = len(enhdf) # total enhancers

    shuf_totals = len(shuffle) # total shuffles

    a = enhdf.loc[enhdf.core_remodeling ==0, "enh_id"].count() # num simple enhancers
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

def get_seg_OR(enhdf, shuffle, sid):

    seg_dict = {}

    for seg_index in enhdf.seg_index.unique():

        enh_totals = len(enhdf) # total enhancers

        shuf_totals = len(shuffle) # total shuffles

        a = enhdf.loc[enhdf.seg_index == seg_index, "enh_id"].count() # num simple enhancers
        b = enh_totals - a # num complex enhancers

        c = shuffle.loc[shuffle.seg_index == seg_index, "enh_id"].count() # num simple shuffle
        d = shuf_totals - c # num complex shuffle

        obs = [[a,b], [c,d]]
        OR, P = stats.fisher_exact(obs)
        table = sm.stats.Table2x2(obs) # get confidence interval
        odds_ci = table.oddsratio_confint()
        ORdf = pd.DataFrame({ "a":[a], "b":[b], "c":[c], "d":[d],
                             "OR":[OR], "P":[P], "ci_lower" :[odds_ci[0]],
                            "ci_upper" :[odds_ci[1]]})

        print(sid, obs, OR, P, round((a/enh_totals), 2), round((c/shuf_totals), 2), seg_index)

        ORdf["yerr"] = ORdf.ci_upper-ORdf.ci_lower
        ORdf['log'] = np.log2(ORdf.OR)
        ORdf['sid'] = sid
        ORdf["seg_index"] = seg_index
        seg_dict[seg_index] = ORdf

    segdf = pd.concat(seg_dict.values())
    return segdf


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


def get_proportion_simple(enhdf, shuffle, sid):

    enh_percentile = len(enhdf.loc[enhdf.seg_index<= 1])

    shuf_percentile = len(shuffle.loc[shuffle.seg_index<= 1].count())

    proportiondf = pd.DataFrame({"sid": [sid],
                                "enh_proportion": [enh_percentile],
                                "shuffle_proportion":[shuf_percentile]})

    return proportiondf

def format_df(sid, fs_list, path):

    sid_df = combine_and_format_sid(sid, fs_list) # combine the datasets only to split them. This is a patch.

    shuffle, enhdf = cleanup_df(sid_df, sid, path) # get shuffles

    return sid_df, shuffle, enhdf

def get_data(sid, shuffle, enhdf):


    ORdf = get_simple_OR(enhdf, shuffle, sid)

    enhLens = get_enh_len(enhdf, sid)

    enhFreq = get_freq_arch(enhdf, sid)

    enhSeg = get_mean_seg(enhdf, shuffle, sid)

    enhSegOR = get_seg_OR(enhdf, shuffle, sid)

    enhProp = get_proportion_simple(enhdf, shuffle, sid)

    return ORdf, enhLens, enhFreq, enhSeg, enhProp, enhSegOR


# In[10]:


# combined

OR ={}
lens = {}
freq = {}
seg = {}
prop = {}
segOR = {}


# In[11]:


for sid, fs_list in sid_fs.items():
    if sid not in OR.keys() and sid not in missing:

        sid_df, shuffle, enhdf = format_df(sid, fs_list, path)
        ORdf, enhLens, enhFreq, enhSeg, enhProp, enhSegOR = get_data(sid, shuffle, enhdf)

        OR[sid] = ORdf
        lens[sid] = enhLens
        freq[sid] = enhFreq
        seg[sid] = enhSeg
        prop[sid] = enhProp
        segOR[sid] = enhSegOR


ORdf = pd.concat(OR.values())
lensdf = pd.concat(lens.values())
freqdf = pd.concat(freq.values())
segdf = pd.concat(seg.values())
propdf = pd.concat(prop.values())
segORdf = pd.concat(segOR.values())


# In[12]:


# format the OR dataframe

ORdf = pd.merge(ORdf, desc_df, how = "left", on = "sid") # add tissue descriptions

ORdf["sid2"] = ORdf.sid + "-" + ORdf.desc # make variable for dataset + tissue descriptions

ORdf["sig"] = ""
ORdf.loc[ORdf.P<0.05, "sig"] = "*"
ORdf["simple"]= "simple"
ORdf.loc[ORdf.OR<1, "simple"] = "complex"
ORdf.head()


# In[13]:


fs_list


# In[14]:


shuffle.head()


# In[25]:


measures = {"OR": ORdf, "LEN": lensdf, "FREQ":freqdf, "SEG": segdf, "PROPORTION":propdf}

for m, d in measures.items():
    outf = "%s%s_98_datasets_simple_trimmed.txt"% (outpath, m)
    d.to_csv(outf, sep = '\t', header = True, index = False)


# In[26]:


sig = ORdf.loc[ORdf.P<0.05]
print("n simple enriched", len(ORdf.loc[ORdf.log >0]), "\nn total", len(ORdf))
print("n sig simple enriched", len(sig.loc[sig.log >0]), "\nn sig total", len(sig))


# In[27]:


ORdf["log2_ci_upper"],ORdf["log2_ci_lower"]  = np.log2(ORdf["ci_upper"]),  np.log2(ORdf["ci_lower"])
ORdf.loc[ORdf.desc.str.contains("CD3")]


# # per tissue sig enrichments box plot

# In[29]:


ORdf.head()


# In[30]:


sns.set("poster")
sns.set_style("white")
fig, ax = plt.subplots(figsize=(8, 40))
ORdf = ORdf.sort_values(by="log", ascending = False)
sns.barplot( y= "desc", x ="log",
            data = ORdf.sort_values(by="log", ascending = False),
           linewidth=2.5, facecolor=(1, 1, 1, 0),
             edgecolor=".2",
            xerr=abs(ORdf["log2_ci_upper"] - ORdf["log2_ci_lower"]),
             ci="sd", errwidth = 5)

ax.set(xlabel = "Fold enrichment (log2-scaled)",
       #xlim = (-4,2),
       ylabel = '')

# annotate graph with significant enrichment #

rects = ax.patches

labels = ORdf.sort_values(by="log", ascending = False)["sig"]

val = 0
for rect, label in zip(rects, labels):
    h = rect.get_height() + val
    w = rect.get_width()
    if w>0:
        w +=0.15
    else:

        w -=0.15

    ax.text(w,h, label,            ha='center', va='bottom')
    val +=1

ax.axvline(0, color = 'k') # add horizontal line
from matplotlib.ticker import MultipleLocator
import matplotlib.ticker as ticker
ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(round(2**x, 1)))
ax.xaxis.set_major_formatter(ticks)
ax.xaxis.set_major_locator(MultipleLocator(0.25))
plt.savefig("%sFigS2.7_ROADMAP_fold_enrichment_310_98_tissues.pdf" % RE, bbox_inches = "tight")


# In[31]:


order = ["simple", "complex"]
fig, ax = plt.subplots(figsize=(6, 6))
sns.boxplot( y= "OR", x ="simple", data = ORdf, linewidth=2.5, palette = palette, order =order)
sns.swarmplot( y= "OR", x ="simple",data = sig, linewidth=2.5, palette = palette, order =order)
ax.set(ylabel = "Fold enrichment\n(log2-scaled)", xlabel = "")
simplen = len(sig.loc[sig.simple == "simple"])
complen = len(sig.loc[sig.simple != "simple"])
labels = ["simple n=%s"%simplen, "complex n=%s" %complen]
ax.set_xticklabels(labels)
ax.axhline(1, color = "grey",)# linestyle = "--")

plt.savefig("%sFigS2.7_ROADMAP_fold_enrichment_310_sig_enrichment_64_tissues.pdf" % RE, bbox_inches = "tight")


# # enh length bar plot

# In[32]:


lensdf.head()


# In[33]:


# drop complex enhancers of human age. Messy data.

indexNames  = lensdf[ (lensdf.mrca_2 ==0) & (lensdf.core_remodeling ==1) ].index
lensdf.drop(indexNames , inplace=True)

fig, ax = plt.subplots(figsize=(8, 8))
sns.barplot( x= "taxon2", y ="enh_len",
            data = lensdf.sort_values(by="mrca_2"),
            hue = "core_remodeling",
           palette = palette)
ax.set(ylabel = "Length(bp)", xlabel = "")
ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
ax.legend().remove()
plt.savefig("%sFigS2.7_ROADMAP_fold_enrichment_310_len_98_tissues.pdf" % RE, bbox_inches = "tight")


# # age freq bar plot

# In[34]:


fig, ax = plt.subplots(figsize=(8, 8))
sns.barplot(y ="arch_freq", x= "taxon2",
            data = freqdf.sort_values(by="mrca_2"),
            hue = "core_remodeling",
           palette = palette)
ax.set(ylabel = "Total Architecture (%)", xlabel = "")
ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
ax.legend().remove()
plt.savefig("%sFigS2.7_ROADMAP_310_mrca_x_percent_98_tissues.pdf" % RE, bbox_inches = "tight")


# # mean seg ratio bar plot

# In[36]:


segORdf["log2"] = np.log2(segORdf.OR)
segORdf.head()
test = segORdf.loc[segORdf.seg_index<8]
fig, ax = plt.subplots(figsize = (6,6))
sns.barplot(x = "seg_index", y = "log2",
            data= test.sort_values(by = "seg_index"),
           linewidth=2.5, facecolor=(1, 1, 1, 0),
             edgecolor=".2")
            #hue = "core_remodeling",
           #palette = palette)
ax.set(xlabel = "No. age segments", ylabel = "No. age segments v. bkgd\n OR (log2-scaled)")
ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(round(2**x, 1)))
ax.yaxis.set_major_formatter(ticks)
ax.yaxis.set_major_locator(MultipleLocator(0.5))
ax.legend().remove()
#ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
#plt.savefig("%sFig2_age_seg_enrichment_98_tissues.pdf"%RE, bbox_inches = 'tight')


# In[ ]:
