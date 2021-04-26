#!/usr/bin/env python
# coding: utf-8

# In[1]:


import glob
import pandas
import matplotlib.pyplot as plt
import numpy as np
import os, sys
from scipy import stats
import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')

# sns colors
arch_colors = [ "amber", "dusty purple", "windows blue","greyish"]
arch_palette = sns.xkcd_palette(arch_colors)
sns.palplot(arch_palette)

# sns graphing preferences
sns.set(color_codes=True)
sns.set(font_scale=1.5)
sns.set_style("white")
sns.despine(bottom=True, left=True)

import datetime
last_run = datetime.datetime.now()
today = (datetime.date.today())
print("last run", datetime.datetime.now())
RE = "/dors/capra_lab/projects/enhancer_ages/linsight/results/"


# In[2]:


linsight_path = "/dors/capra_lab/projects/enhancer_ages/linsight/data/"
fantom_fs = glob.glob("%sFANTOM_*.bed" % linsight_path)

# enh only
fantom_fs = "/dors/capra_lab/projects/enhancer_ages/linsight/data/all_unique_fantom_erna_112_tissue_linsight.bed"


# In[3]:

BUILD = "hg19"
TRIM_LEN = "310"
FIG_ID = "3C"
FRAC = 0.5

def load_syn_gen_bkgd(build):

    F = f"/dors/capra_lab/projects/enhancer_ages/{build}_syn_taxon.bed"
    syngenbkgd = pd.read_csv(F, sep='\t')
    syngenbkgd[["mrca", "mrca_2"]] = syngenbkgd[["mrca", "mrca_2"]].round(3)

    return syngenbkgd


def assign_architecture(df):

    df["code"] = ""
    df.code.loc[(df.core_remodeling == 0)& (df.core == 1)] = "simple"
    df.code.loc[(df.core_remodeling == 1)& (df.core == 1)] = "complex_core"
    df.code.loc[(df.core_remodeling == 1)& (df.core == 0)] = "derived"

    df["arch"] = ""
    df.arch.loc[(df.core_remodeling == 0)] = "simple"
    df.arch.loc[(df.core_remodeling == 1)] = "complex"

    return df


def format_df(fantom_fs, syn_gen_bkgd):

    arch_id = (fantom_fs.split("/")[-1]).split("_")[1]
    print(arch_id)
    df = pandas.read_csv(fantom_fs, sep = '\t', header = None, low_memory=False)

    # rename columns
    df.columns = ['chr_syn','start_syn','end_syn', 'enh_id',
                  'chr_enh', 'start_enh','end_enh',
                  'seg_index', 'core_remodeling', 'core',
                  'mrca', 'code', 'syn_id',
                  "chr_lin",  "start_lin", "end_lin","linsight_score", "overlap"]

    # quantify lengths
    df['enh_len'] = df.end_enh - df.start_enh
    df['syn_len'] = df.end_syn - df.start_syn


    core_age = df.groupby("enh_id")['mrca'].max().reset_index()

    # assign architecture, sub architectures
    df = assign_architecture(df)

    # Format LINSIGHT information
    # exclude the loci that do not overlap a linsight score
    df = df.loc[df.linsight_score != "."]

    df.linsight_score = df.linsight_score.astype(float) # turn linsight scores into floats

    df["linsight_id"] = df.chr_lin + ":" + df.start_lin.map(str) +"-"+ df.end_lin.map(str)

    # make a dataframe only based on linsight scores and enhancer architectures

    base_df = df[["chr_lin", "start_lin", "end_lin", "linsight_score", "code", "arch", "enh_id", "core_remodeling"]].drop_duplicates()

    base_df["lin_len"] = base_df.end_lin - base_df.start_lin

    # apply core age
    base_df = pd.merge(base_df, core_age, how = "left", on = "enh_id")
    base_df.mrca =base_df.mrca.round(3)
    base_df = pd.merge(base_df, syn_gen_bkgd, how = "left", on = "mrca")

    # remove id where architecture wasnt assigned
    base_df = base_df.loc[base_df.code != ""]
    return base_df



def plot_hist(derived, core, simple):


    fig, ax = plt.subplots(figsize = (8, 8))

    sns.distplot(derived, hist_kws=dict(cumulative=True),
                     kde_kws=dict(cumulative=True), label = "derived", kde =False, norm_hist = True, color = "blue")
    sns.distplot(core, hist_kws=dict(cumulative=True),
                     kde_kws=dict(cumulative=True), label = "complex core", kde =False, norm_hist = True, color = "purple")
    sns.distplot(simple, hist_kws=dict(cumulative=True),
                     kde_kws=dict(cumulative=True), label = "simple", kde =False, norm_hist = True, color = "gold")

    k, kp = stats.kruskal(simple, derived, core)

    ax.set_title("FANTOM LINSIGHT scores\nCumulative Distribution")
    ax.set_xlabel("LINSIGHT score\nkruskal = %s, p = %s" % (k, kp))
    ax.set_ylabel("% of enhancer bases")
    ax.legend(bbox_to_anchor=(1.5, 1.0))

    plt.savefig("%sfantom_linsight_architecture.pdf" %(RE), bbox_inches = "tight")

    return k, kp

def get_expanded_arch_dfs(base_df):

    # separate the dataframes by architevture and expand to bp level
    simple_df = base_df.loc[base_df.code.str.contains("simple")]
    simple = np.repeat(simple_df.linsight_score, simple_df.lin_len) # expand linsight value for each simple basepair

    derived_df = base_df.loc[base_df.code.str.contains("derived")]
    derived = np.repeat(derived_df.linsight_score, derived_df.lin_len)

    core_df = base_df.loc[base_df.code.str.contains("core")]
    core = np.repeat(core_df.linsight_score, core_df.lin_len)

    complexenh = pandas.concat([derived, core])

    return simple, derived, core, complexenh


def get_stats(concat):

    medians = concat.groupby("code")["linsight_score"].median().reset_index()
    medians["measurement"] = "median_linsight"
    means = concat.groupby("code")["linsight_score"].mean().reset_index()
    means["measurement"] = "mean_linsight"

    measures = pd.concat([medians, means])
    measures.to_csv(f"{RE}fantom_linsight_score_summary.tsv")

    return measures


def get_counts(df, strat):

    if strat == 1:
        counts = df.groupby(["mrca_2", "core_remodeling"])["enh_id"].count().reset_index()

        # add empty dataframe for complex human enhancers (do not exist by our def)
        empty = pd.DataFrame({"mrca_2": [0.000],
        "core_remodeling": [1],
        "enh_id": [0]
        })

        counts = pd.concat([counts, empty]) # concat the dataframe

        counts =counts.sort_values(by = ["core_remodeling", 'mrca_2']).reset_index()

    else:
        counts = df.groupby("arch")["enh_id"].count().reset_index()
        counts =counts.sort_values(by = "arch", ascending = False).reset_index()

    counts["enh_id"] = counts["enh_id"].astype(int)
    counts = counts.drop(["index"], axis = 1)

    return counts


def plot_figure3(df, order, fig_id, re, trim_len, frac, dataset):

    xlab = ['Homo', 'Prim', 'Euar', 'Bore', 'Euth', 'Ther', 'Mam',
    'Amni', 'Tetr', 'Vert']

    if "98" in dataset:
        title = "98 fantom Datasets"
        # for plotting, get the median overlaps per arch, per dataset

        # for plotting age stratified, only plot 10% of the data
        smalldf = df.sample(frac = 0.05)
        sum_med = smalldf
    else:
        title = dataset
        sum_med = df
        smalldf = df

    # set up plot
    sns.set("poster")
    fig = plt.figure(figsize = (12, 8))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 3])
    ax = plt.subplot(gs[0])

    # plot the first panel
    x, y = "arch", "linsight_score"
    data = sum_med

    splot = sns.barplot( x = x, y = y, data = data,
                palette = palette, order = order,
                ax = ax)

    ax.yaxis.set_major_locator(MultipleLocator(4))

    # get counts for annotation
    STRAT = 0
    counts = get_counts(df, STRAT)
    for n, p in enumerate(splot.patches):
        value = counts.iloc[n]["enh_id"]
        splot.annotate(value,
                       (p.get_x() + p.get_width() / 2.,0.25),
                       ha = 'center', va = 'baseline',
                       size=15,
                       rotation = 90,
                       color = "white",
                       xytext = (0, 1),
                       textcoords = 'offset points'
                       )

    # plot the second panel
    ax2 = plt.subplot(gs[1])

    x, y = "mrca_2", "linsight_score"
    data = smalldf.sort_values(by = "mrca_2")
    hue = "arch"

    STRAT = 1
    agecounts = get_counts(df, STRAT)

    mplot = sns.barplot(x = x, y = y, data = data,
                hue = hue,
                palette = palette,
                ax = ax2)

    for n, p in enumerate(mplot.patches):

        value = agecounts.iloc[n]["enh_id"].astype(int)

        mplot.annotate(value,
                       (p.get_x() + p.get_width() / 2.,0.25),
                       ha = 'center', va = 'baseline',
                       size=15,
                       rotation = 90,
                       color = "white",
                       xytext = (0, 1),
                       textcoords = 'offset points'
                       )
    ax2.set_xticklabels(xlab, rotation = 90)
    ax2.set(xlabel = "", ylabel = "", title = title)

    ax2.legend().remove()

    ax2.yaxis.set_major_locator(MultipleLocator(4))

    ax2lim = ax2.get_ylim()

    ax.set(xlabel = "", ylabel = "linsight score", ylim = ax2lim)
    ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
    sns.set("poster")

    outf = f"{re}fig{fig_id}_pleiotropy-fantom_trim_{trim_len}_noexon_{dataset}_{frac}_mrcas.pdf"

    plt.savefig(outf, bbox_inches = "tight")


#%% expand linsight score estimates across bases

syn_gen_bkgd = load_syn_gen_bkgd(BUILD)

base_df = format_df(fantom_fs, syn_gen_bkgd)


base_df.head()
simple, derived, core, complexnh = get_expanded_arch_dfs(base_df)

k, kp = plot_hist(derived, core, simple)

#%%

simplef = simple.to_frame()
simplef["code"] = "simple"
simplef["arch"] = "simple"

derivedf = derived.to_frame()
derivedf["code"] = "derived"
derivedf["arch"] = "complex"

coref = core.to_frame()
coref["code"] = "complex_core"
coref["arch"] = "complex"

concat = pandas.concat([simplef, derivedf, coref])
concat = concat.sample(frac = 0.05)

# In[14]:

# In[16]:
measures = get_stats(concat)

order = ["simple", "complex_core", "derived"]
fig, ax = plt.subplots(figsize = (8, 8))
sns.barplot(x = "code", y = "linsight_score", data = concat,
            #showfliers = False,
            #notch = True,
            order = order,
           palette = arch_palette)
ax.set_xlabel("\nkruskal = %s, p = %s" % (k, kp))
#ax.set_title("LINSIGHT score by FANTOM architecture age")
ax.set_xlabel("")
ax.set_ylabel("LINSIGHT score")
sns.set("poster")
plt.savefig("%sfantom_linsight_architecture_boxplot.pdf" %(RE), bbox_inches = "tight")



order = ["simple", "complex_core", "derived"]
fig, ax = plt.subplots(figsize = (8, 8))
sns.barplot(x = "taxon2", y = "linsight_score", data = base_df.sort_values(by="mrca_2")
              , hue = "code", palette = arch_palette, hue_order = order)
ax.set_xticklabels(ax.get_xticklabels(), rotation = 90, horizontalalignment = "left")
ax.set_xlabel("Core age")
ax.set_title("LINSIGHT score by FANTOM architecture age")
ax.legend(bbox_to_anchor=(1.5, 1.0))
plt.savefig("%sfantom_linsight_architecture_mrca.pdf" %(RE), bbox_inches = "tight")





# In[25]:
ORDER = ["simple", "complex"]
FIG_ID = "3C"
DATASET = "all_fantom_enh"

plot_figure3(base_df, ORDER, FIG_ID, RE, TRIM_LEN, FRAC, DATASET)


#%%
from matplotlib import gridspec

order = ["simple", "complex"]
order1 = ["simple", "complex"]

fig = plt.figure(figsize = (12, 8))
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 3])
ax0 = plt.subplot(gs[0])

sns.barplot(x = "arch", y = "linsight_score", data = sample,
#            showfliers = False, notch = True,
            order = order1,
           palette = pal,
           ax = ax0)

ax0.set_xlabel("")
ax0.set_xticklabels("")
ax0.set_ylabel("LINSIGHT score")
ax0.set_ylim(0,0.67)
sns.set("poster")
sample2 = base_df.sample(frac =0.03)
ax2 = plt.subplot(gs[1])
sns.barplot(x = "taxon2", y = "linsight_score", data = sample2.sort_values(by="mrca_2"),
              hue = "arch",
              palette = pal,
              hue_order = order,
              #join = False,
              #dodge = 0.25,
             ax = ax2)

ax2.set_xticklabels(ax2.get_xticklabels(),rotation = 90, horizontalalignment = "left")
ax2.set_xlabel("architecture age")

ax2.legend().remove()
ax2.set_xlabel("")
ax2.set_xticklabels("")
#ax2.set_yticklabels("")
ax2.set_ylabel("")
ax2.set_ylim(0,0.67)
ax2.set_xlabel("")
sns.set("poster")
plt.savefig("%sFig3c-JOINT_fantom_linsight_architecture_mrca_enh.pdf" %(RE), bbox_inches = "tight")


# In[112]:


old.head()


# In[121]:


order = ["simple", "complex"]
fig, ax = plt.subplots(figsize= (8,8))
sns.boxplot( x = "arch", y = "linsight_score", data = sample2.loc[sample2.mrca_2> 0.175],
            notch = True, order = order, palette = pal)
old = sample2.loc[sample2.mrca_2> 0.175]
m, mp = stats.mannwhitneyu(old.loc[old.arch.str.contains("complex"), "linsight_score"],
                          old.loc[old.arch.str.contains("simple"), "linsight_score"])
print(m,mp)


# In[114]:


old.groupby("arch")["linsight_score"].mean()


# In[29]:


sample2.head()


# In[52]:


kw_list = []
for i in sample2.mrca_2.unique():
    mrca_list = sample2.loc[sample2.mrca_2 == i, "linsight_score"].to_list()
    kw_list.append(mrca_list)
from scipy.stats import mstats

args=[l for l in kw_list]
stats.mstats.kruskalwallis(*args)


# In[ ]:


core_age = df.groupby(["enh_id"])["mrca"].max().reset_index()
core_age = pandas.merge(core_age, syn_gen_bkgd) # add in taxon2
core_age = core_age[["enh_id", "mrca_2"]]
core_age.columns = ["enh_id", "mrca_2_core"]
core_age.sort_values(by="mrca_2_core").head()


# # Two-way ANOVA does not work because
# Independent observations - yes
# Residue distribution is normal - no (Jarque-bera)
# homogeneity of variance - no (Omnibus)

# In[55]:


sample2.head()


# In[108]:


import statsmodels.api as sm
from statsmodels.formula.api import ols
import statsmodels.stats.multicomp
sample3 = base_df.sample(frac =0.01)
# Fits the model with the interaction term
# This will also automatically include the main effects for each factor
model = ols('linsight_score ~ C(mrca_2)*C(arch)', sample3).fit()

# Seeing if the overall model is significant
print(f"Overall model F({model.df_model: .0f},{model.df_resid: .0f}) = {model.fvalue: .3f}, p = {model.f_pvalue: .4f}")

model.summary()


# In[110]:


subset = sample3[['mrca', 'arch',]].drop_duplicates()
tuples = [tuple(x) for x in subset.to_numpy()]
tuples

resid = model.resid
factor_groups = sample3.groupby(['mrca_2','arch'])

plt.figure(figsize=(6,6));
for values, group in factor_groups:
    i,j = values
    group_num = i*100  # for plotting purposes
    x = [group_num] * len(group)
    plt.scatter(x, resid[group.index],
            s=10)
plt.xlabel('Group');
plt.ylabel('Residuals');


# In[57]:


res = sm.stats.anova_lm(model, typ= 2)
res


# In[ ]:
