#!/usr/bin/env python
# coding: utf-8

# In[1]:


import glob
import pandas

from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator
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

pal = [ "amber", "faded green"]
palette = sns.xkcd_palette(pal)
sns.palplot(palette)

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


linsight_path = "/dors/capra_lab/projects/enhancer_ages/encode/data/ELS_combined_HepG2/ages/linsight/"

# enh only
enh_fs = f"{linsight_path}syn_breaks_ELS_combined_HepG2_ages_x_linsight.bed"



# In[3]:

BUILD = "hg19"
TRIM_LEN = "310"
FRAC = 0.5

def load_syn_gen_bkgd(build):

    F = f"/dors/capra_lab/projects/enhancer_ages/{build}_syn_taxon.bed"
    syngenbkgd = pd.read_csv(F, sep='\t')
    syngenbkgd[["mrca", "mrca_2"]] = syngenbkgd[["mrca", "mrca_2"]].round(3)

    return syngenbkgd


def assign_architecture(df):

    df["code"] = ""
    df.loc[(df.core_remodeling == 0)& (df.core == 1), "code"] = "simple"
    df.loc[(df.core_remodeling == 1)& (df.core == 1), "code"] = "complex_core"
    df.loc[(df.core_remodeling == 1)& (df.core == 0), "code"] = "derived"

    df["arch"] = ""
    df.loc[(df.core_remodeling == 0), "arch"] = "simple"
    df.loc[(df.core_remodeling == 1), "arch"] = "complex"

    # create rank system for simple, core and derived segments
    df["core_remodeling_2"] = df.core_remodeling
    df.loc[df["core"] == 0, "core_remodeling_2"] = 2
    return df


def load_fantom_encode_tfbs_overlap():
    SYN_OVERLAP = '/dors/capra_lab/projects/enhancer_ages/landscape/results/fantom/encode3/ALL_FANTOM_SYN_TF_OVERLAP.tsv'
    syn_tfbs_overlap = pd.read_csv(SYN_OVERLAP, sep = '\t')

    return syn_tfbs_overlap

def load_untrimmed_fantom_ids():
    map_trim = "/dors/capra_lab/projects/enhancer_ages/fantom/data/trimmed-310_all_unique_fantom_erna_112_tissue.bed"
    trim_map = pd.read_csv(map_trim, sep = '\t', header =None)
    cols = ["chr_en"]


def format_df(enh_fs, syn_gen_bkgd):

    arch_id = (enh_fs.split("/")[-1]).split("_")[1]
    print(arch_id)
    df = pandas.read_csv(enh_fs, sep = '\t', header = None, low_memory=False)

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

    base_df = df[["chr_lin", "start_lin", "end_lin", "linsight_score", "code", "arch", "enh_id", "core_remodeling", "core_remodeling_2"]].drop_duplicates()

    base_df["lin_len"] = base_df.end_lin - base_df.start_lin

    # apply core age
    base_df = pd.merge(base_df, core_age, how = "left", on = "enh_id")
    base_df.mrca =base_df.mrca.round(3)
    base_df = pd.merge(base_df, syn_gen_bkgd, how = "left", on = "mrca")

    # remove id where architecture wasnt assigned
    base_df = base_df.loc[base_df.code != ""]

    # column for counting things.
    base_df["counts"] = 1
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


def make_empty_dict(df, groupby_cols):
    # for when i need an empty dictionary with a complete set of architectures and values
    emptydf_dict = {}
    val = 0 # unique identifier
    if len(groupby_cols) == 2:
        for mrca_2 in df[groupby_cols[0]].unique():
            for arch in df[groupby_cols[1]].unique():

                emptydf = pd.DataFrame({ groupby_cols[0]:[mrca_2], groupby_cols[1]:[arch],})
                emptydf_dict[val] = emptydf
                val+=1
    elif len(groupby_cols) == 1:
        for mrca_2 in df[groupby_cols[0]].unique():
            emptydf = pd.DataFrame({ groupby_cols[0]:[mrca_2]})
            emptydf_dict[val] = emptydf
            val+=1
    empty = pd.concat(emptydf_dict.values())
    return empty


def get_counts(df, groupby_cols, groupby_val):

    counts = df.groupby(groupby_cols)[groupby_val].sum().reset_index()

    if "mrca_2" in groupby_cols and "core_remodeling" in groupby_cols:

        empty = make_empty_dict(df, groupby_cols) # make an empty df to fill architectures w/ no counts
        counts = pd.merge(empty, counts, how = "left", on = groupby_cols).fillna(0)

    elif "mrca_2" in groupby_cols and "core_remodeling_2" in groupby_cols:

        empty = make_empty_dict(df, groupby_cols) # make an empty df to fill architectures w/ no counts
        counts = pd.merge(empty, counts, how = "left").fillna(0)

    # change count data to int
    counts[groupby_val] = counts[groupby_val].astype(int)

    # sort and reset the index. Seaborn plots by index value.
    counts = counts.sort_values(by = groupby_cols).reset_index()

    # drop the index column.
    counts = counts.drop(["index"], axis = 1)

    return counts


def plot_annotate_counts(splot, counts_df, groupby_val, height_adjust):

    # annotate plot with counts
    for n, p in enumerate(splot.patches):

        value = counts_df.iloc[n][groupby_val].astype(int)
        #print(n, p, value)
        if height_adjust == 0:
            height_adjust = (p.get_height()-0.03)


        splot.annotate(value,
                       (p.get_x() + p.get_width() / 2.,height_adjust),
                       ha = 'center', va = 'baseline',
                       size=15,
                       rotation = 90,
                       color = "k",
                       xytext = (0, 1),
                       textcoords = 'offset points'
                       )


def plot_figure3(df, fig_id, re, trim_len, frac, dataset, x):

    # for ranking and stratifying by age and architecture, set x1 for counts
    if x == 'arch':
        x1 = "core_remodeling"
        pal = palette
        order_labs = ["simple", "complex"]
        order = [0,1]
    elif x == "code":
        x1 = "core_remodeling_2"
        pal = arch_palette
        order_labs = ["simple", "core", "derived"]
        order = [0,1,2]

    xlab = ['Homo', 'Prim', 'Euar', 'Bore', 'Euth', 'Ther', 'Mam',
    'Amni', 'Tetr', 'Vert']

    title = dataset

    # set up plot
    sns.set("poster")
    fig = plt.figure(figsize = (12, 8))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 3])
    ax = plt.subplot(gs[0])

    # plot the first panel
    y = "linsight_score"
    data = df

    splot = sns.barplot(x = x1, y = y, data = data,
                palette = pal, order = order,
                ax = ax)

    ax.yaxis.set_major_locator(MultipleLocator(0.1))
    ax.set_xticklabels(order_labs, rotation = 90)

    # get counts for annotation
    groupby_cols, groupby_val = [x1], "lin_len"
    countdf = get_counts(df, groupby_cols, groupby_val)
    height_adjust = 0.01
    plot_annotate_counts(splot, countdf, groupby_val, height_adjust)

    # plot the second panel
    ax2 = plt.subplot(gs[1])

    x, y = "mrca_2", "linsight_score"
    data = df.sort_values(by = "mrca_2")
    hue = x1

    # plot mrca-stratified barplot

    mplot = sns.barplot(x = x, y = y, data = data,
                hue = hue,
                palette = pal,
                ax = ax2)

    # add counts and annotate barplot
    groupby_cols, groupby_val = [x1, "mrca_2"], "lin_len"
    countdf = get_counts(df, groupby_cols, groupby_val)
    plot_annotate_counts(mplot, countdf, groupby_val, height_adjust)

    # plot x labels
    ax2.set_xticklabels(xlab, rotation = 90)
    ax2.set(xlabel = "", ylabel = "", title = title)

    ax2.legend().remove()

    ax2.yaxis.set_major_locator(MultipleLocator(0.1))

    ax2lim = ax2.get_ylim()
    # set first plot lim based on second plot lim
    ax.set(xlabel = "", ylabel = "linsight score", ylim = ax2lim)

    outf = f"{re}fig{fig_id}_LINSIGHT_{trim_len}_noexon_{dataset}_{frac}_mrcas.pdf"

    plt.savefig(outf, bbox_inches = "tight")


#%% expand linsight score estimates across bases

syn_gen_bkgd = load_syn_gen_bkgd(BUILD)

syn_tfbs_overlap = load_fantom_encode_tfbs_overlap() # merge infor about ENCODE3 TFBS overlap

base_df = format_df(enh_fs, syn_gen_bkgd)

base_df.head()

base_df = pd.merge(base_df, syn_tfbs_overlap, how = "left")
#%%

# expand linsight counts to per basepair level
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


#%%
measures = get_stats(concat)


#%%
FIG_ID = "3C"
DATASET = "all_fantom_enh"
x = "arch"

plot_figure3(base_df, FIG_ID, RE, TRIM_LEN, FRAC, DATASET, x)

#%% plot only enhancers that overlap TFBS in ENCODE

FIG_ID = "3C"
DATASET = "TFBS_overlap_only_all_fantom_enh"
x = "arch"
data = base_df.loc[base_df.tfoverlap_bin ==1]
data.shape
plot_figure3(data, FIG_ID, RE, TRIM_LEN, FRAC, DATASET, x)


#%%


FIG_ID = "3C"
DATASET = "NO_TFBS_all_fantom_enh"
x = "arch"
no_data = base_df.loc[base_df.tfoverlap_bin ==0]

plot_figure3(no_data, FIG_ID, RE, TRIM_LEN, FRAC, DATASET, x)


#%%

FIG_ID = "3A"
DATASET = "all_fantom_enh"
x = "code"
plot_figure3(base_df, FIG_ID, RE, TRIM_LEN, FRAC, DATASET, x)

#%% evaluate linsight on syntenic architecture w/ evidence for TFBS binding in encode


FIG_ID = "3A"
DATASET = "TFBS_overlap_only_all_fantom_enh"
x = "code"
data = base_df.loc[base_df.tfoverlap_bin ==1]
plot_figure3(data, FIG_ID, RE, TRIM_LEN, FRAC, DATASET, x)


#%% evaluate linsight on syntenic architecture w/ no evidence for TFBS binding


FIG_ID = "3A"
DATASET = "NO_TFBS_all_fantom_enh"
x = "code"
nodata = base_df.loc[base_df.tfoverlap_bin ==0]
plot_figure3(nodata, FIG_ID, RE, TRIM_LEN, FRAC, DATASET, x)


#%%

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
