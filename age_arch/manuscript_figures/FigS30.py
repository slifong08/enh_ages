import glob
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.ticker as ticker
import numpy as np
import os, sys
import pandas as pd
import seaborn as sns
from scipy import stats


PATH = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/results/stats/"
RE = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/results/"

#%% Files
ARCH_STATS = glob.glob(f"{PATH}trim_310_E*_ARCH_metrics.txt")
FREQS = glob.glob(f"{PATH}trim_310_E*_id_mrca_2_freq.txt")
FCS = glob.glob(f"{PATH}trim_310_E*_id_mrca_2_freq.txt")

ARCH_STATS = glob.glob(f"{PATH}trim_310_E*_AGE_ARCH_metrics.txt")
ARCH_FREQS = glob.glob(f"{PATH}trim_310_E*_id_arch_mrca_2_freq.txt")
ARCH_FCS = glob.glob(f"{PATH}trim_310_E*_id_arch_mrca_2_fold_change.txt")

MRCA_OR = glob.glob(f"{PATH}trim_310_E*_summary_age_arch_OR.txt")
SEG_OR = glob.glob(f"{PATH}trim_310_E*_summary_seg_OR.txt")

data = [ARCH_STATS, FREQS, FCS, ARCH_STATS, ARCH_FREQS, ARCH_FCS, MRCA_OR, SEG_OR]
for i in data:

    print(len(i))

#%%
ROADMAP_ID = "trimmed_310"
FIG_ID = "S30"
BUILD = "hg19"

colors = [ "amber", "dusty purple", "windows blue"]
PAL = sns.xkcd_palette(colors)
sns.palplot(PAL)

colors = [ "windows blue"]
DERPAL = sns.xkcd_palette(colors)
sns.palplot(DERPAL)

colors = ["amber", "greyish", "faded green", "slate grey"]
ESPAL = sns.xkcd_palette(colors)
sns.palplot(ESPAL)

colors = ["amber", "faded green"]
EPAL = sns.xkcd_palette(colors)
sns.palplot(EPAL)

colors = ["slate grey", "greyish"]
Eshuf = sns.xkcd_palette(colors)
sns.palplot(Eshuf)

#%%


def make_df(f_list):

    f_dict = {} # collect each file's dataframe

    for f in f_list:
        if "summary_id" not in f:

            sid = "E" + (f.split("E")[1]).split("_")[0] # get the sample id

            df = pd.read_csv(f, sep = '\t')
            df["sid"] = sid
            f_dict[sid] = df

    outdf = pd.concat(f_dict.values())

    return outdf


def load_syn_gen_bkgd(build):

    F = f"/dors/capra_lab/projects/enhancer_ages/{build}_syn_taxon.bed"
    syngenbkgd = pd.read_csv(F, sep='\t')
    syngenbkgd[["mrca", "mrca_2"]] = syngenbkgd[["mrca", "mrca_2"]].round(3)

    return syngenbkgd


def format_df(df):
    #syngenbkgd = load_syn_gen_bkgd(BUILD)
    if "mrca_2" in list(df):
        df.mrca_2 = df.mrca_2.round(3)
    if "arch" in list(df):
        df["core_remodeling"] = 0
        df.loc[df.arch == "complexenh", "core_remodeling"] = 1

    if "OR" in list(df):
        df["yerr"] = df["ci_upper"] - df["ci_lower"]
        df["log2"] = np.log2(df["OR"])
    return df


def get_counts(df, groupby_cols, groupby_val):

    counts = df.groupby(groupby_cols)[groupby_val].sum().reset_index()

    if "mrca_2" in groupby_cols and "core_remodeling" in groupby_cols:

        # add empty dataframe for complex human enhancers (do not exist by our def)
        empty = pd.DataFrame({"mrca_2": [0.000],
        "core_remodeling": [1],
        groupby_val: [0]
        })

        counts = pd.concat([counts, empty]) # concat the dataframe

    # change count data to int
    counts[groupby_val] = counts[groupby_val].astype(int)

    # sort and reset the index. Seaborn plots by index value.
    counts = counts.sort_values(by = groupby_cols).reset_index()

    # drop the index column.
    counts = counts.drop(["index"], axis = 1)

    print(counts.head())
    return counts


def plot_annotate_counts(splot, counts_df, groupby_val, height_adjust):

    # annotate plot with counts
    for n, p in enumerate(splot.patches):

        value = counts_df.iloc[n][groupby_val].astype(int)
        if height_adjust == 0 and p.get_height() > 0:
            height_adjust = (p.get_height()+0.03)


        splot.annotate(value,
                       (p.get_x() + p.get_width() / 2.,height_adjust),
                       ha = 'center', va = 'baseline',
                       size=15,
                       rotation = 90,
                       color = "k",
                       xytext = (0, 1),
                       textcoords = 'offset points'
                       )


def plot_figS2xA_OR_seg(df, fig_id, roadmap_id, countdf, groupby_val):

    x, y = "seg_index", "log2", # plot the median length per dataset
    data = df.sort_values(by = "seg_index")


    fig, ax = plt.subplots(figsize =(6,6))
    splot = sns.barplot(x=x, y=y, data = data,
    linewidth=2.5, facecolor=(1, 1, 1, 0), edgecolor=".2",)
    #yerr =data["yerr"])

    ax.set(ylabel = "Odds Ratio", xlabel = "Number of Age Segments", title = roadmap_id)
    #ax.set_xticklabels(xlabs, rotation = 90)
    ax.legend(bbox_to_anchor = (1,1))
    plt.axhline(0, color = "grey", linewidth = 2.5)
    # annotate w/ counts
    height_adjust = 0.05
    plot_annotate_counts(splot, countdf, groupby_val, height_adjust)

    ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(round(2**x, 2)))

    ax.yaxis.set_major_formatter(ticks)
    ax.yaxis.set_major_locator(MultipleLocator(0.1))


    outf = f"{RE}{fig_id}A_{ROADMAP_ID}_seg_OR.pdf"
    plt.savefig(outf, bbox_inches = "tight")


def plot_figS2xC_len_mrca(arch, fig_id, roadmap_id, countdf, groupby_val, order, hue, palette):

    xlabs = ["Homo", "Prim", "Euar", "Bore", "Euth",
    'Ther', "Mam", "Amni", "Tetr", "Vert"]

    x, y = "mrca_2", "50%", # plot the median length per dataset
    data = arch.loc[arch.val == "enh_len"].sort_values(by = "mrca_2")
    hue = hue

    fig, ax = plt.subplots(figsize =(6,6))
    splot = sns.barplot(x=x, y=y, data = data, hue = hue, palette = palette, hue_order = order)

    ax.set(ylabel = "length (bp)", xlabel = "Sequence Age", title = roadmap_id)
    ax.set_xticklabels(xlabs, rotation = 90)
    ax.legend(bbox_to_anchor = (1,1))

    # annotate w/ counts
    height_adjust = 1
    plot_annotate_counts(splot, countdf, groupby_val, height_adjust = 1)

    outf = f"{RE}{fig_id}C_{ROADMAP_ID}_seg_OR.pdf"
    plt.savefig(outf, bbox_inches = "tight")


def plot_figS2xD_freq_mrca(freq, fig_id, roadmap_id, countdf, groupby_val, order, hue, palette):

    xlabs = ["Homo", "Prim", "Euar", "Bore", "Euth",
    'Ther', "Mam", "Amni", "Tetr", "Vert"]

    x, y = "mrca_2", "freq", # plot the median length per dataset
    data = freq.sort_values(by = "mrca_2")
    hue = hue

    fig, ax = plt.subplots(figsize =(6,6))
    splot = sns.barplot(x=x, y=y, data = data, hue = hue, palette = palette, hue_order = order)

    ax.set(ylabel = "% of architecture", xlabel = "Sequence Age", title = roadmap_id)
    ax.set_xticklabels(xlabs, rotation = 90)
    ax.legend(bbox_to_anchor = (1,1))

    # annotate w/ counts
    height_adjust = 0
    plot_annotate_counts(splot, countdf, groupby_val, height_adjust)
    fig.set_size_inches(10, 10)
    plt.tight_layout()
    outf = f"{RE}{fig_id}D_{ROADMAP_ID}_freq_arch_mrca2.pdf"
    plt.savefig(outf, bbox_inches = "tight")

#%%

ARCH = make_df(ARCH_STATS)
ARCH = format_df(ARCH)
df = ARCH
df.head()

#%%

# get counts for plot
groupby_cols = ["id", "mrca_2"]
groupby_val = "count"
countdf = get_counts(df, groupby_cols, groupby_val)
FIG_ID = "S23"
order = ["enh", "shuf"]
hue = "id"
palette = Eshuf
# plot
plot_figS2xC_len_mrca(df, FIG_ID, ROADMAP_ID, countdf, groupby_val, order, hue, palette)
#%%
ARCH_ENH = ARCH.loc[(ARCH.id == "enh") & (ARCH.val == "enh_len")]

df = ARCH_ENH

# get counts for plot
groupby_cols = ["core_remodeling", "mrca_2"]
groupby_val = "count"
countdf = get_counts(df, groupby_cols, groupby_val)

hue = "core_remodeling"
FIG_ID = "S30"
order = [0,1]
pal = EPAL

# plot
plot_figS2xC_len_mrca(df, FIG_ID, ROADMAP_ID, countdf, groupby_val, order, hue, pal)
#%%
len(ARCH_FREQS)

ARCHFQ = make_df(ARCH_FREQS)
ARCHFQ = format_df(ARCHFQ)
df = ARCHFQ.loc[ARCHFQ['id']== "enh"]
df.head()

#%%

groupby_cols = ["core_remodeling", "mrca_2"]
groupby_val = "counts"
countdf = get_counts(df, groupby_cols, groupby_val)

hue = "core_remodeling"
FIG_ID = "S30"
order = [0,1]
pal = EPAL

plot_figS2xD_freq_mrca(df, FIG_ID, ROADMAP_ID, countdf, groupby_val, order, hue, pal)
#%%
len(SEG_OR)

SEGORS = make_df(SEG_OR)
SEGORS = format_df(SEGORS)
df = SEGORS
df
#%%

groupby_cols = ["seg_index"]
groupby_val = "a"
countdf = get_counts(df, groupby_cols, groupby_val)


#%%
first_five = df.loc[df.seg_index<7]
first_five.loc[first_five.seg_index ==1]["OR"].mean()
plot_figS2xA_OR_seg(first_five, FIG_ID, ROADMAP_ID, countdf, groupby_val,)
#%%
