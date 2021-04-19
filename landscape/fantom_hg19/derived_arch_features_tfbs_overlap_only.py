import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy import stats
import seaborn as sns
import statsmodels
import statsmodels.api as sm
import subprocess


FANTOMPATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/all_fantom_enh/ages/"
FANTOMFILE = "syn_breaks_all_fantom_enh_ages.bed"
FANTOM = os.path.join(FANTOMPATH, FANTOMFILE)

FANTOM_TFBS_ONLY = f"{FANTOMPATH}enh_tfbs_only.txt"

SHUFPATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/shuffle/first_round_breaks"
SHUFFILE = "no-exon_syn_breaks_shuf-all_fantom_enh_ages.bed"
SHUFFILE = "noexon.bed"
SHUFF = os.path.join(SHUFPATH, SHUFFILE)


RE = "/dors/capra_lab/projects/enhancer_ages/landscape/results/fantom/"

#%%
colors = [ "amber", "faded green"]
palette = sns.xkcd_palette(colors)
sns.palplot(palette)

colors = [ "dusty blue", "greyish"]
es = sns.xkcd_palette(colors)
sns.palplot(es)

colors = [ "dusty purple", "grey"]
pur = sns.xkcd_palette(colors)
sns.palplot(pur)

colors = [ "amber", "greyish", "faded green", "grey"]
enhpal = sns.xkcd_palette(colors)
sns.palplot(enhpal)

colors = [ "amber", "greyish",  "dusty purple", "brown grey",  "windows blue", "bluey grey"]
archpal = sns.xkcd_palette(colors)
sns.palplot(archpal)


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


#%%


def add_arch_labels(df):

    df["arch"] = "complex_core"

    df.loc[df.core_remodeling ==0, "arch"] = "simple"
    df.loc[df.core ==0, "arch"] = "complex_derived"

    return df


def format_syndf(enh_age_file):

    syn_cols = ["chr_syn", "start_syn", "end_syn",
    "enh_id",
    "chr", "start", "end",
    "seg_index", "core_remodeling", "core",
    "mrca",]

    syn = pd.read_csv(enh_age_file, sep ='\t', header = None, names = syn_cols)

    syn["syn_id"] = syn.chr_syn + ":" + syn.start_syn.map(str) + "-" + syn.end_syn.map(str)

    syn["syn_len"] = syn.end_syn - syn.start_syn
    syn["enh_len"] = syn.end - syn.start


    # age and taxon file
    syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
    syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t') # read the file
    syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

    syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2"]] # whittle down the df
    syn["mrca"] = syn["mrca"].round(3) # round the ages

    syn = pd.merge(syn, syn_gen_bkgd, how = "left", on = "mrca")

    labeled_syn = add_arch_labels(syn) # add architecture labels

    return labeled_syn


def clean_shuffles(df):

    remove_list = []
    shuf_remove_ids = df.loc[(df["pct_enh"] ==1) & (df.core ==0)]

    df = df.loc[~df.enh_id.isin(shuf_remove_ids)]

    return df


#%%


enh = format_syndf(FANTOM)
enh["id"] = "FANTOM"

tfbs_overlap = pd.read_csv(FANTOM_TFBS_ONLY, sep = '\t')
tfbs_overlap.head()

enh.shape
enh = enh.loc[enh.enh_id.isin(tfbs_overlap.enh_id)]
enh.shape
shuf = format_syndf(SHUFF)

shuf["id"] = "SHUFFLE"


df = pd.concat([enh, shuf])
df.shape



#%% summarize architecture lengths per enhancer


sum_archlen = df.groupby(['id', 'enh_id', 'enh_len', 'core_remodeling', 'core'])["syn_len"].sum().reset_index().drop_duplicates()

sum_archlen = add_arch_labels(sum_archlen) # add architecture labels

 # calculate percent arch per enhancer
sum_archlen["pct_enh"] = sum_archlen.syn_len.divide(sum_archlen.enh_len)


# add oldest ages of enhancer information
enh_mrcas = df.groupby("enh_id")[["mrca_2", "taxon2"]].max().reset_index()
sum_archlen = pd.merge(sum_archlen, enh_mrcas, how = "left", on = "enh_id")

sum_archlen.head()
shuf_remove_ids = sum_archlen.loc[(sum_archlen["mrca_2"] == 0) & (sum_archlen.core ==0), "enh_id"]
sum_archlen = sum_archlen[~sum_archlen["enh_id"].isin(shuf_remove_ids)]
#%%

sum_archlen["arch_"] = "simple"
sum_archlen.loc[sum_archlen.core_remodeling ==1, "arch_"] = "complex"
sum_archlen["arch__"] = sum_archlen["arch_"] + "-" + sum_archlen['id']
sum_archlen["arch__"].unique()

#%% plot syn lengths per age, arch.


x = "mrca_2"
y = "syn_len"
data = sum_archlen.sample(frac = 0.25)
hue = "arch__"


xlabs = ["homo", "prim", "euar", "bore", "euth", "ther", "mam", "amni", "tetr", "vert"]
order = ["simple-FANTOM", 'simple-SHUFFLE', 'complex-FANTOM', 'complex-SHUFFLE']

fig, ax = plt.subplots(figsize = (12,6))
sns.barplot(x=x, y=y, data = data,
            hue = hue,
            palette = enhpal, hue_order = order )

ax.set(ylabel = "syntenic length", xlabel = "")
ax.set_xticklabels(xlabs, rotation = 90)
ax.legend(bbox_to_anchor = (1,1))
outf = "%sTFBS_ONLY_mrca_x_syn_lengths_arch.pdf" % RE

plt.savefig(outf, bbox_inches = 'tight')


#%% plot length of segments


ticklabs = ["simple", "complex\ncore", "complex\nderived"]

x = "arch"
y = "syn_len"
hue = "id"
data = sum_archlen
order = ["simple", "complex_core", "complex_derived"]

fig, ax = plt.subplots(figsize = (6,6))

sns.barplot(x=x, y=y, data = data, hue = hue, palette = es, order = order)
ax.set_xticklabels(ticklabs)
ax.legend(bbox_to_anchor = (1,1))
ax.set(ylabel = "sum length (bp)", xlabel = "")

plt.savefig("%s_TFBS_ONLY_sum_arch_length_all_fantom_shuffle.pdf" %RE, bbox_inches = "tight")


#%% mean lengths

comp = "fantom_tfbs_only_v_shuffle"

lens = sum_archlen.groupby(["id", "arch"])["syn_len"].mean().reset_index()
pcts = sum_archlen.groupby(["id", "arch"])["pct_enh"].mean().reset_index()

outf = f"{RE}{comp}_lens_metrics.tsv"
lens.to_csv(outf, sep = '\t')

outf = f"{RE}{comp}_pcts_metrics.tsv"
pcts.to_csv(outf, sep = '\t')

fantom_core_len = sum_archlen.loc[(sum_archlen.id == "YES_TFBS")& (sum_archlen["arch"] == "complex_core"), "syn_len"]

shuffle_core_len = sum_archlen.loc[(sum_archlen.id == "NO_TFBS")& (sum_archlen["arch"] == "complex_core"), "syn_len"]

fantom_der_len = sum_archlen.loc[(sum_archlen.id == "YES_TFBS")& (sum_archlen["arch"] == "complex_derived"), "syn_len"]

shuffle_der_len = sum_archlen.loc[(sum_archlen.id == "NO_TFBS")& (sum_archlen["arch"] == "complex_derived"), "syn_len"]


core_stat, core_p = stats.mannwhitneyu(fantom_core_len, shuffle_core_len)
# tfbs overlapping cores are longer than non-TFBS overlapping cores
# MannwhitneyuResult(statistic=6404937.5, pvalue=6.22796958050998e-72)

der_stat, der_p = stats.mannwhitneyu(fantom_der_len, shuffle_der_len)
# tfbs overlapping derived regions are longer than non-TFBS overlapping cores
# MannwhitneyuResult(statistic=7359068.0, pvalue=3.2945921502934066e-24)

mwu_df = pd.DataFrame({
"comparison": [comp, comp],
"mwu_analysis": ["core_w_TFBS_v_wo_TFBS", "der_w_TFBS_v_wo_TFBS"],
"stat": [core_stat, der_stat],
"P": [core_p, der_p]
})
outf = f"{RE}{comp}_MWU_core_der_lens.tsv"

mwu_df.to_csv(outf, sep = '\t')
#%%



x = "arch"
y = "pct_enh"
hue = "id"


sns.set("poster")

fig, ax = plt.subplots(figsize=(6,6))

# plot
sns.barplot(x = x, y = y, data = data, order = order,
hue = hue, palette =archpal, errwidth= 10)

ax.set(ylabel = "Percent of enhancer length", xlabel = "")
ax.set_xticklabels(ticklabs)
ax.legend(bbox_to_anchor = (1,1))
plt.savefig("%sTFBS_ONLY_fantom_percent_all.pdf" % RE, bbox_inches = "tight")




#%%
sum_archlen["id2"] = sum_archlen["arch"] + "-" + sum_archlen["id"]


x = "mrca_2"
y = "pct_enh"
hue = "id2"
order2 = ["simple-FANTOM", "simple-SHUFFLE",
"complex_core-FANTOM","complex_core-SHUFFLE",
 "complex_derived-FANTOM", "complex_derived-SHUFFLE"]

sns.set("poster")


fig, ax = plt.subplots(figsize=(12,6))
sns.barplot(x = x, y = y, data = data, hue = hue, hue_order = order2,
palette = archpal, errwidth= 5)

ax.set(ylabel = "Percent of enhancer length", xlabel = "")

ax.set_xticklabels(xlabs, rotation = 90)
ax.legend(bbox_to_anchor = (1,1))
plt.savefig("%sTFBS_ONLY_fantom_percent_arch_mrca_2.pdf" % RE, bbox_inches = "tight")

#%%

x = "mrca_2"
y = "syn_len"
hue = "id2"


sns.set("poster")


fig, ax = plt.subplots(figsize=(12,6))
sns.barplot(x = x, y = y, data = data, hue = hue,
hue_order = order2,  palette = archpal, errwidth= 5)

ax.set(ylabel = "syntenic length", xlabel = "")
ax.set_xticklabels(xlabs, rotation = 90)
ax.legend(bbox_to_anchor = (1,1))
plt.savefig("%sTFBS_ONLY_fantom_syn_len_mrca_2.pdf" % RE, bbox_inches = "tight")
#%%

shuf_remove_ids = sum_archlen.loc[(sum_archlen["mrca_2"] == 0) & (sum_archlen.core ==0), "enh_id"]

test = "chr14:89971043-89971729"
sum_archlen.loc[sum_archlen.enh_id == test]
shuf.loc[shuf.enh_id == test]
testsyn = "chr14:89971596-89971729"
shuf.loc[shuf.syn_id == testsyn]
#%%

enh[["core_remodeling", "enh_id"]].drop_duplicates().groupby("core_remodeling")["enh_id"].count()
"""
core_remodeling
0    13684
1     8683

Total = 22367
% simple = 61.2%
Name: enh_id, dtype: int64
"""

#%%
def MRCA_frequency(catdf, cols, var):

    age_dict = {} # collect age frequency results per dataset
    summary_age_dict = {} # collect summarized age frequencies per dataset

    for n, dataset in enumerate(catdf["id"].unique()):
        # count n enhancers in architecture per age
        test = catdf.loc[catdf["id"] == dataset]

        age = test.groupby(cols)["enh_id"].count().reset_index()

        # rename columns
        age.columns = cols + ["counts"]

        # sum total n enhancers in architecture
        cols_no_var = list(set(cols) - set([var]))
        totals = age.groupby(cols_no_var)["counts"].sum().reset_index()
        # rename columns
        totals.columns = cols_no_var + ["total_id"]

        # merge dataframes
        age = pd.merge(age, totals, how = "left")

        # calculate the % of architecture in each age
        age["freq"] = age.counts.divide(age.total_id)

        age["dataset_name"] = dataset
        age_dict[n] = age

        # summarize frequencies across architectures, before/after eutherian.
        eutherian = age.loc[age["mrca_2"] == 0.19][[ "id", "freq"]]
        eutherian["category"] = "eutherian"

        younger_thaneuth = age.loc[age["mrca_2"] <0.19].groupby(["id"])["freq"].sum().reset_index()
        younger_thaneuth["category"] = "younger than eutherian"

        older_thaneuth = age.loc[age["mrca_2"] >0.19].groupby(["id"])["freq"].sum().reset_index()
        older_thaneuth["category"] = "older than eutherian"

        summarized_freq = pd.concat([eutherian, younger_thaneuth, older_thaneuth])
        summarized_freq["dataset_name"] = dataset

        summary_age_dict[n] = summarized_freq

    # concat age and summarized frequency dataframes
    ages = pd.concat(age_dict.values())
    summarized_freq = pd.concat(summary_age_dict.values())

    # calculate fold-change of enh v. shuf expectation per shuffle


    # select only the enhancer and specific shuffle instance
    enhdf = ages.loc[ages["id"] == "FANTOM"]

    shuf_ = ages.loc[ages["id"] != "FANTOM"]

    merge_cols = list(set(cols) - set(["id"]))

    fc = pd.merge(shuf_, enhdf, how = "left", on =merge_cols)

    # calculate fold changes
    fc["fold_change"] = fc["freq_y"].divide(fc["freq_x"])


    col_id = "_".join(cols)
    outf = f'{RE}{col_id}_freq.txt'
    ages.to_csv(outf, sep = '\t', index = False)

    outf = f'{RE}{col_id}_fold_change.txt'
    fc.to_csv(outf, sep = '\t', index = False)

    outf = f'{RE}summary_{col_id}_freq.txt'
    summarized_freq.to_csv(outf, sep = '\t', index = False)



    return ages, fc

def plot_arch_freq(age_arch_freq, age_freq):
    plots = {"age_arch_tfbs_only" : age_arch_freq, "age_tfbs_only": age_freq}

    for name, frame in plots.items():

        if name == "age_arch_tfbs_only": # arrange order and colors of plot.
            frame["plot_hue"] = frame["arch"].astype(str) + "-" + frame["id"].astype(str)
            order = ["simple-FANTOM", "simple-SHUFFLE",
            "complex_core-FANTOM", "complex_core-SHUFFLE",
            "complex_derived-FANTOM", "complex_derived-SHUFFLE"]
            hue = "plot_hue"

        else:
            order = ["FANTOM", "SHUFFLE"]
            hue = "id"

        if GENOME_BUILD == "hg38":
            xlabs = ["Prim", "Euar", "Bore", "Euth", "Ther", "Mam", "Amni", "Tetr", "Sarg", "Vert"] # set xlabels
        else:
            xlabs = ["Homo", "Prim", "Euar", "Bore", "Euth", "Ther", "Mam", "Amni", "Tetr", "Vert"]

        sns.set("talk")
        fig, ax = plt.subplots(figsize = (6,6))
        x, y = "mrca_2", "freq"
        data = frame

        sns.barplot(x = x, y=y,
        data = data,
        hue = hue,
        hue_order = order,
        palette = archpal)

        ax.set_xticklabels(xlabs, rotation = 90)
        ax.legend(bbox_to_anchor = (1,1))

        outf = f"{RE}{name}_freq_per_age.pdf"

        plt.savefig(outf, bbox_inches= "tight")

def plot_arch_fc(age_arch_fc, age_fc, arch):

    plots = {"age_arch_tfbs_only":age_arch_fc, "age_tfbs_only": age_fc}

    for name, fc in plots.items():

        fc['log2'] = np.log2(fc["fold_change"])

        archs = ["simple", "complex_core", "complex_derived"]

        if name == "age_arch_tfbs_only" and arch = "all":
            order = ["simple", "complex_core", "complex_derived"]
            hue = "arch"
            data = fc

        if name == "age_arch_tfbs_only" and arch in archs:
            order = [arch]
            hue = "arch"
            data = fc.loc[fc.arch == arch]


        else:
            arch = None
            order = ["FANTOM"]
            hue = "id_y"


        if GENOME_BUILD == "hg38":
            xlabs = ["Prim", "Euar", "Bore", "Euth", "Ther", "Mam", "Amni", "Tetr", "Sarg", "Vert"]
        else:
            xlabs = ["Homo", "Prim", "Euar", "Bore", "Euth", "Ther", "Mam", "Amni", "Tetr", "Vert"]

        sns.set("talk")
        fig, ax = plt.subplots(figsize = (6,6))
        x, y = "mrca_2", "log2"


        sns.barplot(x = x, y=y,
        data = data,
        hue = hue,
        hue_order = order,
        palette = PAL)
        ax.set(ylabel = "Fold-Change v. Bkgd\n(log2-scaled)")
        ax.set_xticklabels(xlabs, rotation = 90)
        ax.legend(bbox_to_anchor = (1,1))
        outf = f"{RE}{name}_{arch}_fold_change_per_age.pdf"


    plt.savefig(outf, bbox_inches= "tight")


#%%
cols = ["id", "arch", "mrca_2"]
var = "mrca_2"
age_arch_freq, age_arch_fc = MRCA_frequency(df, cols, var)


cols = ["id", "mrca_2"]
var = "mrca_2"
age_freq, age_fc = MRCA_frequency(df, cols, var)

#%%
GENOME_BUILD = "hg19"
plot_arch_freq(age_arch_freq, age_freq)
plot_arch_fc(age_arch_fc, age_fc, "all")
