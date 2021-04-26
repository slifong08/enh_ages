#!/usr/bin/env python
# coding: utf-8

# In[1]:


import glob
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator
import numpy as np
import os, sys
from scipy import stats
import seaborn as sns
order = ["simple", "complex"]

RE = "/dors/capra_lab/projects/enhancer_ages/fantom/results/for_publication/pleiotropy/"

# sns colors
arch_colors = ["amber", "dusty purple", "windows blue","greyish"]
arch_palette = sns.xkcd_palette(arch_colors)


colors = [ "amber", "faded green", "dusty purple", "windows blue","greyish"]
palette = sns.xkcd_palette(colors)
sns.palplot(arch_palette)


FRAC = 0.5
TRIM_LEN = "310"
BUILD = "hg19"
FIG_ID = "3"

# path to files.

PATH ="/dors/capra_lab/projects/enhancer_ages/fantom/data/download/tfbs/"
#age_fs = glob.glob("%s*_enh_tfbs_density.bed" % path)

PATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/non-genic/no-exon_all_fantom_enh/breaks/"
age_fs = glob.glob(f"{PATH}no-exon_all_fantom_enh_ages_enh_age_arch_summary_matrix.bed")
# get multiintersected files - each has been trimmed to dataset mean,
#exons excluded, multiintersected with 97 other ROADMAP enhancers.

MULTIPATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/non-genic/trimmed/multiintersect/"


multi_fs = glob.glob(f"{MULTIPATH}trim_310_all_fantom_enh_multiintersect_0.5_count.bed")
# raw enhancers have been aged, and number of age segments counted.

print(age_fs, multi_fs)


# other annotations about datasets descriptions

#%% functions

def load_syn_gen_bkgd(build):

    F = f"/dors/capra_lab/projects/enhancer_ages/{build}_syn_taxon.bed"
    syngenbkgd = pd.read_csv(F, sep='\t')
    syngenbkgd[["mrca", "mrca_2"]] = syngenbkgd[["mrca", "mrca_2"]].round(3)

    return syngenbkgd


def load_desc_df():
    desc_file = "/dors/capra_lab/data/fantom/fantom5/facet_expressed_enhancers/sample_id_descriptions.txt"
    desc_df = pd.read_csv(desc_file, sep = '\t', header = None)

    desc_df.columns = ["sid", "desc"]

    return desc_df


def get_sample_dicts(multi_fs, age_fs):

    multi = {}
    age = {}


    for f in multi_fs:

        sid = (f.split("/")[-1]).split(".")[0]
        sid = "all_fantom_enh"
        multi[sid] = f


    for f in age_fs:

        sid = (f.split("/")[-1]).split(".")[0]
        sid = "all_fantom_enh"
        age[sid] = f

    return multi, age


def format_df_get_stats(multi, age, RE, trim_len):

    results_data = {}
    results_stats ={}

    for key, multi_f in multi.items():

        # open the multiintersect df
        multi = pd.read_csv(multi_f, sep = '\t', header = None, usecols = [3, 4, 5, 6, 7])
        m_cols = ["enh_id", "old_len", "trimmed_enh_id",  "mean_len", "count_overlap"]
        multi.columns = m_cols


        # open age file
        age_f = age[key]
        agedf = pd.read_csv(age_f, sep = '\t', usecols = [3, 5, 6, 7, 8, 10, 11]) # open the count_overlap df
        a_cols  = ["enh_id",  "seg_index", "core_remodeling", "arch", "mrca",  "mrca_2", "taxon2"]
        #agedf.columns = a_cols

        # merge age and multi on enh_id (old enhancer id, not trimmed id)
        # old id to preserve the architecture label of the original enhancer
        # while evaluating its length-controlled overlap among other tissues.
        #print(multi.head(), agedf.head())
        enh = pd.merge(multi, agedf, how = "left", on = "enh_id")

        # clean up long enhancers (>10kb), and assign relative simple

        enh = enh.loc[enh.old_len<10000]
        enh["sid"] = key
        print(enh.head())
        # assign dataset id

        results_data["all_fantom_enh"] = enh

        # MWU simple v. complex overlap
        enh.loc[enh["arch"] != "simple", "arch"] = "complex"
        complex_tissue_overlap = enh.loc[enh["arch"] != "simple", "count_overlap"]
        simple_tissue_overlap =enh.loc[enh["arch"] == "simple", "count_overlap"]
        t, pval = stats.mannwhitneyu(complex_tissue_overlap, simple_tissue_overlap)
        #print(DESC, t, pval)

        # get medians
        medians = enh.groupby(['arch'])["count_overlap"].median().reset_index()
        medians["measure"] = "median"

        # get means
        means = enh.groupby(['arch'])["count_overlap"].mean().reset_index()
        means["measure"] = "mean"

        # create a dataframe of stats
        measurements = pd.concat([medians, means])
        measurements["t-stat"],  measurements["p"],  measurements["sid"]= [t, pval, key]

        results_stats[key] = measurements

    # concat dataframes
    res_data = pd.concat(results_data.values())
    res_stats = pd.concat(results_stats.values())

    # write the metrics out
    outf = f"{RE}tissue_overlap_measurements_{trim_len}.tsv"
    res_stats.to_csv(outf, sep = '\t', index =False)

    return res_data, res_stats


def plot_per_dataset_tissue_overlap(res, order, re, trim_len, frac):
    sum_med = res.sample(frac = 0.1)

    sns.set("poster")
    fig, ax = plt.subplots(figsize = (8,40))

    x, y= "count_overlap","arch"
    data = sum_med.sort_values(by = "count_overlap")


    sns.barplot(x = x, y =y, data = data,
           palette = palette, order = order)

    ax.legend().remove()
    ax.set_xlabel("tissue overlap count")
    outf = f"{re}pleiotropy_per_fantom_trimmed_{trim_len}_{frac}.pdf"
    plt.savefig(outf, bbox_inches = "tight")


def plot_per_dataset_delta_overlap(res):

    groupby_df = res.groupby(["sid", "arch"])["count_overlap"].median().reset_index()

    table = pd.pivot(groupby_df,
    index ="sid",
    columns = "arch",
    values = "count_overlap")
    print(table.head())

    table[["simple", "complex"]]= table[["simple", "complex"]].astype(int)
    fig, ax = plt.subplots(figsize = (6,6))

    x, y = "simple", "complex"
    data = table

    g = sns.pointplot(x = x, y = y, data = data, join = False)
    ax.legend(bbox_to_anchor = (1,1)).remove()
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.xaxis.set_major_locator(MultipleLocator(2))
    ax.set(xlim = (-1,14), ylim = (-1,14),
    xlabel = "median simple",
    ylabel = "median complex",
    title = "simple v. complex overla\n per tissue")
    sns.lineplot( [0,14],[0,14], ls = "--" )


def get_counts(res, strat):

    if strat == 1:
        counts = res.groupby(["mrca_2", "core_remodeling"])["enh_id"].count().reset_index()

        # add empty dataframe for complex human enhancers (do not exist by our def)
        empty = pd.DataFrame({"mrca_2": [0.000],
        "core_remodeling": [1],
        "enh_id": [0]
        })

        counts = pd.concat([counts, empty]) # concat the dataframe

        counts =counts.sort_values(by = ["core_remodeling", 'mrca_2']).reset_index()

    else:
        counts = res.groupby("arch")["enh_id"].count().reset_index()
        counts =counts.sort_values(by = "arch", ascending = False).reset_index()

    counts["enh_id"] = counts["enh_id"].astype(int)
    counts = counts.drop(["index"], axis = 1)

    return counts


def plot_figure3(res, order, fig_id, re, trim_len, frac, dataset):

    xlab = ['Homo', 'Prim', 'Euar', 'Bore', 'Euth', 'Ther', 'Mam',
    'Amni', 'Tetr', 'Vert']

    if "98" in dataset:
        title = "98 fantom Datasets"
        # for plotting, get the median overlaps per arch, per dataset

        # for plotting age stratified, only plot 10% of the data
        smallres = res.sample(frac = 0.05)
        sum_med = smallres
    else:
        title = dataset
        sum_med = res
        smallres = res

    # set up plot
    sns.set("poster")
    fig = plt.figure(figsize = (12, 8))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 3])
    ax = plt.subplot(gs[0])

    # plot the first panel
    x, y = "arch", "count_overlap"
    data = sum_med

    splot = sns.barplot( x = x, y = y, data = data,
                palette = palette, order = order,
                ax = ax)

    ax.yaxis.set_major_locator(MultipleLocator(4))

    # get counts for annotation
    STRAT = 0
    counts = get_counts(res, STRAT)
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

    x, y = "mrca_2", "count_overlap"
    data = smallres.sort_values(by = "mrca_2")
    hue = "arch"

    STRAT = 1
    agecounts = get_counts(res, STRAT)

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

    ax.set(xlabel = "", ylabel = "Number of Tissue/Cell Datasets", ylim = ax2lim)
    ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
    sns.set("poster")

    outf = f"{re}fig{fig_id}_pleiotropy-fantom_trim_{trim_len}_noexon_{dataset}_{frac}_mrcas.pdf"

    plt.savefig(outf, bbox_inches = "tight")


#%%

syn_gen_bkgd = load_syn_gen_bkgd(BUILD)
desc_df = load_desc_df()

# make dictionaries of sid + multiintersect or age files
multi, age = get_sample_dicts(multi_fs, age_fs)

#%%

# check dictionary size. Should both be 98
print(len(age.keys()), len(multi.keys()))

#%%
res, res_stats = format_df_get_stats(multi, age, RE, TRIM_LEN)


# mean of the medians of 98 tissues
print(res.groupby("arch")["count_overlap"].mean())
"""
arch        mean of median tissue overlaps in 98 datasets
complex    10.732765
simple     8.435433
p = 0.0
"""
complex = res.loc[(res.arch == "complex")]["count_overlap"]
simple = res.loc[(res.arch == "simple")]["count_overlap"]

tstat, p  = stats.mannwhitneyu(complex, simple)
print(p) #1.2407075684435812e-10
res.head()

#%% plot per dataset tissue overlap stratified by architecture



#%% plot deltas of complex v. simple

plot_per_dataset_delta_overlap(res)


#%% plot figure 3

DATASET = "all_fantom_enh"
plot_figure3(res, order, FIG_ID, RE, TRIM_LEN, FRAC, DATASET)
