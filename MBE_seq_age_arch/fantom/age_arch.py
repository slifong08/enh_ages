import argparse
import glob
from joblib import Parallel, delayed

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.ticker as ticker
import multiprocessing
import numpy as np
import os, sys
import pandas as pd
import subprocess
import seaborn as sns
from scipy import stats
import statsmodels
import statsmodels.api as sm


#%%
PLOT = 0
SAMPLE_ID = "all_fantom_enh"
GENOME_BUILD = "hg19"
FIG_ID = "S2"

ENHPATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/non-genic/"
ENHFILES = glob.glob(f"{ENHPATH}no-exon_all_fantom_enh/breaks/no-exon_all_fantom_enh_ages_enh_age_arch_summary_matrix.bed")
#ENHFILES = glob.glob(f"{ENHPATH}non-genic/combined_enh_shuf_abs_simple/no-exon_E*.bed")


# real shuffles

SHUFFILES = glob.glob(f"{ENHPATH}no-exon_shuffle_fantom_age_arch_summary_matrix.bed")


# make a dir to save results

RE = os.path.join(ENHFILES[0].split("data/")[0], "results/", SAMPLE_ID +"/")
if os.path.exists(RE) == False:
    os.mkdir(RE)

RE
#%% palettes


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


#%% Functions

def get_relative_simple(F, sample_id):
    # truncated columns
    col_names = ["enh_id", "seg_index"]

    df = pd.read_csv(F,
    sep = '\t',
    usecols = [3,5],
    ).drop_duplicates()

    relative_simple = df.seg_index.median()

    return relative_simple

def open_df(F, sid, relative_simple):

    col_names = ["#chr_enh", "start_enh", "end_enh", "enh_id",
     "sample_id", "core_remodeling", "arch",  "seg_index",
     "mrca", "taxon", "mrca_2", "taxon2"]

    if "shuf" in F:
        checkdf =  pd.read_csv(F,
            sep = '\t',
            nrows = 5
            )

        if "#chr_enh" in list(checkdf):
            df = pd.read_csv(F,
            sep = '\t',
            low_memory = False)
        else:
            # truncated columns
            col_names = ["#chr_enh", "start_enh", "end_enh", "enh_id", "sample_id",
            "seg_index", "core_remodeling", "arch", "mrca"]

            col_range = np.arange(0,9)


            df = pd.read_csv(F,
            sep = '\t',
            usecols = col_range,
            header = None,
            low_memory = False
            ).drop_duplicates()

            df.columns = col_names
            df.to_csv(F, sep = '\t', index = None)


    else:

        df = pd.read_csv(F,
        sep = '\t',
        ).drop_duplicates()

        #df.columns = col_names

    # assign relative simple architecture

    print("relative simple # of age segments <", relative_simple)

    df["arch"] = "simple"
    if relative_simple >1:
        df.loc[df.seg_index >= relative_simple, "arch"] = "complexenh"
    elif relative_simple == 1:
        df.loc[df.seg_index > relative_simple, "arch"] = "complexenh"

    # add back core_remodeling annotation
    df["core_remodeling"] = 0
    df.loc[df.arch == "complexenh", "core_remodeling"] = 1

    #print(df.head())
    df.to_csv(F, sep = '\t', index = None)

    # do the rest of the formatting.
    df["sid"] = sid
    df[["start_enh", "end_enh"]] = df[["start_enh", "end_enh"]].astype(int)

    df["enh_len"] = df.end_enh - df.start_enh
    df[["seg_index", "core_remodeling"]] = df[["seg_index", "core_remodeling"]].astype(int)

    SYN_GROUP = f"/dors/capra_lab/projects/enhancer_ages/{GENOME_BUILD}_syn_taxon.bed"
    syn = pd.read_csv(SYN_GROUP, sep = '\t')

    # round all values
    syn[["mrca", "mrca_2"]] = syn[["mrca", "mrca_2"]].round(3)
    df.mrca = df.mrca.round(3)

    df = pd.merge(df, syn, how = "left")

    if "shuf" in F:
        df['id'] = "shuf"
    else:
        df['id'] = "enh"

    dataset_name = (F.split("/")[-1]).split(".")[0]
    print(dataset_name)
    df['dataset_name'] =dataset_name

    return df, dataset_name

def MRCA_frequency(catdf, cols, var, sample_id):

    age_dict = {} # collect age frequency results per dataset
    summary_age_dict = {} # collect summarized age frequencies per dataset

    for n, dataset in enumerate(catdf.dataset_name.unique()):
        # count n enhancers in architecture per age
        test = catdf.loc[catdf.dataset_name == dataset]

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


        age_dict[n] = age

        # summarize frequencies across architectures, before/after eutherian.
        euth_age = catdf.loc[catdf.taxon2 == "Eutheria (105)", "mrca_2"].iloc[0]
        eutherian = age.loc[age["mrca_2"] == euth_age][[ "id", "freq"]]
        eutherian["category"] = "eutherian"

        younger_thaneuth = age.loc[age["mrca_2"] <euth_age].groupby(["id"])["freq"].sum().reset_index()
        younger_thaneuth["category"] = "younger than eutherian"

        older_thaneuth = age.loc[age["mrca_2"] >euth_age].groupby(["id"])["freq"].sum().reset_index()
        older_thaneuth["category"] = "older than eutherian"

        summarized_freq = pd.concat([eutherian, younger_thaneuth, older_thaneuth])


        summary_age_dict[n] = summarized_freq

    # concat age and summarized frequency dataframes
    ages = pd.concat(age_dict.values())
    summarized_freq = pd.concat(summary_age_dict.values())

    # calculate fold-change of enh v. shuf expectation per shuffle


    # select only the enhancer and specific shuffle instance
    enhdf = ages.loc[ages["id"] == "enh"]

    shuf_ = ages.loc[ages["id"] != "enh"]

    merge_cols = list(set(cols) - set(["id"]))

    fc = pd.merge(shuf_, enhdf, how = "left", on =merge_cols)

    # calculate fold changes
    fc["fold_change"] = fc["freq_y"].divide(fc["freq_x"])


    col_id = "_".join(cols)
    outf = f'{RE}{sample_id}_{col_id}_freq.txt'
    ages.to_csv(outf, sep = '\t', index = False)

    outf = f'{RE}{sample_id}_{col_id}_fold_change.txt'
    fc.to_csv(outf, sep = '\t', index = False)

    outf = f'{RE}{sample_id}_summary_{col_id}_freq.txt'
    summarized_freq.to_csv(outf, sep = '\t', index = False)

    return ages, fc

def get_percent_simple(catdf, sample_id):


    count_arch = catdf.groupby(["id", "arch"])["enh_id"].count().reset_index()

    #dELS simple = 58.6%
    simple, complex = count_arch.iloc[1,2], count_arch.iloc[0,2]
    simple_shuf, complex_shuf = count_arch.iloc[3,2], count_arch.iloc[2,2]

    per_simple_enh = simple/(simple+complex)

    #5x shuffle = 59.1%
    per_simple_shuf = (simple_shuf/(simple_shuf+complex_shuf))

    print("% simple enh =", per_simple_enh, "% simple shuf =", per_simple_shuf)

    # quantify by architecgture
    vals = ["enh_len", "seg_index", "mrca_2"]
    val_dict = {}
    for val in vals:
        lens = catdf.groupby(["id", "arch"])[val].describe().reset_index()
        lens["val"] = val
        val_dict[val] = lens


    measures = pd.concat(val_dict.values())
    outf = f"{RE}{sample_id}_ARCH_metrics.txt"
    measures.to_csv(outf, sep = '\t')

    # quantify by age and architecture
    vals = ["enh_len", "seg_index"]
    val_dict = {}
    for val in vals:
        lens = catdf.groupby(["id", "arch", "mrca_2"])[val].describe().reset_index()
        lens["val"] = val
        val_dict[val] = lens
    outf = f"{RE}{sample_id}_AGE_ARCH_metrics.txt"
    measures = pd.concat(val_dict.values())
    measures.to_csv(outf, sep = '\t')

def fdr_correction(collection_dict):

    df = pd.concat(collection_dict.values())

    pvals = df["P"]

    df["rejected"], df["FDR_P"] = statsmodels.stats.multitest.fdrcorrection(pvals, alpha=0.05)
    return df

def or_seg(catdf, sample_id):

    seg_dict = {} # collect results

    for seg_index in catdf.seg_index.unique():
        seg_enh = len(catdf.loc[(catdf.seg_index == seg_index) & (catdf["id"]=="enh")])
        not_seg_enh = len(catdf.loc[(catdf.seg_index != seg_index) & (catdf["id"]=="enh")])
        seg_shuf = len(catdf.loc[(catdf.seg_index == seg_index) & (catdf["id"]=="shuf")])
        not_seg_shuf = len(catdf.loc[(catdf.seg_index != seg_index) & (catdf["id"]=="shuf")])

        a, b, c, d = seg_enh, not_seg_enh, seg_shuf,not_seg_shuf
        obs = [[a,b], [c,d]]

        OR, P = stats.fisher_exact(obs)
        table = sm.stats.Table2x2(obs) # get confidence interval
        odds_ci = table.oddsratio_confint()

        newdf = pd.DataFrame({"seg_index":[seg_index], "a":[a], "b":[b], "c":[c], "d":[d],
                             "OR":[OR], "P":[P], "ci_lower" :[odds_ci[0]],
                            "ci_upper" :[odds_ci[1]]})

        seg_dict[seg_index] = newdf

    ordf = fdr_correction(seg_dict)


    outf = f'{RE}{sample_id}_summary_seg_OR.txt'
    ordf.to_csv(outf, sep = '\t', index = False)

    return ordf.sort_values(by = "seg_index")

def or_age_arch(catdf, arch, sample_id):

    mrca_dict ={}

    enh = catdf.loc[catdf["id"] == "enh"]
    shuffle = catdf.loc[catdf["id"] == "shuf"]

    for mrca_2 in enh.mrca_2.unique():

        # subset dataframes
        in_age_enh = enh.loc[enh.mrca_2 == mrca_2]
        in_age_shuf = shuffle.loc[shuffle.mrca_2 == mrca_2]


        # get counts
        in_arch = len(in_age_enh.loc[in_age_enh.arch==arch])
        not_in_arch = len(in_age_enh.loc[in_age_enh.arch!=arch])
        shuf_in_arch = len(in_age_shuf.loc[in_age_shuf.arch==arch])
        shuf_not_in_arch = len(in_age_shuf.loc[in_age_shuf.arch!=arch])

        # assign 2x2
        a, b, c, d = in_arch, not_in_arch, shuf_in_arch, shuf_not_in_arch

        obs = [[a,b],[c,d]]

        OR, P = stats.fisher_exact(obs)
        table = sm.stats.Table2x2(obs) # get confidence interval
        odds_ci = table.oddsratio_confint()
        newdf = pd.DataFrame({"mrca_2":[mrca_2], "a":[a], "b":[b], "c":[c], "d":[d],
                             "OR":[OR], "P":[P], "ci_lower" :[odds_ci[0]],
                            "ci_upper" :[odds_ci[1]], "core_remodeling_a":[arch]})


        mrca_dict[mrca_2] = newdf

    or_age_arch = fdr_correction(mrca_dict)

    outf = f'{RE}{sample_id}_summary_age_arch_OR.txt'
    or_age_arch.to_csv(outf, sep = '\t', index = False)


    return or_age_arch

def plot_arch_freq(age_arch_freq, age_freq, sample_id):
    plots = {"age_arch" : age_arch_freq, "age": age_freq}

    for name, frame in plots.items():

        if name == "age_arch": # arrange order and colors of plot.
            frame["plot_hue"] = frame["arch"].astype(str) + "-" + frame["id"].astype(str)
            order = ["simple-enh", "simple-shuf",
            "complexenh-enh", "complexenh-shuf"]
            hue = "plot_hue"

        else:
            order = ["enh", "shuf"]
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
        palette = ESPAL)

        ax.set_xticklabels(xlabs, rotation = 90)
        ax.legend(bbox_to_anchor = (1,1))

        outf = f"{RE}{sample_id}{name}_freq_per_age.pdf"

        plt.savefig(outf, bbox_inches= "tight")

def plot_arch_fc(age_arch_fc, age_fc, sample_id):

    plots = {"age_arch":age_arch_fc, "age": age_fc}

    for name, fc in plots.items():

        fc['log2'] = np.log2(fc["fold_change"])

        if name == "age_arch":
            order = ["simple", "complexenh"]
            hue = "arch"

        else:
            order = ["enh"]
            hue = "id_y"


        if GENOME_BUILD == "hg38":
            xlabs = ["Prim", "Euar", "Bore", "Euth", "Ther", "Mam", "Amni", "Tetr", "Sarg", "Vert"]
        else:
            xlabs = ["Prim", "Euar", "Bore", "Euth", "Ther", "Mam", "Amni", "Tetr", "Vert"]

        sns.set("talk")
        fig, ax = plt.subplots(figsize = (6,6))
        x, y = "mrca_2", "log2"
        data = fc

        sns.barplot(x = x, y=y,
        data = data,
        hue = hue,
        hue_order = order,
        palette = EPAL)
        ax.set(ylabel = "Fold-Change v. Bkgd\n(log2-scaled)")
        ax.set_xticklabels(xlabs, rotation = 90)
        ax.legend(bbox_to_anchor = (1,1))

        outf = f"{RE}{FIG_ID}_{sample_id}_{name}_fold_change_per_age.pdf"

    plt.savefig(outf, bbox_inches= "tight")

def plot_len(catdf, sample_id):
    plots = {"age_arch" : catdf, "age": catdf}

    for name, frame in plots.items():

        if name == "age_arch": # arrange order and colors of plot.
            frame["plot_hue"] = frame["arch"].astype(str) + "-" + frame["id"].astype(str)
            order = ["simple-enh", "simple-shuf",
            "complexenh-enh", "complexenh-shuf"]
            hue = "plot_hue"

        else:
            order = ["enh", "shuf"]
            hue = "id"


        if GENOME_BUILD == "hg38":
            xlabs = ["Prim", "Euar", "Bore", "Euth", "Ther", "Mam", "Amni", "Tetr", "Sarg", "Vert"]
        else:
            xlabs = ["Homo", "Prim", "Euar", "Bore", "Euth", "Ther", "Mam", "Amni", "Tetr", "Vert"]
        sns.set("talk")
        fig, ax = plt.subplots(figsize = (6,6))
        x, y = "mrca_2", "enh_len"

        data = frame

        sns.barplot(x = x, y=y,
        data = data,
        hue = hue,
        hue_order = order,
        palette = ESPAL)

        ax.set(ylabel = "length(bp)",
        #ylim = (250, 325)
        )
        ax.set_xticklabels(xlabs, rotation = 90)
        ax.legend(bbox_to_anchor = (1,1))

        outf = f"{RE}{FIG_ID}_{sample_id}_{name}_length_per_age.pdf"

        plt.savefig(outf, bbox_inches= "tight")

def plot_cdf(catdf, sample_id):


    enh_cdf = catdf.loc[catdf["id"] == "enh"]["seg_index"].reset_index()
    shuf_cdf = catdf.loc[catdf["id"] == "shuf"]["seg_index"].reset_index()

    enh_cdf["pct"] = enh_cdf['seg_index'].rank(pct = True)
    shuf_cdf["pct"] = shuf_cdf['seg_index'].rank(pct = True)


    fig, ax = plt.subplots(figsize = (6,6))
    x = "seg_index"
    y = "pct"
    data = enh_cdf
    sns.lineplot(x = x, y=y, data = data, color = "blue", label = "enh")

    data = shuf_cdf
    sns.lineplot(x = x, y=y, data = data, color = "grey", label = "shuf")
    ax.set(xlabel = "number of age segments", ylabel = "cdf")

    outf = f"{RE}{FIG_ID}_{sample_id}_age_segment_cdf.pdf"

    plt.savefig(outf, bbox_inches= "tight")

def plot_or_seg(ordf, sample_id):

    ordf["log2"] = np.log2(ordf["OR"])
    ordf["yerr"] = ordf["ci_upper"] - ordf["ci_lower"]

    x = "seg_index"
    y = 'log2'
    data = ordf.loc[ordf.seg_index<=10].sort_values(by = "seg_index")

    fig, ax = plt.subplots(figsize = (6,6))
    sns.barplot(x = x, y =y, data = data,
    linewidth=2.5, facecolor=(1, 1, 1, 0), edgecolor=".2",
    yerr = data["yerr"])

    ax.set(ylabel= "Fold Change v. Bkgd\n(log2-scaled)",\
     title = "enrichment per age segment",
     xlabel = "number of age segments")#ylim = (-1.2,0.5))

    outf = f"{RE}{FIG_ID}_{sample_id}_age_segment_or.pdf"

    plt.savefig(outf, bbox_inches= "tight")

    print(ordf[["seg_index", "OR", "FDR_P", "rejected"]].sort_values(by = "seg_index"))

def plot_fet_age(or_age_arch, arch, sample_id):

    # format dataframe
    or_age_arch["log2"] = np.log2(or_age_arch["OR"])
    or_age_arch["yerr"] = or_age_arch["ci_upper"] - or_age_arch["ci_lower"]
    or_age_arch.sort_values(by = "mrca_2")

    fig, ax = plt.subplots(figsize = (6,6))
    sns.set("poster")
    sns.set_style("white")

    x = "mrca_2"
    y = "log2"

    data = or_age_arch.loc[or_age_arch.mrca_2>0].sort_values(by = "mrca_2")

    sns.barplot( x=x, y=y, data = data,
    linewidth=2.5, facecolor=(1, 1, 1, 0), edgecolor=".2",
    yerr =data["yerr"])

    if arch == "complexenh":
        other_arch = "simple"
    else:
        other_arch = "complexenh"

    ax.set(ylabel= f"Fold Change\n{arch} v. {other_arch}\n(log2-scaled)",\
     title = f"enrichment per number of age segments", xlabel = "")#ylim = (-1.2,0.5))

    plt.axhline(0, color = "grey", linewidth = 2.5)

    ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(round(2**x, 1)))

    ax.yaxis.set_major_formatter(ticks)
    ax.yaxis.set_major_locator(MultipleLocator(1))

    if GENOME_BUILD == "hg38":
        xlabs = ["Prim", "Euar", "Bore", "Euth", "Ther", "Mam", "Amni", "Tetr", "Sarg", "Vert"]
    else:
        xlabs = ["Prim", "Euar", "Bore", "Euth", "Ther", "Mam", "Amni", "Tetr", "Vert"]

    ax.set_xticklabels(xlabs, rotation = 90)

    outf = f"{RE}{FIG_ID}_{sample_id}_{arch}_odds_per_mrca.pdf"
    plt.savefig(outf, bbox_inches = 'tight')

    print(or_age_arch[["mrca_2", "OR", "FDR_P", "rejected"]].sort_values(by = "mrca_2"))


# run analysis


def run_analyses(enhf, shuffle_list, sample_id, plot):

    relative_simple = get_relative_simple(enhf, sample_id)
    df, sid = open_df(enhf, sample_id, relative_simple)

    shuf_dict = {} # open the shuffle df
    for n, SHUFF in enumerate(shuffle_list):
        shufdf_, dataset_name = open_df(SHUFF, sample_id, relative_simple)
        shuf_dict[n] =shufdf_

    # concatenate enh and shuffle
    shufdf = pd.concat(shuf_dict.values())
    catdf = pd.concat([df, shufdf])
    catdf[["taxon", 'mrca']].drop_duplicates().sort_values(by = "mrca")

    # basic info
    get_percent_simple(catdf, sample_id)

    count_arch = catdf.groupby(["id", "arch"])["enh_id"].count().reset_index()


    # age architecture frequencies and fold changes
    cols = ["id", "arch", "mrca_2"]
    var = "mrca_2"
    age_arch_freq, age_archfc = MRCA_frequency(catdf, cols, var, sample_id)

    # get age freq w.o. arch

    cols = ["id", "mrca_2"]
    var = "mrca_2"

    age_freq, agefc = MRCA_frequency(catdf, cols, var, sample_id)


    # odds of observing segments
    ordf = or_seg(catdf, sample_id)

    # odds of observing archictures per age
    ARCH = "complexenh"
    or_mrca_arch = or_age_arch(catdf, ARCH, sample_id)



    if plot == 1:

        plot_arch_freq(age_arch_freq, age_freq, sample_id)

        # plot the fold change v. shuffled background.

        plot_arch_fc(age_archfc, agefc, sample_id)

        # architecture lengths by ages
        plot_len(catdf, sample_id)

        # cumulative distribution of age segments across enhancers

        plot_cdf(catdf, sample_id)

        plot_or_seg(ordf, sample_id)

        plot_fet_age(or_mrca_arch, ARCH, sample_id)

    return catdf

#%% run analysis

# make dictionary of shuffle files
shuf_id_dict = {}
relative_dict = {}
for SHUFF in SHUFFILES:
    sample_id = get_sample_id(SHUFF)
    shuf_id_dict[sample_id] = SHUFF

SHUFFILES
#%%

#%%
PLOT = 0
SAMPLE_ID = "all_fantom_enh"
ENHF = ENHFILES[0]
data = run_analyses(ENHF, SHUFFILES, SAMPLE_ID, PLOT)
#%%
ENHF

#%%
FS = glob.glob("/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19/download/h3k27ac_plus_h3k4me3_minus_peaks/Hsap_H3K27ac_plus_H3K4me3_minus_E*/non-genic/trimmed/trim_310_E*/breaks/trim_310_E0*_enh_age_arch_summary_matrix.bed")
PLOT = 0
GENOME_BUILD = "hg19"
FIG_ID = "S29"

for ENHF in FS:

    if os.path.exists(RE) == False:
        os.mkdir(RE)
    SAMPLE_ID = (ENHF.split("trimmed/")[1]).split("/breaks")[0]
    print(SAMPLE_ID)
    #RE = os.path.join(ENHF[0].split("data/")[0], "results/", SAMPLE_ID +"/")

    #SAMPLE_ID = "trim_310_E016"

    ENHPATH = "/".join(ENHF.split("/")[:-2])

    # real shuffles
    SHUFPATH = os.path.join(ENHPATH, "shuffle/breaks/")
    SHUFFILES = glob.glob(f"{SHUFPATH}shuf*.bed")
    print( 'shuffles n= ', len(SHUFFILES))

    # make a dir to save results
    RE = os.path.join(ENHF.split("data/")[0], "results/stats/")
    if os.path.exists(RE) == False:
        os.mkdir(RE)

    catdf = run_analyses(ENHF, SHUFFILES, SAMPLE_ID, PLOT)
#%%
FS = glob.glob("/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19/download/h3k27ac_plus_h3k4me3_minus_peaks/Hsap_H3K27ac_plus_H3K4me3_minus_E*/non-genic/no-exon_E*/breaks/no-exon_E*_enh_age_arch_summary_matrix.bed")

FIG_ID = "S29"

for ENHF in FS:

    if os.path.exists(RE) == False:
        os.mkdir(RE)
    SAMPLE_ID = (ENHF.split("non-genic/")[1]).split("/breaks")[0]
    print(SAMPLE_ID)
    #RE = os.path.join(ENHF[0].split("data/")[0], "results/", SAMPLE_ID +"/")

    #SAMPLE_ID = "trim_310_E016"

    ENHPATH = "/".join(ENHF.split("/")[:-2])

    # real shuffles
    SHUFPATH = os.path.join(ENHPATH, "shuffle/breaks/")
    SHUFFILES = glob.glob(f"{SHUFPATH}shuf*.bed")
    print( 'shuffles n= ', len(SHUFFILES))

    # make a dir to save results
    RE = os.path.join(ENHF.split("data/")[0], "results/stats/")
    if os.path.exists(RE) == False:
        os.mkdir(RE)

    catdf = run_analyses(ENHF, SHUFFILES, SAMPLE_ID, PLOT)

RE
