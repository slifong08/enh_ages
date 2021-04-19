import argparse
import glob
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.ticker as ticker
import numpy as np
import os, sys
import pandas as pd
import subprocess
import seaborn as sns
from scipy import stats
import statsmodels
import statsmodels.api as sm


#%%
GENOME_BUILD = "hg19"
FIG_ID = "S29"

ENHPATH ="/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19/download/h3k27ac_plus_h3k4me3_minus_peaks/"
ENHFILES = glob.glob(f"{ENHPATH}Hsap_H3K27ac_plus_H3K4me3_minus_E*/non-genic/no-exon_E*_parallel_breaks_enh_age_arch_summary_matrix.bed")


# real shuffles
SHUFPATH = os.path.join(ENHPATH, "shuffle/")
SHUFFILES = glob.glob(f"{SHUFPATH}shuf-cleaned_Hsap_H3K27ac_plus_H3K4me3_minus_E*_enh_age_arch_summary_matrix.bed")


# make a dir to save results

RE = os.path.join(ENHFILES[0].split("data/")[0], "results/")
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
def get_sample_id(F):

    if 'shuf' in F:
        split1 = "shuf-cleaned_Hsap_H3K27ac_plus_H3K4me3_minus_E"
        split2 = "_enh_age_arch_summary_matrix"
        sample_id =  "E" + (SHUFF.split(split1)[1]).split(split2)[0]
    else:
        split1 = "no-exon_E"
        split2 = "_parallel_breaks_enh_age_arch_summary_matrix"
        sample_id = "E" + (ENHF.split(split1)[1]).split(split2)[0]

    return sample_id


def open_df(F, sid):

    if "shuf" in F:
        checkdf =  pd.read_csv(F,
            sep = '\t',
            nrows = 5
            )
        if "#chr_enh" in list(checkdf):
            df = pd.read_csv(F,
            sep = '\t',)
        else:
            # truncated columns
            col_names = ["#chr_enh", "start_enh", "end_enh", "enh_id", "sample_id",
            "seg_index", "core_remodeling", "arch", "mrca"]
            col_range = np.arange(0,9)


            df = pd.read_csv(F,
            sep = '\t',
            usecols = col_range,
            header = None
            ).drop_duplicates()

            df.columns = col_names
            df.to_csv(F, sep = '\t', index = None)
        print(df.head())

    else:

        # truncated columns
        col_names = ["#chr_enh", "start_enh", "end_enh", "enh_id",
         "sample_id", "seg_index", "arch", "core_remodeling", "mrca"]

        col_range = np.arange(0,10)

        df = pd.read_csv(F,
        sep = '\t',
        ).drop_duplicates()

        df.columns = col_names



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

        age["dataset_name"] = dataset_name
        age_dict[n] = age

        # summarize frequencies across architectures, before/after eutherian.
        eutherian = age.loc[age["mrca_2"] == 0.19][[ "id", "freq"]]
        eutherian["category"] = "eutherian"

        younger_thaneuth = age.loc[age["mrca_2"] <0.19].groupby(["id"])["freq"].sum().reset_index()
        younger_thaneuth["category"] = "younger than eutherian"

        older_thaneuth = age.loc[age["mrca_2"] >0.19].groupby(["id"])["freq"].sum().reset_index()
        older_thaneuth["category"] = "older than eutherian"

        summarized_freq = pd.concat([eutherian, younger_thaneuth, older_thaneuth])
        summarized_freq["dataset_name"] = dataset_name

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
    outf = f'{RE}{sample_id}{sample_id}{col_id}_freq.txt'
    ages.to_csv(outf, sep = '\t', index = False)

    outf = f'{RE}{sample_id}{sample_id}{col_id}_fold_change.txt'
    fc.to_csv(outf, sep = '\t', index = False)

    outf = f'{RE}{sample_id}{sample_id}summary_{col_id}_freq.txt'
    summarized_freq.to_csv(outf, sep = '\t', index = False)



    return ages, fc

def age_frequency(catdf, sample_id):

    age_dict = {} # collect age frequency results per dataset
    fc_dict = {} # collect fold chanfe results
    summary_age_dict = {} # collect summarized age frequencies per dataset

    for dataset in catdf.dataset_name.unique():
        # count n enhancers in architecture per age
        test = catdf.loc[catdf.dataset_name == dataset]

        age = test.groupby(["id", "mrca_2"])["enh_id"].count().reset_index()
        # rename columns
        age.columns = ["id", "mrca_2", "counts"]

        # sum total n enhancers in architecture
        totals = age.groupby(["id"])["counts"].sum().reset_index()
        # rename columns
        totals.columns = ["id", "total_id"]

        # merge dataframes
        age = pd.merge(age, totals, how = "left")

        # calculate the % of architecture in each age
        age["age_freq"] = age.counts.divide(age.total_id)

        age["dataset_name"] = dataset_name
        age_dict[dataset_name] = age

        # summarize frequencies across architectures, before/after eutherian.
        eutherian = age.loc[age["mrca_2"] == 0.19][[ "id", "age_freq"]]
        eutherian["category"] = "eutherian"

        younger_thaneuth = age.loc[age["mrca_2"] <0.19].groupby(["id"])["age_freq"].sum().reset_index()
        younger_thaneuth["category"] = "younger than eutherian"

        older_thaneuth = age.loc[age["mrca_2"] >0.19].groupby(["id"])["age_freq"].sum().reset_index()
        older_thaneuth["category"] = "older than eutherian"

        summarized_freq = pd.concat([eutherian, younger_thaneuth, older_thaneuth])
        summarized_freq["dataset_name"] = dataset_name

        summary_age_dict[dataset_name] = summarized_freq

    # concat age and summarized frequency dataframes
    age = pd.concat(age_dict.values())
    summarized_freq = pd.concat(summary_age_dict.values())

    # calculate fold-change of enh v. shuf expectation per shuffle
    for dataset in catdf.dataset_name.unique():

        test = age.loc[age["id"] == "enh" | age[dataset_name] == dataset] # select only the enhancer and specific shuffle instance
        fc = test.groupby(["mrca_2", "id"])["age_freq"].max().unstack("id").reset_index()

        fc["fold_change"] = fc["enh"].divide(fc["shuf"])

        fc_dict[dataset] = fc

    fc = pd.concat(age_dict.values())


    outf = f'{RE}{sample_id}{sample_id}age_freq.txt'
    age.to_csv(outf, sep = '\t', index = False)

    outf = f'{RE}{sample_id}{sample_id}age_fold_change.txt'
    fc.to_csv(outf, sep = '\t', index = False)

    outf = f'{RE}{sample_id}{sample_id}summary_age_freq.txt'
    summarized_freq.to_csv(outf, sep = '\t', index = False)

    print(summarized_freq)

    return age, fc

def arch_frequency(catdf, sample_id):

    age_arch_dict = {} # collect age frequency results per dataset
    fc_dict = {} # collect fold chanfe results
    summary_age_arch_dict = {} # collect summarized age frequencies per dataset

    for dataset_name in catdf.dataset_name.unique():

        test = catdf.loc[catdf.dataset_name == dataset_name]

        # count n enhancers in architecture per age
        age_arch = test.groupby(["id", "arch", "mrca_2"])["enh_id"].count().reset_index()
        # rename columns
        age_arch.columns = ["id", "arch", "mrca_2", "counts"]
        # sum total n enhancers in architecture
        totals = age_arch.groupby(["id", "arch"])["counts"].sum().reset_index()
        # rename columns
        totals.columns = ["id", "arch", "total_arch"]

        # merge dataframes
        age_arch = pd.merge(age_arch, totals, how = "left")

        # calculate the % of architecture in each age
        age_arch["arch_freq"] = age_arch.counts.divide(age_arch.total_arch)

    outf = f'{RE}{sample_id}age_arch_freq.txt'
    age_arch.to_csv(outf, sep = '\t', index = False)

    # calculate fold-change of enh v. shuf expectation
    fc = age_arch.groupby(["arch", "mrca_2", "id"])["arch_freq"].max().unstack("id").reset_index()

    fc["fold_change"] = fc["enh"].divide(fc["shuf"])

    outf = f'{RE}{sample_id}age_arch_fold_change.txt'
    fc.to_csv(outf, sep = '\t', index = False)

    # summarize frequencies across architectures, before/after eutherian.

    eutherian = age_arch.loc[age_arch["mrca_2"] == 0.19][["arch", "id", "arch_freq"]]
    eutherian["category"] = "eutherian"

    younger_thaneuth = age_arch.loc[age_arch["mrca_2"] <0.19].groupby(["arch", "id"])["arch_freq"].sum().reset_index()
    younger_thaneuth["category"] = "younger than eutherian"

    older_thaneuth = age_arch.loc[age_arch["mrca_2"] >0.19].groupby(["arch", "id"])["arch_freq"].sum().reset_index()
    older_thaneuth["category"] = "older than eutherian"

    summarized_freq = pd.concat([eutherian, younger_thaneuth, older_thaneuth])

    outf = f'{RE}{sample_id}summary_arch_freq.txt'
    summarized_freq.to_csv(outf, sep = '\t', index = False)

    print(summarized_freq)

    return age_arch, fc


def get_percent_simple(catdf):

    median_segs = catdf.groupby("id").seg_index.median()
    print("median number of age segments =", median_segs, "\n")

    count_arch = catdf.groupby(["id", "arch"])["enh_id"].count().reset_index()

    #dELS simple = 58.6%
    simple, complex = count_arch.iloc[1,2], count_arch.iloc[0,2]
    simple_shuf, complex_shuf = count_arch.iloc[3,2], count_arch.iloc[2,2]

    per_simple_enh = simple/(simple+complex)

    #5x shuffle = 59.1%
    per_simple_shuf = (simple_shuf/(simple_shuf+complex_shuf))

    print("% simple enh =", per_simple_enh, "% simple shuf =", per_simple_shuf)

    print(catdf.groupby(["id", "arch"])["enh_len"].describe())


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


    outf = f'{RE}{FIG_ID}_{sample_id}_summary_age_seg_OR.txt'
    ordf.to_csv(outf, sep = '\t', index = False)

    return ordf.sort_values(by = "seg_index")


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

    outf = f'{RE}{FIG_ID}_{sample_id}_summary_arch_age_OR.txt'
    or_age_arch.to_csv(outf, sep = '\t', index = False)


    return or_age_arch


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
        xlabs = ["Homo", "Prim", "Euar", "Bore", "Euth", "Ther", "Mam", "Amni", "Tetr", "Vert"]

    ax.set_xticklabels(xlabs, rotation = 90)


    plt.savefig(f"{RE}{FIG_ID}_{sample_id}_{arch}_odds_per_mrca.pdf", bbox_inches = 'tight')

    print(or_age_arch[["mrca_2", "OR", "FDR_P", "rejected"]].sort_values(by = "mrca_2"))


#%% run analysis

#%% make dictionary of shuffle files
shuf_id_dict = {}

for SHUFF in SHUFFILES:
    sample_id = get_sample_id(SHUFF)
    shuf_id_dict[sample_id] = SHUFF

for ENHF in ENHFILES:

    sample_id = get_sample_id(ENHF)
    SHUFF = shuf_id_dict[sample_id]
    FIG_ID = f"S29-{sample_id}"
    print(sample_id, FIG_ID)

    df, sid = open_df(ENHF, sample_id)

    shufdf_, dataset_name = open_df(SHUFF, sample_id)
    shufdf = shufdf_.sample(frac = 0.01)
    # concatenate enh and shuffle


    catdf = pd.concat([df, shufdf])
    catdf[["taxon", 'mrca']].drop_duplicates().sort_values(by = "mrca")

    # basic info
    get_percent_simple(catdf)

    count_arch = catdf.groupby(["id", "arch"])["enh_id"].count().reset_index()

    #

    # age architecture frequencies and fold changes
    cols = ["id", "arch", "mrca_2"]
    var = "mrca_2"
    age_arch_freq, age_archfc = MRCA_frequency(catdf, cols, var, sample_id)

    #

    cols = ["id", "mrca_2"]
    var = "mrca_2"

    age_freq, agefc = MRCA_frequency(catdf, cols, var, sample_id)

    #

    plot_arch_freq(age_arch_freq, age_freq, sample_id)

    #

    plot_arch_fc(age_archfc, agefc, sample_id)

    # architecture lengths by ages
    plot_len(catdf, sample_id)

    # cumulative distribution of age segments across enhancers

    plot_cdf(catdf, sample_id)

    # odds of observing segments
    ordf = or_seg(catdf, sample_id)

    plot_or_seg(ordf, sample_id)

    # odds of observing archictures per age

    ARCH = "complexenh"
    or_mrca_arch = or_age_arch(catdf, ARCH, sample_id)
    plot_fet_age(or_mrca_arch, ARCH, sample_id)
    break
#%%
