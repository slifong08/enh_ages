or_ageimport argparse
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
'''
arg_parser = argparse.ArgumentParser(description=" describe argparse")

arg_parser.add_argument("enhancers", help='bed file w/ full path')
arg_parser.add_argument("shuffles", help='bed file w/ full path')
arg_parser.add_argument("genome_build", help='hg19 or hg38?')

args = arg_parser.parse_args()

ENHF = args.enhancers
SHUFF = args.shuffles
GENOME_BUILD = args.genome_build

RE = ENHF.split("data/")[0]
'''
#%%
GENOME_BUILD = "hg38"

ENHPATH = "/dors/capra_lab/users/fongsl/tyler/data/"
ENHFILE = "GG-LL_species-specific_OCRs_rank/breaks/GG-LL_species-specific_OCRs_rank_enh_age_arch_summary_matrix_W_SPECIES.bed"

ENHF = os.path.join(ENHPATH, ENHFILE)

# real shuffles
SHUFPATH = os.path.join(ENHPATH, "GG-LL_species-specific_OCRs_rank/shuffle/breaks")
SHUFFILE = "shuf-GG-LL_species-specific_OCRs_rank-0_enh_age_arch_summary_matrix.bed"

# all other OCRs
SHUFPATH = os.path.join(ENHPATH, "GG-LL_subtract_species_specific_OCR/breaks")
SHUFFILE = "GG-LL_subtract_species_specific_OCR_ages_enh_age_arch_summary_matrix.bed"

SHUFF = os.path.join(SHUFPATH, SHUFFILE)

# make a dir to save results

RE = os.path.join(ENHF.split("data/")[0], "results/")
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


def open_df(F):

    if "shuf" in F:
        cols = [0,1,2,3,4,5,6,7,8]
        col_names = ["#chr_enh", "start_enh", "end_enh", "enh_id", "sample_id",
        "seg_index", "core_remodeling", "arch", "mrca"]

        df = pd.read_csv(F,
        sep = '\t',
        header = None,
        usecols = cols,
        names = col_names,
        ).drop_duplicates()


    elif "subtract" in F:
        cols  = ["#chr_enh", "start_enh", "end_enh", "enh_id", "sample_id",
        "seg_index", "core_remodeling", "arch", "mrca"]

        df = pd.read_csv(F,
        sep = '\t',
        ).drop_duplicates()

        df = df[cols]

    else:
        cols = ["#chr_enh", "start_enh", "end_enh", "enh_id", "sample_id",
        "seg_index", "core_remodeling", "arch", "mrca", "species_specific"]

        df = pd.read_csv(F,
        sep = '\t',
        ).drop_duplicates()

        df = df[cols]

    df[["start_enh", "end_enh"]]=df[["start_enh", "end_enh"]].astype(int)

    df["enh_len"] = df.end_enh - df.start_enh
    df[["seg_index", "core_remodeling"]] = df[["seg_index", "core_remodeling"]].astype(int)

    SYN_GROUP = "/dors/capra_lab/projects/enhancer_ages/hg38_syn_taxon.bed"
    syn = pd.read_csv(SYN_GROUP, sep = '\t')

    # round all values
    syn[["mrca", "mrca_2"]] = syn[["mrca", "mrca_2"]].round(3)
    df.mrca = df.mrca.round(3)

    df = pd.merge(df, syn, how = "left")

    if "subtract" in F:
        df['id'] = "subtract"
        df["dataset_name"] = "subtract"
    else:
        df['id'] = "species_specific"
        df["dataset_name"] = "species_specific"

    return df


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

def MRCA_frequency(catdf, cols, var):

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
    enhdf = ages.loc[ages["id"] == "species_specific"]

    shuf_ = ages.loc[ages["id"] == "subtract"]

    merge_cols = list(set(cols) - set(["id"]))

    fc = pd.merge(shuf_, enhdf, how = "left", on =merge_cols)

    # calculate fold changes
    fc["fold_change"] = fc["freq_y"].divide(fc["freq_x"])


    col_id = "_".join(cols)
    outf = f'{RE}{col_id}_freq_ss_v_s.txt'
    ages.to_csv(outf, sep = '\t', index = False)

    outf = f'{RE}{col_id}_fold_change_ss_v_s.txt'
    fc.to_csv(outf, sep = '\t', index = False)

    outf = f'{RE}summary_{col_id}_freq_ss_v_s.txt'
    summarized_freq.to_csv(outf, sep = '\t', index = False)



    return ages, fc

def plot_arch_freq(age_arch_freq, age_freq):
    plots = {"age_arch" : age_arch_freq, "age": age_freq}

    for name, frame in plots.items():

        if name == "age_arch": # arrange order and colors of plot.
            frame["plot_hue"] = frame["arch"].astype(str) + "-" + frame["id"].astype(str)
            order = ["simple-species_specific", "simple-subtract",
            "complexenh-species_specific", "complexenh-subtract"]
            hue = "plot_hue"

        else:
            order = ["species_specific", "subtract"]
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

        outf = f"{RE}{name}_freq_per_age_ss_v_s.pdf"

        plt.savefig(outf, bbox_inches= "tight")


def plot_arch_fc(age_arch_fc, age_fc):

    plots = {"age_arch":age_arch_fc, "age": age_fc}

    for name, fc in plots.items():

        fc['log2'] = np.log2(fc["fold_change"])
        print(list(fc))
        if name == "age_arch":
            order = ["simple", "complexenh"]
            hue = "arch"

        else:
            order = ["species_specific"]
            hue = "id_y"


        if GENOME_BUILD == "hg38":
            xlabs = ["Prim", "Euar", "Bore", "Euth", "Ther", "Mam", "Amni", "Tetr", "Sarg", "Vert"]
        else:
            xlabs = ["Homo", "Prim", "Euar", "Bore", "Euth", "Ther", "Mam", "Amni", "Tetr", "Vert"]

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

        outf = f"{RE}{name}_fold_change_per_age_ss_v_s.pdf"

    plt.savefig(outf, bbox_inches= "tight")


def plot_len(catdf):
    plots = {"age_arch" : catdf, "age": catdf}

    for name, frame in plots.items():

        if name == "age_arch": # arrange order and colors of plot.
            frame["plot_hue"] = frame["arch"].astype(str) + "-" + frame["id"].astype(str)
            order = ["simple-species_specific", "simple-subtract",
            "complexenh-species_specific", "complexenh-subtract"]
            hue = "plot_hue"

        else:
            order = ["species_specific", "subtract"]
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

        outf = f"{RE}{name}_length_per_age_ss_v_s.pdf"

        plt.savefig(outf, bbox_inches= "tight")


def plot_cdf(catdf):


    enh_cdf = catdf.loc[catdf["id"] == "species_specific"]["seg_index"].reset_index()
    shuf_cdf = catdf.loc[catdf["id"] != "species_specific"]["seg_index"].reset_index()

    enh_cdf["pct"] = enh_cdf['seg_index'].rank(pct = True)
    shuf_cdf["pct"] = shuf_cdf['seg_index'].rank(pct = True)


    fig, ax = plt.subplots(figsize = (6,6))
    x = "seg_index"
    y = "pct"
    data = enh_cdf
    sns.lineplot(x = x, y=y, data = data, color = "blue", label = "species_specific")

    data = shuf_cdf
    sns.lineplot(x = x, y=y, data = data, color = "grey", label = "subtract")
    ax.set(xlabel = "number of age segments", ylabel = "cdf")

    outf = f"{RE}age_segment_cdf_ss_v_s.pdf"

    plt.savefig(outf, bbox_inches= "tight")


def fdr_correction(collection_dict):

    df = pd.concat(collection_dict.values())

    pvals = df["P"]

    df["rejected"], df["FDR_P"] = statsmodels.stats.multitest.fdrcorrection(pvals, alpha=0.05)
    return df


def or_seg(catdf):

    seg_dict = {} # collect results

    for seg_index in catdf.seg_index.unique():
        seg_enh = len(catdf.loc[(catdf.seg_index == seg_index) & (catdf["id"]=="species_specific")])
        not_seg_enh = len(catdf.loc[(catdf.seg_index != seg_index) & (catdf["id"]=="species_specific")])
        seg_shuf = len(catdf.loc[(catdf.seg_index == seg_index) & (catdf["id"]=="subtract")])
        not_seg_shuf = len(catdf.loc[(catdf.seg_index != seg_index) & (catdf["id"]=="subtract")])

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


    outf = f'{RE}summary_age_seg_OR_ss_v_s.txt'
    ordf.to_csv(outf, sep = '\t', index = False)

    return ordf.sort_values(by = "seg_index")


def plot_or_seg(ordf):

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

    outf = f"{RE}age_segment_or_ss_v_s.pdf"

    plt.savefig(outf, bbox_inches= "tight")

    print(ordf[["seg_index", "OR", "FDR_P", "rejected"]].sort_values(by = "seg_index"))


def or_age_arch(catdf, arch):

    mrca_dict ={}

    enh = catdf.loc[catdf["id"] == "species_specific"]
    shuffle = catdf.loc[catdf["id"] != "species_specific"]

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

    outf = f'{RE}summary_arch_age_OR_ss_v_s.txt'
    or_age_arch.to_csv(outf, sep = '\t', index = False)


    return or_age_arch


def plot_fet_age(or_age_arch, arch):

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
     title = f"enrichment per number of age segments\nspecies_specific only", xlabel = "")#ylim = (-1.2,0.5))

    plt.axhline(0, color = "grey", linewidth = 2.5)

    ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(round(2**x, 1)))

    ax.yaxis.set_major_formatter(ticks)
    ax.yaxis.set_major_locator(MultipleLocator(1))

    if GENOME_BUILD == "hg38":
        xlabs = ["Prim", "Euar", "Bore", "Euth", "Ther", "Mam", "Amni", "Tetr", "Sarg", "Vert"]
    else:
        xlabs = ["Homo", "Prim", "Euar", "Bore", "Euth", "Ther", "Mam", "Amni", "Tetr", "Vert"]

    ax.set_xticklabels(xlabs, rotation = 90)


    plt.savefig(f"{RE}{arch}_odds_per_mrca_ss_v_s.pdf", bbox_inches = 'tight')

    print(or_age_arch[["mrca_2", "OR", "FDR_P", "rejected"]].sort_values(by = "mrca_2"))


def or_age(catdf):

    mrca_dict ={}

    enh = catdf.loc[catdf["id"] == "species_specific"]
    shuffle = catdf.loc[catdf["id"] != "species_specific"]

    for mrca_2 in enh.mrca_2.unique():

        ANALYSIS = f"{mrca_2} enrichment - species-specific v. background"

        # get counts
        in_arch = enh.loc[enh.mrca_2==mrca_2, "enh_id"].count()
        not_in_arch = enh.loc[enh.mrca_2!=mrca_2, "enh_id"].count()
        shuf_in_arch = shuffle.loc[shuffle.mrca_2==mrca_2, "enh_id"].count()
        shuf_not_in_arch = shuffle.loc[shuffle.mrca_2!=mrca_2, "enh_id"].count()

        # assign 2x2
        a, b, c, d = in_arch, not_in_arch, shuf_in_arch, shuf_not_in_arch

        obs = [[a,b],[c,d]]

        OR, P = stats.fisher_exact(obs)
        table = sm.stats.Table2x2(obs) # get confidence interval
        odds_ci = table.oddsratio_confint()
        newdf = pd.DataFrame({"mrca_2":[mrca_2], "a":[a], "b":[b], "c":[c], "d":[d],
                             "OR":[OR], "P":[P], "ci_lower" :[odds_ci[0]],
                            "ci_upper" :[odds_ci[1]], "analysis":[ANALYSIS]})

        mrca_dict[mrca_2] = newdf

    or_age = fdr_correction(mrca_dict)

    outf = f'{RE}summary_age_OR_ss_v_s.txt'
    or_age.to_csv(outf, sep = '\t', index = False)


    return or_age


def plot_or_age(or_age):

    # format dataframe
    or_age["log2"] = np.log2(or_age["OR"])
    or_age["yerr"] = or_age["ci_upper"] - or_age["ci_lower"]
    or_age.sort_values(by = "mrca_2")


    fig, ax = plt.subplots(figsize = (6,6))
    sns.set("poster")
    sns.set_style("white")

    x = "mrca_2"
    y = "log2"

    data = or_age.loc[or_age.mrca_2>0].sort_values(by = "mrca_2")

    sns.barplot(
    x=x, y=y, data = data,
    linewidth=2.5, facecolor=(1, 1, 1, 0), edgecolor=".2",
    yerr =data["yerr"])



    ax.set(ylabel= f"OR species-specific v. Bkgd\n(log2-scaled)",\
     title = f"age enrichment", xlabel = "fdr<5%")

    plt.axhline(0, color = "grey", linewidth = 2.5)

    ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(round(2**x, 1)))

    ax.yaxis.set_major_formatter(ticks)
    ax.yaxis.set_major_locator(MultipleLocator(0.2))

    if GENOME_BUILD == "hg38":
        xlabs = ["Prim", "Euar", "Bore", "Euth", "Ther", "Mam", "Amni", "Tetr", "Sarg", "Vert"]
    else:
        xlabs = ["Homo", "Prim", "Euar", "Bore", "Euth", "Ther", "Mam", "Amni", "Tetr", "Vert"]

    ax.set_xticklabels(xlabs, rotation = 90)


    plt.savefig(f"{RE}OR_age_ss_v_s.pdf", bbox_inches = 'tight')

    print(or_age[["mrca_2", "OR", "FDR_P", "rejected"]].sort_values(by = "mrca_2"))

#%% run analysis

df = open_df(ENHF)

shufdf = open_df(SHUFF)


#%% concatenate enh and shuffle


catdf = pd.concat([df, shufdf])
catdf[["taxon", 'mrca']].drop_duplicates().sort_values(by = "mrca")
catdf.head()

# basic info

get_percent_simple(catdf)

count_arch = catdf.groupby(["id", "arch"])["enh_id"].count().reset_index()
count_arch
#%%


# age architecture frequencies and fold changes
cols = ["id", "arch", "mrca_2"]
var = "mrca_2"
age_arch_freq, age_archfc = MRCA_frequency(catdf, cols, var)

age_arch_freq.mrca_2.unique()
#age_arch_freq, age_archfc = arch_frequency(catdf)
#age_freq, agefc = age_frequency(catdf)
#%%

cols = ["id", "mrca_2"]
var = "mrca_2"

age_freq, agefc = MRCA_frequency(catdf, cols, var)
#%%
age_freq

plot_arch_freq(age_arch_freq, age_freq)


#%%
plot_arch_fc(age_archfc, agefc)
#%%
# architecture lengths by ages
plot_len(catdf)

#%%
# cumulative distribution of age segments across enhancers
plot_cdf(catdf)

#%%
# odds of observing segments
ordf = or_seg(catdf)

plot_or_seg(ordf)

#%% odds of observing archictures per age
ARCH = "complexenh"
or_mrca_arch = or_age_arch(catdf, ARCH)
plot_fet_age(or_mrca_arch, ARCH)
#%%
catdf["species_specific"].unique()
shared = catdf.loc[catdf["species_specific"].isna()]
shared.groupby("arch")["enh_id"].count()
hu_gain = catdf.loc[catdf["species_specific"] == "GM12878_specific"]


#%%
analysis={}

ANALYSIS = "HU-specific gains v. subtracted_background"
gains = hu_gain.groupby("arch")["enh_id"].count().reset_index()
shares = shared.groupby("arch")["enh_id"].count().reset_index()

a, b = gains.iloc[0,1], gains.iloc[1,1]
c , d = shares.iloc[0,1], shares.iloc[1,1]
obs = [[a, b], [c, d]]

OR, P = stats.fisher_exact(obs)
table = sm.stats.Table2x2(obs) # get confidence interval
odds_ci = table.oddsratio_confint()

gaindf = pd.DataFrame({"a":[a], "b":[b], "c":[c], "d":[d],
                     "OR":[OR], "P":[P], "ci_lower" :[odds_ci[0]],
                    "ci_upper" :[odds_ci[1]], "comparison": [ANALYSIS],
                    })

analysis[ANALYSIS] = gaindf
#%%
ANALYSIS = "HU-specific losses v. subtracted_background"
hu_loss = catdf.loc[catdf["species_specific"] == "LCL8664_specific"]
losses = hu_loss.groupby("arch")["enh_id"].count().reset_index()
losses
a, b = losses.iloc[0,1], losses.iloc[1,1]
obs = [[a, b], [c, d]]

OR, P = stats.fisher_exact(obs)
table = sm.stats.Table2x2(obs) # get confidence interval
odds_ci = table.oddsratio_confint()

lossdf = pd.DataFrame({"a":[a], "b":[b], "c":[c], "d":[d],
                     "OR":[OR], "P":[P], "ci_lower" :[odds_ci[0]],
                    "ci_upper" :[odds_ci[1]], "comparison": [ANALYSIS],
                    })

analysis[ANALYSIS] = lossdf

#%%
or_age = or_age(catdf)
or_age.sort_values(by = "mrca_2")
plot_or_age(or_age)
