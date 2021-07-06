import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy import stats
import seaborn as sns
import statsmodels
import statsmodels.api as sm
import subprocess

CL = "K562"
GENOME_BUILD = "hg38"

cCREPATH = f"/dors/capra_lab/projects/enhancer_ages/encode/data/ELS_combined_{CL}/ages/"
cCREFILE = f"syn_breaks_ELS_combined_{CL}_ages.bed"
cCRE = os.path.join(cCREPATH, cCREFILE)

#cCRE_TFBS_ONLY = f"{cCREPATH}enh_tfbs_only.txt"

SHUFPATH = f"/dors/capra_lab/projects/enhancer_ages/encode/data/ELS_combined_{CL}/shuffle/ages/"
SHUFFILE = f"syn_breaks_shuf-ELS_combined_{CL}.bed"
SHUFF = os.path.join(SHUFPATH, SHUFFILE)


RE = f"/dors/capra_lab/projects/enhancer_ages/landscape/results/cCRE/{CL}/"
SHUFF
#%%
colors = [ "amber", "faded green"]
yg = sns.xkcd_palette(colors)
sns.palplot(yg)

colors = [ "dusty blue", "greyish"]
bgy = sns.xkcd_palette(colors)
sns.palplot(bgy)

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
blue = sns.xkcd_palette(colors)
sns.palplot(blue)

colors = [ "dusty purple"]
purple = sns.xkcd_palette(colors)
sns.palplot(purple)

colors = [ "amber"]
amber = sns.xkcd_palette(colors)
sns.palplot(amber)

colors = ["amber", "greyish", "faded green", "slate grey"]
ESPAL = sns.xkcd_palette(colors)
sns.palplot(ESPAL)

colors = ["faded green",  "greyish"]
EPAL = sns.xkcd_palette(colors)
sns.palplot(EPAL)


#%%


def add_arch_labels(df):

    df["arch"] = "complex_core"

    df.loc[df.core_remodeling ==0, "arch"] = "simple"
    df.loc[df.core ==0, "arch"] = "complex_derived"

    return df
def reEval_PrimComplex(enh):

    # get all the complex enhancers w/ primate core ages
    prComEnhID = enh.loc[(enh.core ==1) &
    (enh.core_remodeling ==1) &
    (enh.taxon2.str.contains("Primate"))]["enh_id"].unique()

    # get all the complex enhancer ids where there is a real human derived sequence
    pr_complex = enh.loc[(enh.enh_id.isin(prComEnhID)) &
    (enh.core_remodeling == 1) &
    (enh.core ==0) &
    (enh.mrca ==0),
    ]["enh_id"]


    # i'm going to reassign any primate complex enhancer
    # where derived regions are from other primates
    # get the set of primate complex enhancers w/ primate derived sequences
    # and rename them as simple enhancers
    pr_simple = set(prComEnhID) - set(pr_complex)

    # reassign core and core remodeling columns
    enh.loc[enh.enh_id.isin(pr_simple), "core"] = 1
    enh.loc[enh.enh_id.isin(pr_simple), "core_remodeling"] = 0
    return enh

def format_syndf(enh_age_file, build):

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
    syn_gen_bkgd_file = f"/dors/capra_lab/projects/enhancer_ages/{build}_syn_gen_bkgd.tsv"
    syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t') # read the file
    syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

    syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2"]] # whittle down the df
    syn["mrca"] = syn["mrca"].round(3) # round the ages

    syn = pd.merge(syn, syn_gen_bkgd, how = "left", on = "mrca")
    syn = reEval_PrimComplex(syn)
    syn[syn_cols].to_csv(enh_age_file, sep = "\t", header = False, index = False)
    labeled_syn = add_arch_labels(syn) # add architecture labels

    return labeled_syn


def clean_shuffles(df):

    remove_list = []
    shuf_remove_ids = df.loc[(df["pct_enh"] ==1) & (df.core ==0)]

    df = df.loc[~df.enh_id.isin(shuf_remove_ids)]

    return df


def get_mwu(dataset_col, id_col, test_dif, df, outf):

    results_dict = {}

    for dataset in df[dataset_col].unique(): # simple complex core derived
        print(dataset)
        val_list = [] # two lists for mwu
        median_list = [] # median of list,
        mean_list = [] # mean of list

        id_col_list = df[id_col].unique() # get the variables to compare lists
        comp_name = f"{dataset_col}_{id_col_list[0]}_v_{id_col_list[1]}" # for annotation


        test = df.loc[df[dataset_col] == dataset] # filter dataset on variable

        for comp in id_col_list:
            print(comp)
            values = list(test.loc[(test[id_col] == comp), test_dif])
            val_list.append(values)
            median_list.append(np.median(values))
            mean_list.append(np.mean(values))

        print("val list len", len(val_list))
        stat_, p_ = stats.mannwhitneyu(val_list[0], val_list[1])

        mwu_df = pd.DataFrame({
        "comparison": [comp_name],
        "dataset": [dataset],
        "stat": [stat_],
        "P": [p_],
        f"median{id_col_list[0]}": [median_list[0]],
        f"median{id_col_list[1]}": [median_list[1]],
        f"mean{id_col_list[0]}": [mean_list[0]],
        f"mean{id_col_list[1]}": [mean_list[1]],
        })
        results_dict[dataset] = mwu_df

    results = pd.concat(results_dict.values()) # concat results
    results.to_csv(outf, sep = '\t') # write the results

    return results


def plot_mrca_stratified(x, y, data, hue, palette, outf, build):

    order_dict = {
    "arch" : ["simple", "complex_core", "complex_derived"],
    "mrca_2": list(data.mrca_2.sort_values().unique())
    #"mrca_2" : [0.0, 0.126, 0.131,  0.152, 0.175, 0.308, 0.38, 0.49, 0.656, 0.957]
    }

    hue_order_dict = {
    "arch__": ["simple-cCRE", 'simple-SHUFFLE', 'complex-cCRE', 'complex-SHUFFLE'],
    "id":["cCRE", "SHUFFLE"],
    "id2": ["simple-cCRE", "simple-SHUFFLE", "complex_core-cCRE",
    "complex_core-SHUFFLE", "complex_derived-cCRE", "complex_derived-SHUFFLE"]
    }

    ylab_dict = {
    "syn_len": "syntenic length",
    "pct_enh": "percent\nenhancer sequence"
    }

    xlab_dict = {
    "mrca_2": ["homo", "prim", "euar", "bore", "euth", "ther", "mam", "amni", "tetr", "sarg", "vert"],
    "arch": ["simple", "complex\ncore", "complex\nderived"]
    }

    figsize_dict = {
    "mrca_2": (12,6),
    "arch":(6,6)
    }

    xlabs = xlab_dict[x]
    ylab =  ylab_dict[y]
    order = order_dict[x]
    hue_order = hue_order_dict[hue]
    figsize = figsize_dict[x]



    fig, ax = plt.subplots(figsize = figsize)
    sns.barplot(x=x, y=y, data = data,
                hue = hue,
                palette = palette,
                order = order,
                hue_order = hue_order
                )

    ax.set(ylabel = ylab, xlabel = "")
    ax.set_xticklabels(xlabs, rotation = 90)
    ax.legend(bbox_to_anchor = (1,1))

    plt.savefig(outf, bbox_inches = 'tight')


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
    enhdf = ages.loc[ages["id"] == "cCRE"]

    shuf_ = ages.loc[ages["id"] != "cCRE"]

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


def plot_arch_freq(age_arch_freq, age_freq, build):
    plots = {"age_arch" : age_arch_freq, "age": age_freq}

    for name, frame in plots.items():

        if name == "age_arch": # arrange order and colors of plot.
            frame["plot_hue"] = frame["arch"].astype(str) + "-" + frame["id"].astype(str)
            order = ["simple-cCRE", "simple-SHUFFLE",
            "complex_core-cCRE", "complex_core-SHUFFLE",
            "complex_derived-cCRE", "complex_derived-SHUFFLE"]
            hue = "plot_hue"
            p = archpal

        else:
            order = ["cCRE", "SHUFFLE"]
            hue = "id"
            p = EPAL


        if build == "hg38":
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
        palette = p)

        ax.set_xticklabels(xlabs, rotation = 90)
        ax.legend(bbox_to_anchor = (1,1))

        outf = f"{RE}{name}_freq_per_age.pdf"

        plt.savefig(outf, bbox_inches= "tight")


def plot_arch_fc(age_arch_fc, age_fc, build, arch):
    print(build, arch)
    plots = {"age_arch":age_arch_fc, "age_tfbs": age_fc}
    color_dict ={"simple": amber, "complex_core": purple, "complex_derived": blue}
    archs = ["simple", "complex_core", "complex_derived"]

    # don't plot the entire dataset if you just want to plot FC of one architecture.
    if arch in archs:
        plots.pop("age_tfbs")
        print(plots.keys())

    for name, fc in plots.items():

        fc['log2'] = np.log2(fc["fold_change"])


        if name == "age_arch" and arch in archs:
            order = [arch] # set order
            hue = "arch" # set hue
            data = fc.loc[fc.arch == arch] # filter dataset
            p = color_dict[arch] # get the color corresponding to architecture.

        elif name == "age_arch" and arch == "all":
            order = ["simple", "complex_core", "complex_derived"]
            hue = "arch"
            data = fc
            p = PAL

        elif name =="age_tfbs":
            arch = None
            order = ["FANTOM"]
            hue = "id_y"
            data = fc
            p = EPAL

        backbone_df = age_arch_fc["mrca_2"].copy().drop_duplicates().reset_index()
        data = pd.merge(backbone_df, data, how = "left").fillna(0).sort_values(by = "mrca_2")
        data["counts_y"] = data["counts_y"].astype(int)

        if build == "hg38":
            xlabs = ["Homo", "Prim", "Euar", "Bore", "Euth", "Ther", "Mam", "Amni", "Tetr", "Sarg", "Vert"]

        else:
            xlabs = ["Homo", "Prim", "Euar", "Bore", "Euth", "Ther", "Mam", "Amni", "Tetr", "Vert"]

        sns.set("talk")
        fig, ax = plt.subplots(figsize = (6,6))
        x, y = "mrca_2", "log2"


        splot = sns.barplot(x = x, y=y,
        data = data,
        hue = hue,
        hue_order = order,
        palette = p)


        for n, p in enumerate(splot.patches):
            value = data.iloc[n]["counts_y"]
            splot.annotate(value,
                           (p.get_x() + p.get_width() / 2.,0.05),
                           ha = 'center', va = 'baseline',
                           size=15,
                           rotation = 90,
                           color = "k",
                           xytext = (0, 1),
                           textcoords = 'offset points'
                           )

        ax.set(ylabel = "Fold-Change v. Bkgd\n(log2-scaled)")
        ax.set_xticklabels(xlabs, rotation = 90)
        ax.legend(bbox_to_anchor = (1,1))
        outf = f"{RE}{name}_{arch}_fold_change_per_age.pdf"


    plt.savefig(outf, bbox_inches= "tight")


#%%


enh = format_syndf(cCRE, GENOME_BUILD)
enh["id"] = "cCRE"

shuf = format_syndf(SHUFF, GENOME_BUILD)
shuf["id"] = "SHUFFLE"


df = pd.concat([enh, shuf])
df.shape
df.drop_duplicates().shape
enh.groupby(["core_remodeling", "core", "mrca_2"])["enh_id"].count()
shuf.groupby(["core_remodeling", "core"])["enh_id"].count()

enh[["core_remodeling","enh_id"]].drop_duplicates().groupby("core_remodeling").count()
#%% summarize architecture lengths per enhancer


sum_archlen = df.groupby(['id', 'enh_id', 'enh_len', 'core_remodeling', 'core'])["syn_len"].sum().reset_index().drop_duplicates()

sum_archlen = add_arch_labels(sum_archlen) # add architecture labels

 # calculate percent arch per enhancer
sum_archlen["pct_enh"] = sum_archlen.syn_len.divide(sum_archlen.enh_len)


# add oldest ages of enhancer information
enh_mrcas = df.groupby("enh_id")[["mrca_2", "taxon2"]].max().reset_index()
sum_archlen = pd.merge(sum_archlen, enh_mrcas, how = "left", on = "enh_id")


shuf_remove_ids = sum_archlen.loc[(sum_archlen["mrca_2"] == 0) & (sum_archlen.core ==0), "enh_id"]
sum_archlen = sum_archlen[~sum_archlen["enh_id"].isin(shuf_remove_ids)]

#%% summarize simple v. complex architectures only

sum_archlen["arch_"] = "simple"
sum_archlen.loc[sum_archlen.core_remodeling ==1, "arch_"] = "complex"
sum_archlen["arch__"] = sum_archlen["arch_"] + "-" + sum_archlen['id']

# add column for plotting purposes
sum_archlen["id2"] = sum_archlen["arch"] + "-" + sum_archlen["id"]


#%% plot syn lengths per age, arch.

x = "mrca_2"
y = "syn_len"
data = sum_archlen#.sample(frac = 0.25)
hue = "arch__"
palette = enhpal
outf = f"{RE}ALL_mrca_x_syn_lengths_arch.pdf"

plot_mrca_stratified(x, y, data, hue, palette, outf, GENOME_BUILD)


#%% plot syn length of segments w/o mrca


x = "arch"
y = "syn_len"
hue = "id"
palette = EPAL
outf = f"{RE}ALL_sum_arch_length_all_cCRE_shuffle.pdf"

plot_mrca_stratified(x, y, data, hue, palette, outf, GENOME_BUILD)


#%% stratify syn length by age

x = "mrca_2"
y = "pct_enh"
hue = "id2"
palette = archpal
outf = f"{RE}ALL_cCRE_percent_arch_mrca_2.pdf"
plot_mrca_stratified(x, y, data, hue, palette, outf, GENOME_BUILD)


#%% plot percent architecture w/o age


x = "arch"
y = "pct_enh"
hue = "id"
outf = f"{RE}ALL_cCRE_percent_all.pdf"

plot_mrca_stratified(x, y, data, hue, palette, outf, GENOME_BUILD)


#%% stratify percent enh by age

x = "mrca_2"
y = "syn_len"
hue = "id2"
palette = archpal
outf = f"{RE}ALL_cCRE_syn_len_mrca_2.pdf"
plot_mrca_stratified(x, y, data, hue, palette, outf, GENOME_BUILD)


#%% mean lengths

comp = "all_cCRE_v_shuffle"

lens = sum_archlen.groupby(["id", "arch"])["syn_len"].mean().reset_index()
pcts = sum_archlen.groupby(["id", "arch"])["pct_enh"].mean().reset_index()

outf = f"{RE}{comp}_lens_metrics.tsv"
lens.to_csv(outf, sep = '\t')

outf = f"{RE}{comp}_pcts_metrics.tsv"
pcts.to_csv(outf, sep = '\t')

#%% do mwu on architecture, comparing cCRE v. SHUFFLE ids

dataset_col = "arch" # first level - simple, complex core, complex derived
id_col = "id" # do mwu separating on this variable
test_dif = "syn_len" # calculate mwu from these variables' values


outf = f"{RE}{comp}_MWU_arch_lens.tsv"
mwu_df = get_mwu( dataset_col, id_col, test_dif, sum_archlen, outf)

#%% do mwu on architecture, comparing cCRE v. cCRE ids

dataset_col = "arch" # first level - simple, complex core, complex derived
id_col = "arch" # do mwu separating on this variable
test_dif = "syn_len" # calculate mwu from these variables' values


outf = f"{RE}{comp}_MWU_arch_v_arch_lens.tsv"
mwu_df = get_mwu( dataset_col, id_col, test_dif, sum_archlen, outf)


#%% get arch counts

arch_count = enh[["core_remodeling", "enh_id"]].drop_duplicates().groupby("core_remodeling")["enh_id"].count().reset_index()
arch_count.columns = ["core_remodeling", "count_arch"]
arch_count["total"] = arch_count.count_arch.sum()
arch_count["percent"] = arch_count.count_arch.divide(arch_count.total)

outf = f"{RE}{comp}_arch_count.tsv"
arch_count.to_csv(outf, sep = '\t')


#%% GET MRCA FREQUENCIES

cols = ["id", "arch", "mrca_2"]
var = "mrca_2"
age_arch_freq, age_arch_fc = MRCA_frequency(df, cols, var)


cols = ["id", "mrca_2"]
var = "mrca_2"
age_freq, age_fc = MRCA_frequency(df, cols, var)


#%%

plot_arch_freq(age_arch_freq, age_freq, GENOME_BUILD)
plot_arch_fc(age_arch_fc, age_fc, GENOME_BUILD, "all")
plot_arch_fc(age_arch_fc, age_fc, GENOME_BUILD, "complex_derived")
plot_arch_fc(age_arch_fc, age_fc, GENOME_BUILD,  "complex_core")
plot_arch_fc(age_arch_fc, age_fc, GENOME_BUILD, "simple")

#age_arch_fc.loc[age_arch_fc.arch == "complex_core"]
#%%
df.head()
#%%
der = df.loc[df.core == 0]
der_enh = der.loc[der.id == "cCRE", "mrca_2"]
der_shuf = der.loc[der.id == "SHUFFLE", "mrca_2"]
stats.mannwhitneyu(der_enh, der_shuf)
der.groupby("id")["mrca_2"].median()
"""
id
cCRE     0.175
SHUFFLE    0.152
Name: mrca_2, dtype: float64
"""
der.groupby("id")["mrca_2"].mean()
"""
id
cCRE     0.206010
SHUFFLE    0.172115
Name: mrca_2, dtype: float64
"""
