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
FANTOMFILE = "All_FANTOM_x_ENCODE_hg19.bed"
FANTOM = os.path.join(FANTOMPATH, FANTOMFILE)

SHUFPATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/shuffle/first_round_breaks"
SHUFFILE = "no-exon_syn_breaks_shuf-all_fantom_enh_ages.bed"
SHUFF = os.path.join(SHUFPATH, SHUFFILE)


RE = "/dors/capra_lab/projects/enhancer_ages/landscape/results/fantom/tfbs_overlap_analysis/"

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

#%%


def add_arch_labels(df):

    df["arch"] = "complex_core"

    df.loc[df.core_remodeling ==0, "arch"] = "simple"
    df.loc[df.core ==0, "arch"] = "complex_derived"

    return df


def remove_tf_cols(df):
    df_labels = list(df) # get the columns

    tf_df_labels = ['chr_t', 'start_t', 'end_t','tf_id',
    'peak_len', 'tf','cl','tf_overlap'] # get the columns related to the encode TFBS info

    # drop the TF info.
    # We want to summarize at the syntenic block and enhancer level.
    # Not individual TFs.

    df = df.drop(tf_df_labels, axis = 1)

    # get the new columns
    new_labels = list(set(df_labels)- set(tf_df_labels))

    # subtract the column that we want to groupby
    new_labels_no_tfbin = list(set(new_labels)- set(["tf_bin"]))

    # get the max TF_bin value to demarcate whether syn/enh do or do not overlap TFBS
    df = df.groupby(new_labels_no_tfbin)["tf_bin"].max().reset_index()

    return df


def format_syndf(enh_age_file):

    syn = pd.read_csv(enh_age_file, sep ='\t')

    syn["syn_id"] = syn["#chr_syn"] + ":" + syn.start_syn.map(str) + "-" + syn.end_syn.map(str)

    syn["syn_len"] = syn.end_syn - syn.start_syn
    syn["enh_len"] = syn.end - syn.start
    syn = syn.loc[syn.syn_len >5]
    syn = syn.loc[syn.enh_len >=100]

    syn["tf"] = syn.tf_id.apply(lambda x: x.split("_")[0])

    syn["tf_bin"] = 0
    syn.loc[syn["tf_overlap"]>5, "tf_bin"] = 1 # if more than 5bp overlap, count the TF overlap


    # age and taxon file
    syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_taxon.bed"
    syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t') # read the file
    syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

    syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2"]] # whittle down the df
    syn["mrca"] = syn["mrca"].round(3) # round the ages

    syn = pd.merge(syn, syn_gen_bkgd, how = "left", on = "mrca")

    labeled_syn = add_arch_labels(syn) # add architecture labels

    labeled_syn = remove_tf_cols(labeled_syn)

    return labeled_syn


def clean_shuffles(df):

    remove_list = []
    shuf_remove_ids = df.loc[(df["pct_enh"] ==1) & (df.core ==0)]

    df = df.loc[~df.enh_id.isin(shuf_remove_ids)]

    return df

def summarize_arch(df):

    sum_archlen = df.groupby(['id', 'enh_id', 'enh_len', 'core_remodeling', 'core'])["syn_len"].sum().reset_index().drop_duplicates()

    sum_archlen = add_arch_labels(sum_archlen) # add architecture labels

     # calculate percent arch per enhancer
    sum_archlen["pct_enh"] = sum_archlen.syn_len.divide(sum_archlen.enh_len)


    # add oldest ages of enhancer information
    enh_mrcas = df.groupby("enh_id")[["mrca_2", "taxon2"]].max().reset_index()
    sum_archlen = pd.merge(sum_archlen, enh_mrcas, how = "left", on = "enh_id")


    # rename the archiectures for plotting.

    sum_archlen["arch_"] = "simple"
    sum_archlen.loc[sum_archlen.core_remodeling ==1, "arch_"] = "complex"
    sum_archlen["arch__"] = sum_archlen["arch_"] + "-" + sum_archlen['id']
    sum_archlen["arch__"].unique()

    return sum_archlen

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
    enhdf = ages.loc[ages["id"] == "YES_TFBS"]

    shuf_ = ages.loc[ages["id"] != "YES_TFBS"]

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

def plot_arch_freq(age_arch_freq, age_freq, comp):
    plots = {"YN_TFBS_age_arch" : age_arch_freq, "YN_TFBS_age": age_freq}

    for name, frame in plots.items():

        if name == "YN_TFBS_age_arch": # arrange order and colors of plot.
            frame["plot_hue"] = frame["arch"].astype(str) + "-" + frame["id"].astype(str)

            hue = "plot_hue"
            order = ["simple-YES_TFBS", "simple-NO_TFBS",
            "complex_core-YES_TFBS","complex_core-NO_TFBS",
            "complex_derived-YES_TFBS", "complex_derived-NO_TFBS"]

        else:
            order = ["YES_TFBS", "NO_TFBS"]
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

        outf = f"{RE}{name}_{comp}_age_freq.pdf"

        plt.savefig(outf, bbox_inches= "tight")

def plot_arch_fc(age_arch_fc, age_fc, comp):

    plots = {"YN_TFBS_age_arch":age_arch_fc,} #"YN_TFBS_age": age_fc}

    for name, fc in plots.items():

        fc['log2'] = np.log2(fc["fold_change"])
        print(list(fc))
        if name == "YN_TFBS_age_arch":
            order = ["simple", "complex_core", "complex_derived"]
            hue = "arch"

        else:
            order = ["YES_TFBS"]
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
        palette = PAL)
        ax.set(ylabel = "Fold-Change Y_TFBS_overlap v. N_TFBS_overlap\n(log2-scaled)")
        ax.set_xticklabels(xlabs, rotation = 90)
        ax.legend(bbox_to_anchor = (1,1))

        outf = f"{RE}{name}_{comp}fc.pdf"

    plt.savefig(outf, bbox_inches= "tight")


#%%


enh = format_syndf(FANTOM)
enh["id"] = "YES_TFBS"

#%% find syn_id and enh_id for regions that do not overlap TFBS in ENCODE
no_syn_overlap = enh.loc[enh.tf_bin == 0, ["enh_id", "core_remodeling", "syn_id"]] # find all syntenic blocks that do not overlap TFBS

no_syn_overlap.count() #14808 syntenic blocks do not overlap TFBS
simple = no_syn_overlap[["core_remodeling", "enh_id"]].drop_duplicates()
simple.groupby(["core_remodeling"])["enh_id"].count()


no_syn_overlap.enh_id.unique().size # 10591 enhancers have syntenic blocks that do not overlap TFBS
no_syn_overlap_enh_id = list(no_syn_overlap.enh_id.unique())

# find only enh_ids - entire enhancer does not overlap any TFBS
enh_tf = enh.groupby("enh_id")["tf_bin"].sum().reset_index()
no_enh_overlap = enh_tf.loc[enh_tf.tf_bin == 0, "enh_id"] # find all syntenic blocks that do not overlap TFBS
no_enh_overlap.size #6375 enhancer do not have any overlap w/ TFBS

comparison_dict = {
"syn_tfbs_only": no_syn_overlap_enh_id,
"enh_tfbs_only": no_enh_overlap
}

simple = enh.loc[enh.enh_id.isin(no_enh_overlap)].drop_duplicates()
simple = simple[["core_remodeling", "enh_id"]].drop_duplicates()
simple.groupby(["core_remodeling"])["enh_id"].count()


#%%

for comp, overlap_list in comparison_dict.items():

    # "shuf" is an old convention and I am reusing old code.
    shuf = enh.loc[enh.enh_id.isin(overlap_list)]
    # Shuf variables WERE mapped to shuffle datasets,
    # but here ARE mapped to enhancer datasets that do not overlap enhancers.
    shuf["id"] = "NO_TFBS"

    tf_overlap = enh.loc[~enh.enh_id.isin(overlap_list)]

    outf = f"{FANTOMPATH}{comp}.txt"
    tf_overlap["enh_id"].to_csv(outf, sep  = '\t', index = False)
    df = pd.concat([tf_overlap, shuf])


    # summarize architecture lengths per enhancer
    # get rid of the TFBS annotations. They create redundant rows for enhancer arch information in the df
    sum_archlen = summarize_arch(df)

    # plot syn lengths per age, arch.


    x = "mrca_2"
    y = "syn_len"
    data = sum_archlen.sample(frac = 0.25)
    hue = "arch__"


    xlabs = ["homo", "prim", "euar", "bore", "euth", "ther", "mam", "amni", "tetr", "vert"]
    order = ["simple-YES_TFBS", 'simple-NO_TFBS', 'complex-YES_TFBS', 'complex-NO_TFBS']

    sns.set("poster")
    fig, ax = plt.subplots(figsize = (12,6))
    sns.barplot(x=x, y=y, data = data,
                hue = hue,
                palette = enhpal, hue_order = order )

    ax.set(ylabel = "syntenic length", xlabel = "")
    ax.set_xticklabels(xlabs, rotation = 90)
    ax.legend(bbox_to_anchor = (1,1))
    outf = f"{RE}YN_TFBS_{comp}_mrca_x_syn_lengths_arch.pdf"

    plt.savefig(outf, bbox_inches = 'tight')


    # plot length of segments


    ticklabs = ["simple", "complex\ncore", "complex\nderived"]

    x = "arch"
    y = "syn_len"
    hue = "id"
    data = sum_archlen
    order = ["simple", "complex_core", "complex_derived"]
    hue_order = ["YES_TFBS", "NO_TFBS"]

    fig, ax = plt.subplots(figsize = (6,6))

    sns.barplot(x=x, y=y, data = data, hue = hue,
    palette = es, order = order, hue_order=hue_order)
    ax.set_xticklabels(ticklabs)
    ax.legend(bbox_to_anchor = (1,1))
    ax.set(ylabel = "sum length (bp)", xlabel = "")
    outf = f"{RE}YN_TFBS_{comp}_sum_arch_length_all_fantom_shuffle.pdf"
    plt.savefig(outf, bbox_inches = "tight")


    # calculate some stats and quantify some measures on core v. derived lengths, percents


    lens = sum_archlen.groupby(["id", "arch"])["syn_len"].mean().reset_index()
    pcts = sum_archlen.groupby(["id", "arch"])["pct_enh"].mean().reset_index()

    outf = f"{RE}YN_TFBS_{comp}_lens_metrics.tsv"
    lens.to_csv(outf, sep = '\t')

    outf = f"{RE}YN_TFBS_{comp}_pcts_metrics.tsv"
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
    outf = f"{RE}YN_TFBS_{comp}_MWU_core_der_lens.tsv"

    mwu_df.to_csv(outf, sep = '\t')


    #

    x = "arch"
    y = "pct_enh"
    hue = "id"


    sns.set("poster")

    fig, ax = plt.subplots(figsize=(6,6))

    # plot
    sns.barplot(x = x, y = y, data = data, order = order, hue_order = hue_order,
    hue = hue, palette =archpal, errwidth= 10)

    ax.set(ylabel = "Percent of enhancer length", xlabel = "")
    ax.set_xticklabels(ticklabs)
    ax.legend(bbox_to_anchor = (1,1))
    outf = f"{RE}YN_TFBS_{comp}_fantom_percent_all.pdf"
    plt.savefig(outf, bbox_inches = "tight")

    #
    sum_archlen["id2"] = sum_archlen["arch"] + "-" + sum_archlen["id"]


    x = "mrca_2"
    y = "pct_enh"
    hue = "id2"
    order2 = ["simple-YES_TFBS", "simple-NO_TFBS",
    "complex_core-YES_TFBS","complex_core-NO_TFBS",
     "complex_derived-YES_TFBS", "complex_derived-NO_TFBS"]

    sns.set("poster")


    fig, ax = plt.subplots(figsize=(12,6))
    sns.barplot(x = x, y = y, data = data, hue = hue, hue_order = order2,
    palette = archpal, errwidth= 5)

    ax.set(ylabel = "Percent of enhancer length", xlabel = "")

    ax.set_xticklabels(xlabs, rotation = 90)
    ax.legend(bbox_to_anchor = (1,1))

    outf = f"{RE}YN_TFBS_{comp}_percent_arch_mrca_2.pdf"
    plt.savefig(outf, bbox_inches = "tight")


    #

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

    outf = f"{RE}YN_TFBS_{comp}_syn_len_arch_mrca_2.pdf"
    plt.savefig(outf, bbox_inches = "tight")


    cols = ["id", "arch", "mrca_2"]
    var = "mrca_2"
    age_arch_freq, age_arch_fc = MRCA_frequency(df, cols, var)


    cols = ["id", "mrca_2"]
    var = "mrca_2"
    age_freq, age_fc = MRCA_frequency(df, cols, var)


    GENOME_BUILD = "hg19"
    plot_arch_freq(age_arch_freq, age_freq, comp)
    plot_arch_fc(age_arch_fc, age_fc, comp)

    age_arch_fc
