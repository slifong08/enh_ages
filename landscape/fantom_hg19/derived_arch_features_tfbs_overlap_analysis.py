import matplotlib.pyplot as plt
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


#%%


enh = format_syndf(FANTOM)
enh["id"] = "YES_TFBS"
#%%

enh.head()
#%% find syn_id and enh_id for regions that do not overlap TFBS in ENCODE
no_syn_overlap = enh.loc[enh.tf_bin == 0, ["enh_id", "syn_id"]]

no_syn_overlap.count() #14808 syntenic blocks do not overlap TFBS

no_syn_overlap.enh_id.unique().size # 10591 enhancers have syntenic blocks that do not overlap TFBS

# find only enh_ids - entire enhancer does not overlap any TFBS
enh_tf = enh.groupby("enh_id")["tf_bin"].sum().reset_index()

no_enh_overlap = enh_tf.loc[enh_tf.tf_bin == 0, "enh_id"]
no_enh_overlap.size #6375 enhancer do not have any overlap w/ TFBS
no_enh_overlap
comparison_dict = {
"syn_tfbs_only": list(no_syn_overlap.enh_id.unique()),
"enh_tfbs_only": no_enh_overlap
}

#%%
for comp, overlap_list in comparison_dict.items():

    # "shuf" is an old convention and I am reusing old code.
    shuf = enh.loc[enh.enh_id.isin(overlap_list)]
    # Shuf variables WERE mapped to shuffle datasets,
    # but here ARE mapped to enhancer datasets that do not overlap enhancers.
    shuf["id"] = "NO_TFBS"

    tf_overlap = enh.loc[ ~enh.enh_id.isin(overlap_list)]
    tf_overlap.enh_id.unique().size #22367 enhancers overlap TFBS
    outf = f"{FANTOMPATH}{comp}.txt"
    tf_overlap["enh_id"].to_csv(outf, sep  = '\t', index = False)
    df = pd.concat([tf_overlap, shuf])


    # summarize architecture lengths per enhancer
    # get rid of the TFBS annotations. They create redundant rows for enhancer arch information in the df


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
    outf = f"{RE}TFBS_OVERLAP_{comp}_mrca_x_syn_lengths_arch.pdf"

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
    outf = f"{RE}TFBS_OVERLAP_{comp}_sum_arch_length_all_fantom_shuffle.pdf"
    plt.savefig(outf, bbox_inches = "tight")


    # calculate some stats and quantify some measures on core v. derived lengths, percents


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
    outf = f"{RE}TFBS_OVERLAP_{comp}_fantom_percent_all.pdf"
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

    outf = f"{RE}TFBS_OVERLAP_{comp}_percent_arch_mrca_2.pdf"
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

    outf = f"{RE}TFBS_OVERLAP_{comp}_syn_len_arch_mrca_2.pdf"
    plt.savefig(outf, bbox_inches = "tight")
