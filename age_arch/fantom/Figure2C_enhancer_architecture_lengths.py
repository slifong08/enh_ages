import glob
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.ticker as ticker
import numpy as np
import os, sys
from scipy import stats
import seaborn as sns
import statsmodels
import statsmodels.api as sm

RE ="/dors/capra_lab/projects/enhancer_ages/fantom/results/for_publication/age_breaks/"

es_colors = [ "slate grey","greyish"]
es_pal = sns.xkcd_palette(es_colors)
sns.palplot(es_pal)


PATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/"

enh = "%snon-genic/no-exon_all_fantom_enh/breaks/no-exon_all_fantom_enh_ages_enh_age_arch_summary_matrix.bed" % PATH

shuf = "%snon-genic/no-exon_shuffle_fantom_age_arch_summary_matrix.bed" % PATH


#%% functions

def get_syn_gen():
    syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
    syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t') # read the file
    syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

    syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2", "mya", "mya2"]] # whittle down the df

    return syn_gen_bkgd


def load_df(f):

    print(f)

    cols = ["#chr_enh", "start_enh", "end_enh", "enh_id", "id",
        "seg_index", "core_remodeling", "arch", "mrca",
        "taxon", "mrca_2", "taxon2", "mya", "mya2"]

    df = pd.read_csv(f, sep = '\t')
    print(df.head())
    #df = df.rename(columns={'datatype': 'id', 'chr_enh': '#chr_enh'})
    df["enh_len"] = df.end_enh - df.start_enh

    return df


def generate_shuf_id(shuffle, final_merge):

    enh_shape = final_merge.shape[0] # get length of enh dataframe

    no_shufs = round(shuffle.shape[0]/enh_shape, 0) # find number of shuffled datasets to label

    sampleD = {} # collect sample results here.

    for n in np.arange(no_shufs):
        sample = shuffle.sample(enh_shape, replace = False) # sample the shuffled dataframe space
        sample["shuf_id"] = n # add a pseudo shuffle id
        sampleD[n] = sample # add sample to the collection dictionary

    new_shuf = pd.concat(sampleD.values())
    return new_shuf


def get_lens(df, syn_gen_bkgd):

    max_val = "mrca_2"

    if "shuf_id" in list(df):
        cols = ["enh_id", "core_remodeling", "enh_len", "shuf_id", "arch"]
        dataset = "SHUFFLE"
    else:

        cols = ["enh_id", "core_remodeling","enh_len", "arch"]
        dataset = "FANTOM"

    lens = df.groupby(cols)[max_val].max().reset_index()
    print(lens[max_val].unique())
    lens.mrca_2 = lens.mrca_2.round(3)

    lens = pd.merge(lens, syn_gen_bkgd[[ "mrca_2","taxon2", "mya2"]], how = "left", on = "mrca_2")
    lens = lens.drop_duplicates()

    lens["datatype"]= dataset
    print(lens.shape)

    return lens


def get_len_freq(df):

    x = df.shape[0] # how many enhancers are there?

    len_freq = df.groupby(["mrca_2", "core_remodeling"])["enh_id"]\
    .count().divide(x).round(3).reset_index()

    if "shuf_id" in list(df):
        freq = "shuf_freq"
    else:
        freq = "enh_freq"

    len_freq.columns = ["mrca_2", "core_remodeling", freq]

    return len_freq


def plot_fig2c_arch_len(enh_lens, shuf_len):

    xlabs = [
    "Homo",
    "Prim", "Euar", "Bore", "Euth", "Ther"
    , "Mam", "Amni", "Tetr", "Vert"]

    e_colors = [ "amber", "faded green"]
    e_pal = sns.xkcd_palette(e_colors)
    s_colors = [ "greyish", "slate grey"]
    s_pal = sns.xkcd_palette(s_colors)

    hue_order = ["FANTOM", "Shuffle"]
    fig,(ax1) = plt.subplots(figsize = (8, 8))
    order =["Simple", "Complexenh"]

    x = "taxon2"
    y = "enh_len"
    hue = "arch"
    data = enh_lens.sort_values(by = "mrca_2")
    sns.set("poster")
    sns.barplot(x = x, y = y, data = data,
     ax = ax1, hue = hue,  palette = e_pal,
     estimator = np.median)#showfliers=False)

    simple_lens = enh_lens.loc[enh_lens.arch.str.contains("imple"), "enh_len"]
    shuf_simple_lens = shuf_len.loc[shuf_len.arch.str.contains("imple"), "enh_len"]

    # MWU on simple v. shuffle lengths
    ms, msp = stats.mannwhitneyu(simple_lens, shuf_simple_lens)
    print("simple", ms, msp)

    # MWU on complex v. shuffle lengths
    complex_lens = enh_lens.loc[enh_lens.arch.str.contains("omplex"), "enh_len"]
    shuf_complex_lens = shuf_len.loc[shuf_len.arch.str.contains("omplex"), "enh_len"]

    mc, mcp = stats.mannwhitneyu(complex_lens,shuf_complex_lens)

    print("complex", mc, mcp)

    ax1.set(ylabel= "Enhancer Length (bp)", ylim = (190,400), xlabel = "")
    ax1.set_xticklabels(xlabs, rotation = 90, horizontalalignment = "left")


    ax1.get_legend().remove()
    plt.savefig("%sfig2c-Fantom_ENH_MRCA_x_LEN_ENH.pdf" % RE, bbox_inches = "tight")


def plot_Figs9_lengths_per_age(enh_lens, shuf_len):

    xlabs = [
        "Homo",
        "Prim", "Euar", "Bore", "Euth", "Ther"
        , "Mam", "Amni", "Tetr", "Vert"]
    len_concat = pd.concat([enh_lens, shuf_len])

    len_concat["datatype2"] = len_concat.datatype + "-"+ len_concat.arch
    len_concat["datatype2"].unique()


    e_colors = [ "amber","greyish", "faded green", "slate grey"]
    e_pal = sns.xkcd_palette(e_colors)
    hue_order = ['FANTOM-simple','SHUFFLE-simple','FANTOM-complexenh','SHUFFLE-complexenh']


    fig,(ax1) = plt.subplots(figsize = (8, 8))
    order =["Simple", "Complexenh"]
    sns.barplot(y = "enh_len", x = "taxon2",
                data = len_concat.sort_values(by = "mrca_2"), ax = ax1,
                hue = "datatype2",  palette = e_pal,
                hue_order = hue_order,
                estimator = np.median)#showfliers=False)


    ax1.set(ylabel = "Enhancer Length (bp)", xlabel = "")
    ax1.set_xticklabels(xlabs, rotation = 90, horizontalalignment = "left")

    ax1.legend().remove()
    plt.savefig("%sfigS9A-Fantom_ENH_MRCA_x_LEN_ENH.pdf" % RE, bbox_inches = "tight")

    # plot LM
    e_colors = [ "amber","greyish", "faded green", "slate grey"]
    e_pal = sns.xkcd_palette(e_colors)
    hue_order = ['FANTOM-simple','SHUFFLE-simple','FANTOM-complexenh','SHUFFLE-complexenh']

    order =["Simple", "Complexenh"]
    sns.lmplot(y = "enh_len", x = "mya2", data = len_concat,\
                hue = "datatype2",  palette = e_pal,
                hue_order = hue_order, scatter = False)
    plt.savefig("%sfigS9B-LM_enh_len.pdf" % RE, bbox_inches = "tight")


#%% load dataframe


shuffle_ = load_df(shuf)

final_merge = load_df(enh)


shuffle = generate_shuf_id(shuffle_, final_merge)

syn_gen_bkgd = get_syn_gen()
#%%


shuf_len = get_lens(shuffle, syn_gen_bkgd)
enh_lens = get_lens(final_merge, syn_gen_bkgd)

#%% sample 1/10th of shuffle dataset and combine w/ enhancers for comparison


sample_shuf = shuf_len.sample(frac = 0.25)

lens = pd.concat([enh_lens, sample_shuf])


#%% calculate frequency of obs and expected enhancer lengths from fantom, shuffle

sample_shuf.mrca_2.unique()

enh_len_freq = get_len_freq(enh_lens)
final_merge.mrca.unique()
shuf_len_freq = get_len_freq(sample_shuf)


# calculat difference between exp and observed frequency
len_freq = pd.merge(enh_len_freq, shuf_len_freq, how = 'left' )
len_freq

#%%

enh_lens = pd.merge(enh_lens, len_freq, how = "left")
shuf_len = pd.merge(shuf_len, len_freq, how = "left")
shuf_len.head()

#%% plot enhancer architecture length per age
plot_fig2c_arch_len(enh_lens, shuf_len)

""" RESULTS enhancer lengths v. expected shuffle lengths for simple, complex
simple 15249058638.5 2.627822101951775e-05
complex 6146319551.5 1.3788416832439736e-06
"""
#%% Supplemental Figure 9 - FANTOM/SHUFFLE enhancer age x arch length
plot_Figs9_lengths_per_age(enh_lens, shuf_len)
#%%
outf = f"{RE}mrca2_arch_counts.tsv"
mrca_arch = final_merge.groupby([ "arch", "taxon2", "mrca_2"])["enh_id"].count().reset_index().sort_values(by = ["mrca_2"])
mrca_arch_p = pd.pivot(data = mrca_arch, columns = ["taxon2", "mrca_2"], index = "arch", values = "enh_id")#.sort_values(by = "mrca_2")
mrca_arch_p
mrca_arch_p.to_csv(outf, sep = '\t', index = False)
#%%
final_merge.groupby("arch")["enh_id"].count()
#%%
seg_arch = final_merge.groupby([ "arch", "seg_index"])["enh_id"].count().reset_index()
seg_arch
