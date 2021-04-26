import os, sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression
import seaborn as sns


PATH = "/dors/capra_lab/users/fongsl/tyler/data/CON_ACC/"
RE = "/dors/capra_lab/users/fongsl/tyler/results/CON_ACC/"
dataset_list = ["species_specific_10k", "all"]
multiz_list = [20, 30, 100]

fs = {}
for msa in multiz_list:
    for dataset in dataset_list:

        chr_num = "chr21"

        f = f"{PATH}{dataset}/multiz{msa}way/{chr_num}_con_acc.bed"
        names = ["#chr", "file", "file_type", "start", "end", "log_p", "multiz", "??", "id"]
        df = pd.read_csv(f, sep = '\t', header = None, names = names)

        df["multiz"] = f"{msa}_way"
        df["dataset"] = dataset

        key = f"{msa}way-{dataset}"
        df["dataset_id"] = key
        df["enh_id"] = df["#chr"] + ":" + df["start"].map(str) + "-" + df["end"].map(str)

        fs[key] = df

def separate_shared_specific(df):

    all_ = df.loc[df.dataset == "all", "enh_id"].drop_duplicates() # get all enh_id
    ss_ = df.loc[df.dataset != "all", "enh_id"].drop_duplicates() # ge species_specific enh_id

    shared_ = list(set(all_)-(set(ss_))) # get shared (non-species_specific) enh_ids

    n_shared, n_specific, n_all = len(shared), len(ss_), len(all_)
    n_list = [n_shared, n_specific, n_all]
    print(n_shared, n_specific, n_all)

    # get dataframes
    shared_df = df.loc[df.enh_id.isin(shared_)]
    specific_df = df.loc[df.enh_id.isin(ss_)]

    # get logp results per multiz
    log_p_dict = {}
    ways = list(df.multiz.unique())

    for way in ways:

        shared_p = shared_df.loc[shared_df.multiz == way, "log_p"]
        specific_p = specific_df.loc[specific_df.multiz == way, "log_p"]
        log_p_dict[way] = [shared_p, specific_p]

    return shared_df, specific_df, log_p_dict, n_list


def plot_joint(x, y, data, RE):

    sns.set("talk")
    g = sns.jointplot(data = data, x= x,  y = y, kind = "hist", marginal_ticks =True)

    # run pearson's R
    X, Y = table[x], table[y]
    r, p = stats.pearsonr(X, Y)


    g.ax_joint.annotate(f'pearsons r= {r:.3f}, p = {p:.3f}',
                        xy=(0.05, -0.25), xycoords='axes fraction',
                        ha='left', va='center',
                        bbox={'boxstyle': 'round', 'fc': 'white', 'ec': 'navy'})

    #g.ax_joint.scatter(X, Y)
    g.set_axis_labels(xlabel=x, ylabel=y, size=15)

    # save the plot
    outf = f"{RE}{x}_v_{y}.pdf"
    plt.savefig(outf, bbox_inches= "tight")


def plot_dist_p(log_p_dict, way):

    shared_p, specific_p = log_p_dict[way]
    result_stat, p = stats.mannwhitneyu(shared_p, specific_p)
    shared_med, spec_med= round(shared_p.median(), 3), round(specific_p.median(), 3)
    xlabel = f"phyloP logp\n(-) = less, (+) = more conservation\nmwu p = {p}\nmedian logp shared ={shared_med}, specific = {spec_med}"


    fig, ax = plt.subplots()
    sns.distplot(shared_p, hist = False, label = f"shared n = {n_list[0]}")
    sns.distplot(specific_p, hist = False, label = f"specific n = {n_list[1]}")
    ax.set(title = f"chr21 - {way}",
    xlabel = xlabel,
    ylabel = "KDE Density")
    ax.legend(bbox_to_anchor =(1,1))

    outf = f"{RE}ch21_shared_v_spec_{way}_kde.pdf"
    plt.savefig(outf, bbox_inches = "tight")


    fig, ax = plt.subplots()
    sns.distplot(shared_p, kde = False, label = f"shared n = {n_list[0]}")
    sns.distplot(specific_p, kde = False, label = f"specific n = {n_list[1]}")
    ax.set(title = f"chr21 - {way}",
    xlabel = xlabel,
    ylabel = "N")
    ax.legend(bbox_to_anchor =(1,1))

    outf = f"{RE}ch21_shared_v_spec_{way}_hist.pdf"
    plt.savefig(outf, bbox_inches = "tight")


#%%


df = pd.concat(fs.values())


ways = list(df.multiz.unique())


table = pd.pivot(df, index = "id", columns = "multiz", values = "log_p")

# split up df
shared_df, specific_df, log_p_dict, n_list = separate_shared_specific(df)


#%%



#%% compare 20 v. 30 way


x = "20_way"
y = "30_way"
data = table
plot_joint(x, y, data, RE)


#%% compare 20 v. 100 way


x = "20_way"
y = "100_way"
data = table
plot_joint(x, y, data, RE)


#%% compare 30 v. 100 way


x = "30_way"
y = "100_way"
plot_joint(x, y, data, RE)

#%% MWU p shared v. specific



for way in ways:

    plot_dist_p(log_p_dict, way)
    print(way)

#%%
