import os, sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression
import seaborn as sns


PATH = "/dors/capra_lab/users/fongsl/tyler/data/CON_ACC/"
RE = "/dors/capra_lab/users/fongsl/tyler/results/CON_ACC/"
dataset_list = ["phastCons", "hars", "all", 'hu_specific', "rhe_specific",
 "all/subtract_te", 'hu_specific/subtract_te', "rhe_specific/subtract_te"]
multiz_list = [20, 30, 100]

fs = {}
for msa in multiz_list:
    for dataset in dataset_list:

        chr_num = "chr21"

        f = f"{PATH}{dataset}/multiz{msa}way/{chr_num}_con_acc.bed"
        if os.path.exists(f) == True:
            names = ["#chr", "file", "file_type", "start", "end", "log_p", "multiz", "??", "id"]
            df = pd.read_csv(f, sep = '\t', header = None, names = names)

            df["multiz"] = f"{msa}_way"
            df["dataset"] = dataset

            key = f"{msa}way-{dataset}"
            df["dataset_id"] = key
            df["enh_id"] = df["#chr"] + ":" + df["start"].map(str) + "-" + df["end"].map(str)

            fs[key] = df

def separate_shared_specific(df):

    id_dict = {}
    n_dict = {}

    for dataset in df.dataset.unique(): # get ids in each dataset
        ids = df.loc[df.dataset == dataset, "enh_id"].drop_duplicates()
        id_dict[dataset] = list(ids) # get all enh_id
        n_dict[dataset] = len(ids) # count the number of ids in the list

    sub_from_shared = []
    all_ = []
    for key, value in id_dict.items():

        if key != "all":
            sub_from_shared = sub_from_shared + value # add all the ids into one list

        else:
            all_ = id_dict[key]


    shared_ = list(set(all_) - set(sub_from_shared)) # get shared (non-species_specific) enh_ids

    n_dict["shared"] = len(shared_)
    id_dict["shared"] = shared_
    df.loc[df.enh_id.isin(shared_), "dataset"] = "shared" # rename the dataset as shared


    return n_dict, df


def plot_joint(x, y, data, RE):

    sns.set("talk")
    g = sns.jointplot(data = data, x= x,  y = y, kind = "hist", marginal_ticks = True)

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


def plot_dist_p(df, ndict):


    test = df.loc[(df.dataset != "all")& (df.dataset != "species_specific_10k") & (df.log_p<10)]

    x = "log_p"
    hue = "dataset"
    data = test
    col = "multiz"

    sns.displot(data = data, x = x, hue = hue, col = col, kind = "hist")
    plt.xlim(-5,5)
    #plt.title(f"chr21_{way}")
    #plt.legend(bbox_to_anchor =(1,1))

    outf = f"{RE}ch21_hist.pdf"
    plt.savefig(outf, bbox_inches = "tight")

    sns.displot(data = data, x = x, col = col, hue = hue, kind = "ecdf")
    plt.xlim(-5,5)
    #plt.legend(bbox_to_anchor =(1,1))
    outf = f"{RE}ch21_cdf.pdf"
    plt.savefig(outf, bbox_inches = "tight")

    medians = test.groupby(["multiz", "dataset"])["log_p"].median().reset_index()

    for way in test.multiz.unique():
        test2 = test.loc[test.multiz == way]
        shared = test2.loc[test2.dataset == "shared", "log_p"]
        hu_specific = test2.loc[test2.dataset == "hu_specific", "log_p"]
        result_stat, p = stats.mannwhitneyu(shared, hu_specific)
        print(way,"shared v. hu-specific", p)

        shared = test2.loc[test2.dataset == "shared", "log_p"]
        hu_specific = test2.loc[test2.dataset == "hu_specific/subtract_te", "log_p"]
        result_stat, p = stats.mannwhitneyu(shared, hu_specific)
        print(way,"shared v. hu-specific/subtract_te", p)
    print(medians)
    return medians

 #%%


df = pd.concat(fs.values())


ways = list(df.multiz.unique())


#table = pd.pivot(df, index = "id", columns = "multiz", values = "log_p")

#%%
df.head()
df.dataset.unique()
# split up df
#%%
n_dict, newdf = separate_shared_specific(df)

#%%


"""
#%% compare 20 v. 30 way
small_df = df[["enh_id", "multiz", "log_p"]].drop_duplicates() # drop redundant enh_id in different datasets.
table = small_df.pivot(index = "enh_id", columns = "multiz", values = "log_p")
table.head()
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
"""
#%% MWU p shared v. specific

plot_dist_p(newdf, n_dict)
