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

        fs[key] = df

#%%

df = pd.concat(fs.values())

table = pd.pivot(df, index = "id", columns = "multiz", values = "log_p")

#%%


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
