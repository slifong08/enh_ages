import matplotlib.pyplot as plt
import numpy as np
import os, sys
import pybedtools as pb
import pandas as pd
import subprocess
import seaborn as sns
from scipy import stats
from sklearn.preprocessing import StandardScaler



AGEPATH = "/dors/capra_lab/projects/enhancer_ages/kircher19/GRCh37_ALL_elements/ages/"
AGEFILE = "syn_breaks_GRCh37_ALL_elements_ages.bed"
AGEF = os.path.join(AGEPATH, AGEFILE)

ACTIVITYPATH = "/dors/capra_lab/projects/enhancer_ages/kircher19"
ACTIVITYFILE = "no_header_reduced.bed"
ACTIVITYF = os.path.join(ACTIVITYPATH, ACTIVITYFILE)

OUTFILE = "GRch37_ALL_x_syn.bed"

OUTF = os.path.join(ACTIVITYPATH, OUTFILE)

RE = os.path.join(ACTIVITYPATH, "results")

#%% FUNCTIONS


def add_arch_labels(df):

    df["arch"] = "complex_core"

    df.loc[df.core_remodeling ==0, "arch"] = "simple"
    df.loc[df.core ==0, "arch"] = "complex_derived"

    return df

def add_annot_labels(df):

    promoters = ["F9", 'FOXE1', "GP1BB", "HBB", "HNF4A", "LDLR", "MSMB",
    "PKLR",'HBG1', "TERT", 'GP1BA',
    ] #'LDLR.2','TERT-GAa', 'TERT-HEK', 'TERT-GSc', 'TERT-GBM','PKLR-24h', 'PKLR-48h',

    enhancers = ["BCL11A", "IRF4", "IRF6", "MYC", "RET", "SORT1",
     "TCF7L2", "UC88", "ZFAND3", "ZRS",
     'MYCrs6983267', 'MYCrs11986220'] # "SORT1-flip", 'SORT1.2','ZRSh-13', 'ZRSh-13h2',

    df["annot"] = ""
    df.loc[df.Element.isin(promoters), "annot"] = "promoter"
    df.loc[df.Element.isin(enhancers), "annot"] = "enh"

    df["annot2"] = df.annot + "-" + df.arch

    return df


def format_df(df):

    cols = ["chr", "start", "end", "Value", "pval",
    "Element", "chr_syn", "start_syn", "end_syn",
    "enh_id", "chr_enh", "start_enh", "end_enh",
    "seg_index", "core_remodeling", "core", "mrca"]

    df.columns = cols

    df["abs_val"] = df["Value"].abs() # get absolute value of change
    df["syn_len"] =  df.end_syn -df.start_syn # get syn length

    df = df.loc[df.syn_len>=6] # ignore any values where syn length is less than 6 bp long.


    # control for the fact that length increases random probability that a
    # variant will have high change in activity.


    df["abs_val_persynlen"] = df["abs_val"].divide(df["syn_len"]) # divide abs val by syn length
    df["val_persynlen"] = df["Value"].divide(df["syn_len"]) # divide val by syn length

    return df

def standard_scaler(df, Element):

    values = ["Value", "abs_val", "val_persynlen", "abs_val_persynlen"]

    test = df.loc[df.Element == Element, values]

    for col in values:
        col_zscore = col + '_zscore'
        x = test[col]
        test[col_zscore] = stats.zscore(x)

    df = pd.merge(df, test, how = "left")
    return df


#%%


if os.path.exists(OUTF) == False:
    act = pb.BedTool(ACTIVITYF).sort()
    age = pb.BedTool(AGEF).sort()

    act.intersect(age, wa = True, wb = True).saveas(OUTF)
#%%



df = pd.read_csv(OUTF, sep = '\t', header = None)
df = format_df(df)
df = add_arch_labels(df)
df = add_annot_labels(df)
list(df)
for e in df.Element.unique():
    df = standard_scaler(df, e)
df.head()

#%%

annotations = ['enh', 'promoter']


for annot in df.annot.unique():
    if annot in annotations:

        order = ["simple", "complex_core", "complex_derived"]
        print(annot)
        x = "Element"
        hue = "arch"
        data = df.loc[df.annot == annot]

        fig, (ax1, ax2, ax3, ax4) = plt.subplots( ncols = 4, figsize = (24,12))
        sns.set("poster")

        y = "Value_zscore"
        sns.boxplot( y=x, x=y, data = data, hue = hue, hue_order = order,
         showfliers = False, notch = True, ax = ax1)
        ax1.legend().remove()
        ax1.axvline(0)

        y = "abs_val_zscore"
        sns.barplot( y=x, x=y, data = data, hue = hue,hue_order = order, ax = ax2,
        estimator = np.median)
        ax2.legend().remove()
        ax2.set_yticklabels("")

        y = "val_persynlen_zscore"
        sns.boxplot( y=x, x=y, data = data, hue = hue, hue_order = order,
         showfliers = False, notch = True, ax = ax3)
        ax3.legend().remove()
        ax3.axvline(0)
        ax3.set_yticklabels("")

        y = "abs_val_persynlen_zscore"
        sns.barplot( y=x, x=y, data = data, hue = hue,hue_order = order, ax = ax4,
        estimator = np.median)
        ax4.legend(bbox_to_anchor = (1,1))
        ax4.set(title = f"{annot}")
        ax4.set_yticklabels("")

        outf = f"{RE}/{annot}.pdf"
        plt.savefig(outf, bbox_inches = "tight")
