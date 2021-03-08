import matplotlib.pyplot as plt
import numpy as np
import os, sys
import pybedtools as pb
import pandas as pd
import subprocess
import seaborn as sns
from scipy import stats

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

def bed_cmd(a, b, outf):

    cmd = f"bedtools intersect -a {a} -b {b} -wao > {outf}"

    subprocess.call(cmd, shell = True)
    print(cmd)

    return outf


def add_arch_labels(df):

    df["arch"] = "complex_core"

    df.loc[df.core_remodeling ==0, "arch"] = "simple"
    df.loc[df.core ==0, "arch"] = "complex_derived"

    return df

def  add_annot_labels(df):

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
#%%

act = pb.BedTool(ACTIVITYF).sort()
age = pb.BedTool(AGEF).sort()

act.intersect(age, wa = True, wb = True).saveas(OUTF)
#%%
cols = ["chr", "start", "end", "Value", "pval",
"Element", "chr_syn", "start_syn", "end_syn",
"enh_id", "chr_enh", "start_enh", "end_enh",
"seg_index", "core_remodeling", "core", "mrca"]


df = pd.read_csv(OUTF, sep = '\t', header = None)
df.columns = cols
df = add_arch_labels(df)
df = add_annot_labels(df)

df["abs_val"] = df["Value"].abs()
df.head()
#%%
for annot in df.annot.unique():
    order = ["simple", "complex_core", "complex_derived"]
    print(annot)
    x = "Element"
    y = "Value"
    hue = "arch"
    data = df.loc[df.annot == annot]

    fig, (ax1, ax2) = plt.subplots( ncols = 2, figsize = (12,12))
    sns.set("poster")
    sns.boxplot( y=x, x=y, data = data, hue = hue, hue_order = order,
     showfliers = False, notch = True, ax = ax1)

    ax1.legend().remove()
    #ax1.title(f"{annot}")
    ax1.axvline(0)

    y = "abs_val"
    sns.barplot( y=x, x=y, data = data, hue = hue,hue_order = order, ax = ax2,
    estimator = np.median)
    ax2.legend(bbox_to_anchor = (1,1))
    ax2.set(title = f"{annot}")
    ax2.set_yticklabels("")
    outf = f"{RE}/{annot}.pdf"
    plt.savefig(outf, bbox_inches = "tight")
#%%
data
