#!/usr/bin/env python

# intersect enhancers with TEs and measure architecture fold enrichment per age.
# hg19 genome build required. 

import argparse
import datetime
import glob

from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

import numpy as np
import os, sys
import pandas as pd
import pybedtools as pb
import seaborn as sns
from scipy import stats
from statsmodels.stats import multitest
import subprocess

"""
arg_parser = argparse.ArgumentParser(description=" describe argparse")

arg_parser.add_argument("enhancers", help='bed file w/ full path')

args = arg_parser.parse_args()

ENHF = args.enhancers

ENHPATH = "/".join(ENHF.split("/")[:-1])
SID = (ENHF.split("/")[-1]).split(".bed")[0]
"""


ENHPATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/"
ENHF = f"{ENHPATH}trimmed-310_all_unique_fantom_erna_112_tissue.bed"
SID = (ENHF.split("/")[-1]).split(".")[0]


#%% function to make directories for outputs


def make_data_path(path):
    if os.path.exists(path) == False:
        os.mkdir(path)


#%% Paths and files to load.

# directory to stash intersection data
OUTPATH = os.path.join(ENHPATH, "te")
make_data_path(OUTPATH)

# directory to stash results
RE = os.path.join(OUTPATH, "results")
make_data_path(RE)

# path to all the TE files
TE_PATH = "/dors/capra_lab/users/abraha1/projects/transposable_elements/data/"
TE_FAM = f"{TE_PATH}hg19_TE_counts_wlin.txt"
DFAM_FILE = f"{TE_PATH}dfam/hg38_dfam.nrph.hits_genomicCoordinates_tidy.bed"
REPEATMASKER = f"{TE_PATH}filtered_formated_hg19fromhg38-TE_coords.tsv"


#%% FUNCTIONS


def get_syn_gen_bkgd():

    SYN_GEN_BKGD_FILE = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_taxon.bed"
    syn_gen_bkgd= pd.read_csv(SYN_GEN_BKGD_FILE, sep = '\t') # read the file
    syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages


    # handle the TEs that have an old age
    old = pd.DataFrame({"mrca": [0.39], "taxon": ["older_than_mammalia"], "mrca_2": [0.39],
                            "taxon2": ["older_than_mammalia"]})

    syn_gen_bkgd = pd.concat([syn_gen_bkgd, old])

    return syn_gen_bkgd


def format_te_file(te_fam, syn_gen_bkgd):
    # FORMAT TE FAMILY FILE to match syntenic background

    # open file as df
    famdf = pd.read_csv(te_fam, sep = '\t')

    # rename columns
    famdf.columns = ["te_fam", "fam", "taxon", "count"]

    # make taxon annotations consistent w/ my syn_gen_bkgd file
    famdf.taxon.loc[famdf.taxon == "Primates"] = "Primate"
    famdf.taxon.loc[famdf.taxon == "theria"] = "Theria"
    famdf.taxon.loc[famdf.taxon == "Euarchonta"] = "Euarchontoglires"
    famdf.taxon.loc[famdf.taxon == "Hominoidea"] = "Hominidae"
    famdf.taxon.loc[famdf.taxon == "Rodentia"] = "Euarchontoglires"
    famdf.taxon.loc[famdf.taxon == "Muridae"] = "Euarchontoglires"
    famdf.taxon.loc[famdf.taxon == "Homo_sapiens"] = "Homo sapiens"
    famdf.taxon.loc[famdf.taxon == "old"] = "older_than_mammalia"

    # merge w/ age annotations
    famdf = pd.merge(famdf, syn_gen_bkgd, how = "left", on = "taxon") # get ages

    # rename columns again.
    famdf.columns = ["te", "fam", "taxon_te", "count_te",  "mrca_te",  "mrca_2_te",  "taxon2_te"]

    return famdf


def bed_intersection(enhf, repeatmasker, outpath, sid):


    outf = "%s%s_x_te.bed" % (outpath, sid)

    if os.path.exists(outf) == False: # have you done the intersection?

        f = pb.BedTool(enhf).sort()
        r = pb.BedTool(repeatmasker).sort()

        f.intersect(r, wao = True).saveas(outf)

    return outf


def cleandf(outpath, outf):

    # make sure file is tab separated.
    outid = (outf.split("/")[-1]).split(".")[0]
    newf = "%scleaned_%s.bed"% (outpath, outid)
    cmd = '''awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' %s > %s''' % (outf, newf)

    subprocess.call(cmd, shell=True)
    return newf


def format_df(intersected_f, sid):

    df = pd.read_csv(intersected_f, sep = '\t', header = None)


    df.columns = ["chr_enh", "start_enh", "end_enh", "old_enh_len", "old_enh_id",
                 "mrca", "core_remodeling",
                 "chr_te", "start_te", "end_te", "te_fam", "len_te"]


    # add annotations
    df["enh_id"] = df.chr_enh + ":" + df.end_enh.map(str) +"-"+ df.start_enh.map(str) # enh id
    df["enh_len"] = df["end_enh"] - df["start_enh"] # enh length

    df["arch"] = "" # code the enhancers
    df["arch"].loc[df["core_remodeling"] ==0] = "simple"
    df["arch"].loc[df["core_remodeling"] ==1] = "complexenh"
    df["dataset"] = sid

    # cleanup
    df = df.loc[df.old_enh_len >= 6]
    df = df.loc[df.chr_enh != "chrX"]
    df = df.drop_duplicates() # drop duplicate
    df = df.loc[df["core_remodeling"] != 0.175] # idk what this is, but excluding...

    # formatting
    df = df.replace(".", 0) # replace all ".", with 0's (when TE doesn't overlap enh)
    df.len_te = df.len_te.fillna(0)
    df.len_te = df.len_te.astype(int) # change the datatype
    df.mrca = df.mrca.round(3) # round the ages
    df.core_remodeling = df.core_remodeling.astype(float) # change the datatype.

    # add summarized MRCA values.
    df = pd.merge(df, syn_gen_bkgd, how = "left", on = "mrca")

    return df


def format_te_overlaps(df):

    # create a binary column for TE overlaps
    df["te_bin"] = 0 # te overlap binary
    df.te_bin.loc[df["len_te"]>= 6] = 1 # re-write all the len_te that are less than 10bp long.

    # format the TE, TE family annotations
    df["fam"] = df.te_fam.loc[df.te_fam != 0].apply(lambda x: (x.split("_")[-1]))
    df["te"] = df.te_fam.loc[df.te_fam != 0].apply(lambda x: "".join(x.split("_")[0]))

    return df


def merge_enh_te_info(enhdf, famdf):

    # save enhancer age info to re-merge after merging enh, te info
    enh_mrca = enhdf[["enh_id","mrca", "mrca_2", "taxon2"]].drop_duplicates()
    enh_mrca.columns = ["enh_id", "core_mrca", "core_mrca_2", "taxon2"]

    # merge enhancer+TE intersection with TE family information

    enh_te = pd.merge(enhdf[["enh_id", "old_enh_id", "te", "fam",
                             "te_fam","core_remodeling",
                             "te_bin", "len_te", "mrca", "mrca_2"]],
                             famdf, how = "left")\
                             .drop_duplicates()

    # merge enh+TE information with enhancer ages information
    enh_te = pd.merge(enh_te, enh_mrca, how = "left", on = "enh_id")

    # annotate when the TE origin is older, younger, or same age as enhancer origin
    enh_te["class2"] ="No TE"
    enh_te["class2"].loc[(enh_te.core_mrca_2 == enh_te.mrca_2) &(enh_te.te_bin ==1)] = "TE same age"
    enh_te["class2"].loc[(enh_te.core_mrca_2<enh_te.mrca_2)  &(enh_te.te_bin ==1)] = "TE older"
    enh_te["class2"].loc[(enh_te.core_mrca_2>enh_te.mrca_2)  &(enh_te.te_bin ==1)] = "TE younger"


    return enh_te


def get_numbers(baseFam):

    # total TE that overlap enhancer
    num_te = len(baseFam.loc[baseFam.te_bin==1]) # 16454 TEs overlap 13277 enhancers
    print("number of TEs overlapping enhancers", num_te)

    # total enhancers that overlap TE
    num_enh_w_te = len(baseFam.loc[baseFam.te_bin==1].enh_id.unique()) #13277/30434 enhancers overlap TEs
    print("number of enhancers w/ overlapping TEs", num_enh_w_te)

    num_enh = len(baseFam.enh_id.unique())
    print("number of enhancers", num_enh)


def prep_df_for_OR(enh_te):

    # get enhancer_id, arch, age, and TE overlap (some enhancers overlap multiple TEs, so we need to reduce the datafram )
    base = enh_te.groupby(["enh_id", "core_remodeling", "mrca_2"])["te_bin"].max().reset_index()

    # merge enh, arch, oldest age, bin w/ TE families per fam
    # left join
    baseFam = pd.merge(base, enh_te[["enh_id", "fam"]].drop_duplicates(), how = "left")

    get_numbers(baseFam) # count te, enh ooverlaps

    baseTE = baseFam.loc[baseFam.te_bin==1]

    return baseTE, baseFam


def fdr_formatting(results_dict):

    df = pd.concat(results_dict.values()) # concat results dict

    fdr_bool, fdr_pval =  multitest.fdrcorrection(df.pvalue, alpha=0.10, method='indep', is_sorted=False) # FDR 5%

    df["fdr_pval"] = fdr_pval # add FDR to df
    df["-log10p"] = np.log10(df["fdr_pval"])*-1 # -log10(p)

    df["log2_odds"]= np.log2(df["odds"]) # log2(OR)

    return df


def get_mrca_OR(enh_te, baseTE):

    all_TE_mrca_dict = {} # collect all te/mrca results here.

    fams = baseTE.fam.unique() # get list of TE families
    mrcas = baseTE.mrca_2.unique() # list of mrcas

    for mrca in mrcas: #

        mrca_dict = {}

        for fam in fams[1:]:

            test = baseTE.loc[baseTE.mrca_2 == mrca]

            simple_Yfam_Ymrca = len(test.loc[(test.core_remodeling ==0) & (test.fam == fam)])
            simple_Nfam_Ymrca = len(test.loc[(test.core_remodeling ==0) & (test.fam != fam)])

            complex_Yfam_Ymrca = len(test.loc[(test.core_remodeling ==1) & (test.fam == fam)])
            complex_Nfam_Ymrca = len(test.loc[(test.core_remodeling ==1) & (test.fam != fam)])

            # set up 2x2 table
            a,b = complex_Yfam_Ymrca, simple_Yfam_Ymrca
            c,d = complex_Nfam_Ymrca, simple_Nfam_Ymrca

            obs = [[a,b], [c,d]]

            if a >= 5 or b >= 5:

                odds, pvalue = stats.fisher_exact(obs)

                df = pd.DataFrame({"fam": [fam],
                "odds": [odds],
                "pvalue":[pvalue],
                "test_simple": [simple_Yfam_Ymrca],
                "test_complex": [complex_Yfam_Ymrca],
                "test_simple_notFam":[simple_Nfam_Ymrca],
                "test_complex_notFam":[complex_Nfam_Ymrca],
                "mrca_2":[mrca]})

                mrca_dict[fam]=df # collect results per fam per age

        fdr_df = fdr_formatting(mrca_dict) # FDR correction at 10% FDR

        all_TE_mrca_dict[mrca] = fdr_df # collect results of all fams (fdr_df) per age

    all_TE_mrca_OR = pd.concat(all_TE_mrca_dict.values()) # concat all results across all ages

    return all_TE_mrca_OR


def format_OR_results(ORdf, syn_gen_bkgd):


    # add taxon2 info
    ORdf = pd.merge(ORdf, syn_gen_bkgd[["mrca_2", "taxon2"]], how = "left", on = "mrca_2").drop_duplicates()

    # MANIPULATION: replace -log10p for cases where p = 0 (super significant)
    ORdf.loc[abs(ORdf["fdr_pval"]) == 0, "-log10p"] = 1 #

    # MANIPULATION: articifically set odds to 2 when odds are positively inf
    # (when TE is ONLY in complex, but NEVER in SIMPLE)

    ORdf.loc[ORdf["log2_odds"] == np.inf, "log2_odds"] = 2

    return ORdf


def get_pval_table_for_plotting(ORdf_formatted):

    # pivot - rows = fams, ages = cols, values = -log10p
    pval = ORdf_formatted.sort_values(by = "mrca_2").pivot(index = "fam", columns = "mrca_2", values="-log10p")

    # fill in all the insignificant log10 pvals (ie the log10p vales that are 0 or None)
    pval = pval.fillna(0)

    # drop the first column, cannot compare human simple v. complex bc no human complex enhancers
    pval = pval.drop([0.0], axis = 1)

    # like pval pivot, but for log2_odds as table values
    # pivot - rows = fams, ages = cols, values = log2_odds
    logodds = ORdf_formatted.sort_values(by = "mrca_2").pivot(index = "fam", columns = "mrca_2", values="log2_odds")
    logodds = logodds.drop([0.0], axis = 1)

    # drop all the TE families with no differences in log odds of simple/complex
    #in at least 5 age categories.
    logodds = logodds.replace(0, np.nan) # replace all zero odds with None
    logodds = logodds.dropna(thresh = 5) # drop all rows where TE families have more than 5 logodds = zero values.
    logodds = logodds.fillna(0) # refill smaller table the missing values w/ zeros


    pval = pval.loc[pval.index.isin(logodds.index)] # match TE fams in pval table to logodds table
    pval_asterisks = pval.apply(lambda x: ["*" if y >=1 else "" for y in x]) # asterisks dataframe

    return pval, logodds, pval_asterisks


def plot_heatmaps(ORdf_formatted,log2odds, pval, pval_asterisks, outfile):

    # prepare to plot tables
    cm = log2odds.fillna(0)
    cmp = pval.fillna(0)

    xtick = ["Prim", "Euar", "Bore", "Euth", "Ther", "Mam", "Amni", "Tert", "Vert"]

    sns.set(rc={'axes.facecolor':'grey', 'figure.facecolor':'white'})
    sns.set(font_scale=1.5)

    # first, make the dendrogram to get the clustering of TE families.
    cg = sns.clustermap(cm,
                mask = cm == 0, # mask all the 0 values.
                col_cluster = False, center = 0,
                annot = True,
                annot_kws={"size": 16},
                fmt = "0.1f",
                cmap="YlGn", robust = True,
                xticklabels = xtick,
                linewidths=0.5,figsize = (14,14))

    #plt.savefig(f"{RE}te_fam_cluster_enrichment.pdf", bbox_inches = 'tight')

    # reindex tables by dendrogram
    new_index = cg.dendrogram_row.reordered_ind

    log2odds = log2odds.reset_index().reindex(new_index).set_index("fam") # reindex pvals by dendrogram
    pval = pval.reindex(log2odds.index) # reindex pvals by dendrogram
    pval_asterisks = pval_asterisks.reindex(log2odds.index) # reindex pvals by dendrogram

    pval = pval.fillna(0)

    labels =  pval_asterisks


    # Create a custom diverging colormap
    amber = '#feb308'
    faded_green = '#7bb274'
    offwhite = '#ffffe4'

    a = mcolors.to_rgba(amber)
    f = mcolors.to_rgba(faded_green)
    o = mcolors.to_rgba(offwhite)

    colors = [a, o, f]  # R -> G -> B
    cmap_name = 'my_list'
    cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=128)

    # set other plotting parameters
    grid_kws = {"height_ratios": (.9, .1), "hspace":1.1}

    # plot pvalue asterisks
    f, (ax, cbar_ax) = plt.subplots(2, gridspec_kw=grid_kws, figsize=(7,14))
    ax = sns.heatmap(pval,
        mask = pval < 1 ,
        ax = ax,
        annot = labels,
        annot_kws={"size": 30, "color": "black", "fontweight":'bold'},
        fmt = '',
        cbar= False)

    # plot log2 odds colors
    ax = sns.heatmap(log2odds,
        center =0,
        ax=ax,
        mask = log2odds ==0,
        cbar_ax = cbar_ax,
        robust = True,
        cbar_kws={"orientation": "horizontal",
        "label": "Log2-fold enrichment\nSimple Complex"},
        cmap = cm,
        xticklabels = xtick,
        linewidths=0.5)


    ax.set(xlabel = "", ylabel = "")

    plt.savefig(outfile, bbox_inches = 'tight')


#%% run functions


syn_gen_bkgd = get_syn_gen_bkgd() # get summarized age information
syn_gen_bkgd.head()
famdf = format_te_file(TE_FAM, syn_gen_bkgd) # get TE family dataframe w/ages

outf = bed_intersection(ENHF, REPEATMASKER, OUTPATH, SID) # intersect enhancers w/ TE
newf = cleandf(OUTPATH, outf) # clean up intersection file

#%% Format enhdf

enhdf = format_df(newf, SID)
enhdf = format_te_overlaps(enhdf)

enh_te = merge_enh_te_info(enhdf, famdf)


### only test enhancers that overlap TEs ###


# calculate enrichment of TE in complex, simple arch per mrca
baseTE, baseFam = prep_df_for_OR(enh_te)

ORdf = get_mrca_OR(enh_te, baseTE)

ORdf_formatted = format_OR_results(ORdf, syn_gen_bkgd)

pval, log2odds, pval_asterisks  = get_pval_table_for_plotting(ORdf_formatted)

outfile = f"{RE}/fig5d-te_fam_cluster_enrichment_pval.pdf"

plot_heatmaps(ORdf_formatted,log2odds, pval, pval_asterisks, outfile)


#%%
### test all enhancers, regardless whether they overlap TEs ###


ORdf = get_mrca_OR(enh_te, baseFam)

ORdf_formatted = format_OR_results(ORdf, syn_gen_bkgd)

pval, log2odds, pval_asterisks  = get_pval_table_for_plotting(ORdf_formatted)

outfile = f"{RE}/fig5d-all_enh_te_fam_cluster_enrichment_pval.pdf"

plot_heatmaps(ORdf_formatted,log2odds, pval, pval_asterisks, outfile)
