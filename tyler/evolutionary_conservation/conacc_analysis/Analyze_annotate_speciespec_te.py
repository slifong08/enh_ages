import argparse
import matplotlib.pyplot as plt
import numpy as np
import os, sys
import pandas as pd
from scipy import stats
import seaborn as sns
import subprocess

###
# goal
###

# take all estimates of conservation and acceleration
# from tyler's atac-seq dataset.
#
# and annotate with
# hu-specific/rhe-specific accessibility
# hg38 TE overlap from repeatmasker.


#%% ARGPARSE arguments.
"""
arg_parser = argparse.ArgumentParser(description=" describe argparse")

arg_parser.add_argument("bedfile", help ='bed file w/ full path')
arg_parser.add_argument("-br","--branches", help ='hg38, rheMac3')
arg_parser.add_argument("-m", "--multiz", help ='20-, 30-, 100-way multiple sequence alignments in hg38')


# PARSE THE ARGUMENTS
args = arg_parser.parse_args()

F = args.bedfile # the bedfile
BRANCH = args.branches # the branches to test.
PATH = "/".join(F.split("/")[:-1]) + "/" # the path
MSA_WAY = args.multiz # multiple sequence alignment.
"""

branches = ["hg38-rheMac8", "hg38", "rheMac8"]
#BRANCH = "hg38-rheMac8"

MSA_WAY = 30

TE_F = "/dors/capra_lab/data/transposable_elements/repeatmasker/hg38.bed"

SPECIES_DA_F = "/dors/capra_lab/users/fongsl/tyler/data/GG-LL_species-specific_OCRs_rank/GG-LL_species-specific_OCRs_rank.bed"

SELF_F = "/dors/capra_lab/data/ucsc/hg38/self/hg38_repeats_self_coor.bed"


#%% FUNCTIONS

def get_pathVars(branch, msa_way):

    conacc_f = f"/dors/capra_lab/users/fongsl/tyler/data/CON_ACC/all/multiz30way_{branch}/all_con_acc.bed"
    re = f"/dors/capra_lab/projects/enhancer_ages/tyler/results/CON_ACC/{msa_way}way_{branch}/"
    path = "/".join(conacc_f.split("/")[:-1]) + "/" # the path

    if os.path.exists(re) == False:
        os.mkdir(re)

    return conacc_f, re, path


def remove_doubletabs(f, path):
    os.chdir(path) # go to the directory
    temp = "temp.bed" # temp file
    cmd = f"tr -s '\t' < {f} > {temp} && mv {temp} {f}"
    subprocess.call(cmd, shell = True)


def format_f(conacc_f):

    cols = ["#chr", "b", "bf", "start", "end", "conacc", "?", "??", "id"]
    test = pd.read_csv(conacc_f, sep  = '\t', nrows = 5)

    if "start" not in list(test): # do you need to rename the columns?
        df = pd.read_csv(conacc_f, sep  = '\t', header = None, names = cols)
        df = df.drop(["b", "bf", "?", "??"], axis = 1) # drop the columns you don't need
        df.to_csv(conacc_f, sep = '\t', index = False) # save this formatted file.

    else: # or you've renamed them already and just need to read the df.
        df = pd.read_csv(conacc_f, sep  = '\t')

    return df


def make_id(df, chr, start, end, idname):

    df[idname] = df[chr] + ":" + df[start].map(str) + "-" + df[end].map(str)

    return df


def test_colnames_add_id(outf, cols, chr, start, end, idname):

    # test to see if the columns have already been named.
    test = pd.read_csv(outf, sep = '\t', nrows = 5)

    if "start" not in list(test):

        # open df
        df = pd.read_csv(outf, sep = '\t', header = None, )
        df.columns = cols

        # assign new id

        df = make_id(df, chr, start, end, idname)

        df = df.drop([chr, start, end], axis = 1) # drop the extra cols

        df.to_csv(outf, sep = '\t', index = False) # save the modified file.


def intersect_te(conacc_f, path, te):


    #te = pd.read_csv(TE, sep = '\t', skiprows = 1)

    # get info to write new files
    sample_id = (conacc_f.split("/")[-1]).split(".bed")[0]
    outTE = f"{path}{sample_id}_TE.bed"

    if os.path.exists(outTE) == False: # if you haven't done this intersection already...

        cmd = f"bedtools intersect -a {conacc_f} -b {te} -wao > {outTE}"
        subprocess.call(cmd, shell = True)

        remove_doubletabs(outTE, path)

    # format the output file

    # col names
    cols = ["#chr", "start", "end", "CON_ACC", "id",
    "chr_te", "start_te", "end_te", "te", "te_fam", "len_te_overlap"]

    # add id for te locus
    idname = "te_id"

    #specify locus
    chr, start, end = "chr_te", "start_te", "end_te"

    # add column names, id, drop chr, start, end, and save the new formatted file.
    test_colnames_add_id(outTE, cols, chr, start, end, idname)

    return outTE


def intersect_species_specific(conacc_f, path, species_da_f):

    sample_id = (conacc_f.split("/")[-1]).split(".bed")[0]
    outSP = f"{path}{sample_id}_species.bed"

    if os.path.exists(outSP) == False:
        # intersect dataframe with species-specific classes.
        cmd = f"bedtools intersect -a {conacc_f} -b {species_da_f} -wao > {outSP}"
        subprocess.call(cmd, shell = True)
        remove_doubletabs(outSP, path) # remove tabs

    # format the output dataframe
    # col names
    cols = ['#chr', 'start', 'end', 'CON_ACC', 'id',
    "chr_sp", "start_sp", "end_sp", "species", "spoverlap"]

    # add id for sp locus
    idname = "sp_id"

    #specify locus
    chr, start, end = "chr_sp", "start_sp", "end_sp"

    # add column names, id, drop chr, start, end, and save the new formatted file.
    test_colnames_add_id(outSP, cols, chr, start, end, idname)

    return outSP


def intersect_self(conacc_f, path, self_f):

    sample_id = (conacc_f.split("/")[-1]).split(".bed")[0]
    outSLF = f"{path}{sample_id}_self.bed"

    if os.path.exists(outSLF) == False:
        cmd = f"bedtools intersect -a {conacc_f} -b {self_f} -wao > {outSLF}"
        subprocess.call(cmd, shell = True)
        print(cmd)
        remove_doubletabs(outSLF, path) # remove tabs

    # format the output dataframe
    # col names
    cols = ['#chr', 'start', 'end', 'CON_ACC', 'id',
    "chr_slf", "start_slf", "end_slf", "slfOverlap"]

    # add id for te locus
    idname = "slf_id"

    #specify locus
    chr, start, end = "chr_slf", "start_slf", "end_slf"

    # add column names, id, drop chr, start, end, and save the new formatted file.
    test_colnames_add_id(outSLF, cols, chr, start, end, idname)

    return outSLF


def compile_df(outte, outsp, outslf):
    dfte = pd.read_csv(outte, sep = '\t')
    dfsp = pd.read_csv(outsp, sep = '\t')
    dfslf = pd.read_csv(outslf, sep = '\t')

    df = pd.merge(dfte, dfsp, how = "left")
    df = pd.merge(df, dfslf, how = "left")

    # drop cols
    df = df.drop(["spoverlap", "slf_id", "sp_id"], axis = 1).drop_duplicates()

    # te_bin
    # make a binary for TE overlap

    te = df.groupby("id")["len_te_overlap"].sum().reset_index()

    te["te_bin"] = 0

    # assign te_bin ==1 to all enhancer where TEs overlap 1+ bp
    te.loc[te["len_te_overlap"] >0, "te_bin"] = 1
    te = te.drop(["len_te_overlap"], axis = 1)

    # drop the TE cols, we don't need them
    drop_TE_cols = ['te', 'te_fam', 'len_te_overlap', 'te_id']
    df = df.drop(drop_TE_cols, axis = 1).drop_duplicates()
    df = pd.merge(df, te, how = "left")

    # fill in the shared elements
    df.loc[df["species"] == ".", "species"] = "shared"

    # create a
    df["sp_te"] = df["species"] + "-" + df["te_bin"].map(str)

    # sometimes, an element can overlap more than one self region.
    # sum these multi-overlapping rows together.
    df = df.groupby(['#chr', 'start', 'end','CON_ACC',
    'id', 'species', 'te_bin','sp_te'])["slfOverlap"].sum().reset_index()

    # create a binary for self-chain overlap
    df["slf_bin"] = 0
    df.loc[df["slfOverlap"]>0, "slf_bin"] = 1
    # drop the self overlap length column
    df.drop(["slfOverlap"], axis = 1)

    df["sp_te_slf"] = df["species"] + "-" + df["te_bin"].map(str) + "-" +  df["slf_bin"].map(str)
    df["te_slf"] = df["te_bin"].map(str) +"-" +  df["slf_bin"].map(str)
    df["sp_slf"] = df["species"] + "-" +  df["slf_bin"].map(str)

    return  df


def plot_cdf(x, data, hue, title, pal):
    sns.set("talk")
    n = data.groupby([hue])["id"].count().reset_index()

    g = sns.displot(data, x=x, hue = hue,  kind ="ecdf", palette = pal)
    g.set( title = title,
    xlim = (-0.5,0.5),
    ylim = (0.05, 0.70),
    xlabel = f"CON_ACC\n{n}"
    )

    outf = f"{RE}{title}_{MSA_WAY}way_{BRANCH}.pdf"
    plt.savefig(outf, bbox_inches = "tight")

#%% run functions
PATH
#%%
for BRANCH in branches:
    print(BRANCH)
    CONACC_F, RE, PATH = get_pathVars(BRANCH, MSA_WAY)

    df = format_f(CONACC_F) # format the file

    outTE = intersect_te(CONACC_F, PATH, TE_F) # intersect w/ repeatmasker TE

    outSP = intersect_species_specific(CONACC_F, PATH, SPECIES_DA_F) # intersect species specific accessibility calls.

    outSLF = intersect_self(CONACC_F, PATH, SELF_F) # intersect self accessibility calls.

    #
    newdf = compile_df(outTE, outSP, outSLF)
    newdf.info()
    #

    (newdf.groupby("species")["id"].count())

    # TE stratified
    x = "CON_ACC"
    data = newdf.sort_values(by = "sp_te")
    hue = "te_bin"
    title = "TE_stratified"
    pal = "Set1"
    plot_cdf(x, data, hue, title, pal)

    # everything stratified
    hue = "sp_te"
    title = "all"
    pal = "tab20"
    data = newdf.sort_values(by = "sp_te")
    plot_cdf(x, data, hue, title, pal)

    # no TE stratified

    data = newdf.loc[newdf.te_bin ==0]
    title = "TE_excluded"
    pal = "Set1"
    plot_cdf(x, data, hue, title, pal)


    # TE only stratified

    data = newdf.loc[newdf.te_bin ==1]
    title = "TE_only"
    pal = "Set1"
    plot_cdf(x, data, hue, title, pal)


    # Species stratified
    data = newdf
    title = "Species_only"
    hue = "species"
    pal = "Set1"
    plot_cdf(x, data, hue, title, pal)

    # Self stratified
    data = newdf
    title = "Self"
    hue = "slf_bin"
    pal = "tab10"
    plot_cdf(x, data, hue, title, pal)

    # Species, self, and TE stratified
    data = newdf.sort_values(by = "sp_te_slf")
    title = "Sp_Te_Self"
    hue = "sp_te_slf"
    pal = "tab20c"
    plot_cdf(x, data, hue, title, pal)

    # self, and TE stratified
    data = newdf.sort_values(by = "te_slf")
    title = "Te_Self"
    hue = "te_slf"
    pal = "tab10"
    plot_cdf(x, data, hue, title, pal)

    #Species, self stratified
    data = newdf.sort_values(by = "sp_slf")
    title = "Sp_Self"
    hue = "sp_slf"
    pal = "tab10"
    plot_cdf(x, data, hue, title, pal)

    # self, te stratified
    TeSlf = newdf.loc[
    (newdf.te_bin == 1) &
    (newdf.slf_bin == 1)
    ].drop_duplicates()

    TeSlf.describe()

    # no self, no te stratified
    noTeSlf = newdf.loc[
    (newdf.te_bin == 0) &
    (newdf.slf_bin == 0)
    ].drop_duplicates()
    noTeSlf.quantile(0.01)

    noTeSlf.describe()

    con_noSlfTe = noTeSlf.CON_ACC
    con_SlfTe = TeSlf.CON_ACC

    fig, ax =plt.subplots()

    sns.histplot(con_noSlfTe, ax = ax, label = "no_slf_te", color = "r")
    sns.histplot(con_SlfTe, ax = ax, label = "slf_te")
    ax.legend()

    print("self v. non-self", BRANCH)
    print(stats.mannwhitneyu(con_noSlfTe, con_SlfTe))
    print(con_noSlfTe.median(),con_SlfTe.median())
