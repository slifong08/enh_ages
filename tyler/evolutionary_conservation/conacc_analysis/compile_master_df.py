import argparse
import matplotlib.pyplot as plt
import numpy as np
import os, sys
import pandas as pd
from scipy import stats
import seaborn as sns
import statsmodels
import statsmodels.api as sm
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
# mappability from liftover
# sequence identity from biopython


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
branches = ["hg38", "rheMac8", "hg38-rheMac8"]
PRUNE = True
FDR = True
MSA_WAY = 30


TE_F = "/dors/capra_lab/data/transposable_elements/repeatmasker/hg38.bed"

SPECIES_DA_F = "/dors/capra_lab/users/fongsl/tyler/data/GG-LL_species-specific_OCRs_rank/GG-LL_species-specific_OCRs_rank.bed"

SELF_F = "/dors/capra_lab/data/ucsc/hg38/self/hg38_repeats_self_coor.bed"


#%% FUNCTIONS
def get_vars(branch, msa_way, FDR, prune):

    if prune == False:
        CONACC_F =f"/dors/capra_lab/users/fongsl/tyler/data/CON_ACC/all/multiz30way_{branch}/all_con_acc.bed"
        PATH = "/".join(CONACC_F.split("/")[:-1]) + "/" # the path
        if FDR == True:
            RE = f"/dors/capra_lab/projects/enhancer_ages/tyler/results/CON_ACC/{msa_way}way_{branch}_FDR/"
        else:
            RE = f"/dors/capra_lab/projects/enhancer_ages/tyler/results/CON_ACC/{msa_way}way_{branch}/"

        if os.path.exists(RE) == False:
            os.mkdir(RE)
    elif prune == True:
        CONACC_F =f"/dors/capra_lab/users/fongsl/tyler/data/CON_ACC/prune/multiz30way_{branch}/all_con_acc.bed"
        PATH = "/".join(CONACC_F.split("/")[:-1]) + "/" # the path
        if FDR == True:
            RE = f"/dors/capra_lab/projects/enhancer_ages/tyler/results/CON_ACC/prune/{msa_way}way_{branch}_FDR/"
        else:
            RE = f"/dors/capra_lab/projects/enhancer_ages/tyler/results/CON_ACC/prune/{msa_way}way_{branch}/"

        if os.path.exists(RE) == False:
            os.mkdir(RE)

    return CONACC_F, PATH, RE

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
        #print(cmd)
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

def combine_annotations(df):

    df["sp_te_slf"] = df["species"] + "-" + df["te_bin"].map(str) + "-" +  df["slf_bin"].map(str)
    df["te_slf"] = df["te_bin"].map(str) +"-" +  df["slf_bin"].map(str)
    df["sp_slf"] = df["species"] + "-" +  df["slf_bin"].map(str)
    return df

def merge_liftover():


    multiple_maps = {
    "Hg38-rheMac8": "/dors/capra_lab/users/fongsl/tyler/data/liftover/hg38_to_RheMac8_multimap.txt",
    "RheMac8-hg38": "/dors/capra_lab/users/fongsl/tyler/data/liftover/rheMac8_to_Hg38_multimap.txt"
    }

    results = {}
    for comp, file in multiple_maps.items():

        # open the file
        df = pd.read_csv(file, header = None)

        # format the annotation column
        count_col = f"{comp}_map_count"
        df = df.rename(columns = {0:"annot"})

        # get multimap counts
        df[count_col] = df["annot"].apply(lambda x: x.split(" Peak_")[0]).map(int)

        # get peak id column
        df["id"] = df["annot"].apply(lambda x: "Peak_" + x.split("Peak_")[1]).map(str)
        df = df.drop(["annot"], axis = 1) # drop annot column
        results[comp] = df

    merged = pd.merge(results["Hg38-rheMac8"], results["RheMac8-hg38"], on = "id")

    return merged

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

    # combine annotations
    df = combine_annotations(df)

    # format the ID column
    df.id = df.id.apply(lambda x: x.split('"')[1])

    # merge liftover multi-map counts
    merged = merge_liftover()
    df2 = pd.merge(df, merged)

    return  df

def fdr_correction(df):

    # reverse the -log10p calculation to get actual p values
    # take the absolute value of the log10p conacc (the sign only indicates if element is accelerated or conserved)
    #df["conacc_abs"] = abs(df["conacc"])
    df["p_conacc"] = 10**(-1*(abs(df["CON_ACC"])))
    pvals = df["p_conacc"]

    df["reject_null"], df["FDR_P"] = statsmodels.stats.multitest.fdrcorrection(pvals, alpha=0.05)

    return df

def plot_cdf(x, data, hue, title, pal):

    sns.set("talk")

    n = data.groupby([hue])["id"].count().reset_index()
    nzoom_acc = data.loc[data[x]<0].groupby([hue])["id"].count().reset_index()
    nzoom_con = data.loc[data[x]>0].groupby([hue])["id"].count().reset_index()

    fig, axes =plt.subplots(1,3, figsize = (18,6))

    ax = axes[0]
    sns.ecdfplot(data,
    x=x,
    hue = hue,
    palette = pal,
    ax = ax)

    ax.set( title = title + f"\n{BRANCH} branch",
    xlim = (data[x].min(),data[x].max()),
    #ylim = (0.05, 0.70),
    xlabel = f"CON_ACC\n{n}"
    )


    ax = axes[1]
    sns.ecdfplot(data.loc[data[x]<0],
    x=x,
    hue = hue,
    legend=False,
    palette = pal,
    ax = ax)

    ax.set(title = title + f"\n{BRANCH} branch, CON_ACC<0",
    xlim = (data[x].min(),0),
    #ylim = (0.05, 0.70),
    xlabel = f"CON_ACC\n{nzoom_acc}"
    )


    ax = axes[2]
    sns.ecdfplot(data.loc[data[x]>0],
    x=x,
    hue = hue,
    legend=False,
    palette = pal,
    ax = ax)

    ax.set(title = title + f"\n{BRANCH} branch, CON_ACC>0",
    xlim = (0, data[x].max()),
    #ylim = (0.05, 0.70),
    xlabel = f"CON_ACC\n{nzoom_con}"
    )

    #print("min and max",data[x].min(),data[x].max())

    outf = f"{RE}{title}_{MSA_WAY}way_{BRANCH}.pdf"
    plt.savefig(outf, bbox_inches = "tight")

def get_medians(hue, x, data, title):
    meds = data.groupby(hue)[x].median().reset_index()

    melted = pd.melt(meds, id_vars=[x], value_vars=[hue])
    melted["title"] = title

    return melted

#%% run functions

for BRANCH in branches:

    CONACC_F, PATH, RE = get_vars(BRANCH, MSA_WAY, FDR, PRUNE)

    df = format_f(CONACC_F) # format the file

    outTE = intersect_te(CONACC_F, PATH, TE_F) # intersect w/ repeatmasker TE

    outSP = intersect_species_specific(CONACC_F, PATH, SPECIES_DA_F) # intersect species specific accessibility calls.

    outSLF = intersect_self(CONACC_F, PATH, SELF_F) # intersect self accessibility calls.

    # the amalgamated df
    newdf_ = compile_df(outTE, outSP, outSLF)
    newdf_ = fdr_correction(newdf_)


    if FDR == True:
        newdf = newdf_.loc[newdf_.reject_null==True]

    elif FDR == False:

        newdf = newdf_.loc[newdf_.p_conacc < 0.05]
    else:
        newdf = newdf_

    # dictionary of analyses to do
    # title :[data, hue, pal]
    analysis = {
    "TE_stratified": [newdf.sort_values(by = "sp_te"), "te_bin", "Set1"],
    "all" : [newdf.sort_values(by = "sp_te"), "sp_te", "tab20"],
    "TE_excluded":[newdf.loc[newdf.te_bin ==0], "sp_te", "Set1"],
    "TE_included":[newdf.loc[newdf.te_bin ==1], "sp_te", "Set1"],
    "Species_only":[newdf, "species", "Set1"],
    "Self":[newdf, "slf_bin", "tab10"],
    "Sp_Te_Self":[newdf.sort_values(by = "sp_te_slf"), "sp_te_slf", "tab20c"],
    "Te_Self":[newdf.sort_values(by = "te_slf"), "te_slf", "tab10"],
    "Self_excluded":[newdf.loc[newdf.slf_bin ==0], "sp_slf", "tab10"],
    "Self_included":[newdf.loc[newdf.slf_bin ==1], "sp_slf", "tab10"],
    }

    #
    median_results = pd.DataFrame()

    x = "CON_ACC"

    for title, vals in analysis.items():

        data, hue, pal = vals[0], vals[1], vals[2]
        plot_cdf(x, data, hue, title, pal)
        melted = get_medians(hue, x, data, title)
        median_results = median_results.append(melted)

        # save the medians file
        median_results["branch"]  = BRANCH
        median_results['msa_way'] = MSA_WAY
        median_out = f'{RE}medians_{BRANCH}_{MSA_WAY}.tsv'
        median_results.to_csv(median_out, sep = '\t', index = False)




#%%
