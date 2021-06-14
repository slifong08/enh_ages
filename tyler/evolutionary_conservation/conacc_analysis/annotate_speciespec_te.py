import argparse
import numpy as np
import os, sys
import pandas as pd
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
F ="/dors/capra_lab/users/fongsl/tyler/data/CON_ACC/all/multiz30way_hg38-rheMac8/all_con_acc.bed"
BRANCH = "hg38-rheMac8"
MSA_WAY = 30
PATH = "/".join(F.split("/")[:-1]) + "/" # the path


#%% FUNCTIONS

def format_f(F):

    cols = ["#chr", "b", "bf", "start", "end", "conacc", "?", "??", "id"]
    test = pd.read_csv(F, sep  = '\t', nrows = 5)

    if "start" not in list(test): # do you need to rename the columns?
        df = pd.read_csv(F, sep  = '\t', header = None, names = cols)
        df = df.drop(["b", "bf", "?", "??"], axis = 1) # drop the columns you don't need
        df.to_csv(F, sep = '\t', index = False) # save this formatted file.
    else: # or you've renamed them already and just need to read the df.
        df = pd.read_csv(F, sep  = '\t')

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


def intersect_te(F, PATH):
    TE = "/dors/capra_lab/data/transposable_elements/repeatmasker/hg38.bed"
    #te = pd.read_csv(TE, sep = '\t', skiprows = 1)

    # get info to write new files
    sample_id = (F.split("/")[-1]).split(".bed")[0]
    outTE = f"{PATH}{sample_id}_TE.bed"

    if os.path.exists(outTE) == False: # if you haven't done this intersection already...

        cmd = f"bedtools intersect -a {F} -b {TE} -wao > {outTE}"
        subprocess.call(cmd, shell = True)

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


def intersect_species_specific(outTE, PATH):

    SPECIES_F = "/dors/capra_lab/users/fongsl/tyler/data/GG-LL_species-specific_OCRs_rank/GG-LL_species-specific_OCRs_rank.bed"

    sample_id = (outTE.split("/")[-1]).split(".bed")[0]
    outTE_SP = f"{PATH}{sample_id}_species.bed"

    if os.path.exists(outTE_SP) == False:
        cmd = f"bedtools intersect -a {outTE} -b {SPECIES_F} -wao > {outTE_SP}"
        subprocess.call(cmd, shell = True)

    # format the output dataframe
    # col names
    cols = ['#chr', 'start', 'end', 'CON_ACC', 'id',
 'te', 'te_fam', 'len_te_overlap', 'te_id',
 "chr_sp", "start_sp", "end_sp", "", "species", "spoverlap"]

    # add id for te locus
    idname = "sp_id"

    #specify locus
    chr, start, end = "chr_sp", "start_sp", "end_sp"

    # add column names, id, drop chr, start, end, and save the new formatted file.
    test_colnames_add_id(outTE_SP, cols, chr, start, end, idname)

    return outTE_SP


#%% run functions

df = format_f(F) # format the file

outTE = intersect_te(F, PATH) # layer TE info
list(pd.read_csv(outTE, sep = '\t'))
outTE
outTE_SP = intersect_species_specific(outTE, PATH) # layer species specific info
#%%
pd.read_csv(outTE_SP, sep = '\t')
