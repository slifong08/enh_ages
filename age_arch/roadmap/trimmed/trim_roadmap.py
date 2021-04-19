# intro

# created 2019-08-05
# sarahfong

# Goal - remove sex chromosomes, exons, enhancers longer than 10kb
# trim roadmap enhancers to mean lengths

import glob
import numpy as np
import os, sys
import pandas as pd
import datetime


#  Trim function


def make_chr_list():
    n = list(np.arange(1, 23))
    #n.append("X")

    chr_list = []
    for num in n:
        chrn = "chr" + str(num)
        chr_list.append(chrn)

    return chr_list

def trim(bedfile, trim_len, sid):

    cols = ["#chr", "start", "end"] # name columns
    df = pd.read_csv(bedfile,
        sep ='\t',
        header = None,
        usecols = [0,1,2],
        names = cols ) # open file


    df["old_id"] = df["#chr"] + ":" + df.start.map(str) + "-" + df.end.map(str)
    df["old_len"]= df.end - df.start # calculate enhancer length

    # clean up
    df = df.loc[df.old_len<10000]
    chr_list = make_chr_list()
    df = df.loc[df["#chr"].isin(chr_list)] # autosomes only

    trim_sid = str(trim_len)
    if trim_len == "mean": # take the mean

        trim_len = int(df["old_len"].mean()) # mean enhancer length in dataset

    df["midpoint"] = (df.start + (df.old_len)/2).astype(int) # identify the midpoint of each enhancer

    df["new_len"] = trim_len

    df["start_new"] = ((df.midpoint- (trim_len/2)).round(0)).astype(int) # calculate new start as the midpoint - (mean length/2)

    df["end_new"] = ((df.midpoint + (trim_len/2)).round(0)).astype(int)

    df["new_id"] = df["#chr"] + ":" + df.start_new.map(str) + "-" + df.end_new.map(str)

    trimmed = df[["#chr", "start_new", "end_new",
    "old_id", "old_len",
    "new_id", "new_len"]].drop_duplicates()

    # write the file
    outpath = "/".join(bedfile.split("/")[:-1]) + "/trimmed/"
    outfile = f"{outpath}trim_{trim_sid}_{sid}.bed"

    if os.path.exists(outpath) == False: # make a path
        os.mkdir(outpath)

    trimmed.to_csv(outfile, sep = '\t', index = None)

    print("trimmed and made this file:", outfile)

    return trimmed

#%% Trim them bed files!

PATH = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19/download/h3k27ac_plus_h3k4me3_minus_peaks/"

test_list = glob.glob(f"{PATH}Hsap_H3K27ac_plus_H3K4me3_minus_E*/non-genic/no-exon_E*.bed")

trim_lens =["mean", 310]

for infile in test_list:

    sid = (infile.split("no-exon_")[1]).split(".")[0]

    if "summary_matrix" not in sid:

        print(sid)

        for trim_len in trim_lens:

            trim(infile, trim_len, sid)


#%% trim the all roadmap enh bedfile
f = f"{PATH}all_roadmap_enh/no-exon_all_roadmap_enh.bed"
sid = "no-exon_all_roadmap_enh"

for trim_len in trim_lens:

    trim(f, trim_len, sid)
