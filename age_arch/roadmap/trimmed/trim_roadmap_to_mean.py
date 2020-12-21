# intro

# created 2019-08-05
# sarahfong

# Goal - remove sex chromosomes, exons, enhancers longer than 10kb
# trim roadmap enhancers to mean lengths



import glob
import os, sys
import pandas as pd
import datetime



trim_len = "mean"
trim_len = 310

#  Trim function

def trim(bedfile, trim_len):

    df = pd.read_csv(bedfile, sep ='\t', header = None, usecols = [0,1,2]) # open file
    df.columns = ["chr", "start", "end"] # name columns

    df["id"] = df.chr + ":" + df.start.map(str) + "-" + df.end.map(str)
    df["old_len"]= df.end - df.start # calculate enhancer length

    # clean up
    df = df.loc[df.old_len<10000]
    df = df.loc[df.chr != "chrX"]

    if trim_len == "mean": # take the mean

        trim_len =int(df["old_len"].mean()) # mean enhancer length in dataset
        print(trim_len)

    df["midpoint"] = (df.start + (trim_len)/2).astype(int) # identify the midpoint of each enhancer

    df["new_len"] = trim_len

    df["start_new"] =((df.midpoint - (trim_len/2)).round(0)).astype(int) # calculate new start as the midpoint - (mean length/2)

    df["end_new"] = ((df.midpoint + (trim_len/2)).round(0)).astype(int)

    trimmed = df[["chr", "start_new", "end_new", "id", "old_len", "new_len"]].drop_duplicates()

    return trimmed

path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19/download/h3k27ac_plus_h3k4me3_minus_peaks/"

# -a
test_list = glob.glob("%sHsap_H3K27ac_plus_H3K4me3_minus_E*/cleaned_Hsap_H3K27ac_plus_H3K4me3_minus_E*/shuffle/breaks/*summary_matrix.bed" % path)


# make a output directory


for infile in test_list:
    fpath = "/".join(infile.split("/")[:-1])
    outpath = "%strimmed/" % fpath

    if os.path.exists(outpath) == False:
        os.mkdir(outpath)

    sid = (infile.split("/")[-4]).split("_")[-1]
    print(sid)
    # Entire Enhancer #
    # Trim enhancers to mean lengths #


    trimmed_df = trim(infile, trim_len)
    trimmed_df.to_csv("%strimmed-%s_%s.bed" % (outpath,  trim_len, sid), sep = '\t', header = False, index = False)
