import pandas as pd
import os, sys

bedfile = sys.argv[1] # must have 6 fields: chr, start, end, tf, reads, cell_line
trim_len = 30 # trim to 30 bp.


def trim(bedfile, trim_len):

    outpath = "/".join(bedfile.split("/")[:-1])
    outfile = outpath + "/trimmed_" + bedfile.split("/")[-1]

    df = pd.read_csv(bedfile, sep ='\t', header = None) # open file

    df.columns = ["chr", "start", "end", "tf", "reads", "cell_line"] # name columns

    df["old_enh_id"] = df.tf + "_" + df.chr + ":" + df.start.map(str) + "-"+ df.end.map(str)

    df["old_len"]= df.end - df.start # calculate enhancer length

    df["new_len"] = trim_len

    # BUG WAS HERE. GOT OLD_LEN AND TRIM_LEN VARS mixed.
    df["midpoint"] = (df.start + (df.old_len)/2).astype(int) # identify the midpoint of each enhancer

    df["start_new"] = ((df.midpoint- (trim_len/2)).round(0)).astype(int) # calculate new start as the midpoint - (mean length/2)

    df["end_new"] = ((df.midpoint + (trim_len/2)).round(0)).astype(int)

    trimmed = df[["chr", "start_new", "end_new", "old_enh_id",
    "old_len", "new_len", "tf", "reads", "cell_line"]].drop_duplicates()


    trimmed.to_csv(outfile, sep = '\t', header = None, index = None)

    print("trimmed and made this file:", outfile)


# Trim that bed file!


trim(bedfile, trim_len)
