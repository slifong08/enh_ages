import glob
import pandas as pd
import os, sys

#bedfile = sys.argv[1] # must have 6 fields: chr, start, end, tf, reads, cell_line
bedfiles = glob.glob("/dors/capra_lab/users/fongsl/tyler/data/GG-LL_*/*.bed")
trim_len = 466 #291 # trim to 291 bp.
bedfiles

def trim(bedfile, trim_len):

    outpath = "/".join(bedfile.split("/")[:-1])
    outfile = outpath + f"/trimmed_{trim_len}_" + bedfile.split("/")[-1]

    df = pd.read_csv(bedfile, sep ='\t', header = None, usecols = [0,1,2,3]) # open file

    df.columns = ["#chr", "start", "end", "id",] # name columns

    df["old_enh_id"] = df.id + "_" + df["#chr"] + ":" + df.start.map(str) + "-"+ df.end.map(str)

    df["old_len"]= df.end - df.start # calculate enhancer length

    df["new_len"] = trim_len

    # BUG WAS HERE. GOT OLD_LEN AND TRIM_LEN VARS mixed.
    df["midpoint"] = (df.start + (df.old_len)/2).astype(int) # identify the midpoint of each enhancer

    df["start_new"] = ((df.midpoint- (trim_len/2)).round(0)).astype(int) # calculate new start as the midpoint - (mean length/2)

    df["end_new"] = ((df.midpoint + (trim_len/2)).round(0)).astype(int)

    trimmed = df[["#chr", "start_new", "end_new", "old_enh_id",
    "old_len", "new_len", "id"]].drop_duplicates()


    trimmed.to_csv(outfile, sep = '\t', index = None)

    print("trimmed and made this file:", outfile)


#%% Trim that bed file!

for bedfile in bedfiles:
    if "trimmed" not in bedfile:
        trim(bedfile, trim_len)
