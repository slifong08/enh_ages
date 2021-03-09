import pandas as pd
import os, sys

#bedfile = sys.argv[1] #infile
BEDPATH = "/dors/capra_lab/data/encode/encode3_hg38/TF/"
BEDFILE = "trimmed_encRegTfbsClusteredWithCells.hg38.bed"
bedfile = os.path.join(BEDPATH, BEDFILE)


bedfile
#%% function for splitting file on cell line


def split_on_cell_line(cell_line, df, outfile):

    cl_df = df.loc[df.cell_line.str.contains(cell_line)]

    cl_df.to_csv(outfile, sep = '\t', header = False, index = False)
    print("made this:", outfile)


#%% read the file

path = "/".join(bedfile.split("/")[:-1])

df = pd.read_csv(bedfile, sep = '\t', header = None, usecols = [0,1,2,3, 4, 6, 8])

cols = ["chr", "start_new", "end_new", "old_enh_id",
"old_len",  "tf", "cell_line"]

df.columns = cols # rename the columns
df["tf"] = df['old_enh_id'].apply(lambda x: x.split("_")[0])


#%% process per cell line. There should be 129 cell lines total.
df.head()

cls = df.cell_line.unique() # array of unique cell lines, and overlaps

for cell_line in cls:

    if "/" in cell_line:
        c = "_".join(cell_line.split("/")) # do some formatting.
    else:
        c = cell_line

    if "," not in c: # split only on single cell lines

        outfile = "%s/cells/%s.bed" % (path, c)
        split_on_cell_line(cell_line, df, outfile)
