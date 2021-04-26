import pandas as pd
import os, sys

bedfile = sys.argv[1] #infile
BEDPATH = "/dors/capra_lab/data/encode/encode3_hg38/TF/liftOver_hg19/"
BEDFILE = "trimmed_encRegTfbsClusteredWithCells.liftOver.to.hg19.bed"
bedfile = os.path.join(BEDPATH, BEDFILE)
BEDOUT = "trimmed_encRegTfbsClusteredWithCells.liftOver.to.hg19_rename.bed"
bedout = os.path.join(BEDPATH, BEDOUT)

bedfile
#%% function for splitting file on cell line


def split_on_cell_line(cell_line, df, outfile):

    cl_df = df.loc[df.cell_line.str.contains(cell_line)]

    cl_df.to_csv(outfile, sep = '\t', header = False, index = False)
    print("made this:", outfile)


#%% read the file

path = "/".join(bedfile.split("/")[:-1])

df = pd.read_csv(bedfile, sep = '\t', header = None)

cols = ["#chr", "start_new", "end_new", "old_enh_id",
"old_len",  "tf", "cell_line"]

df.columns = cols # rename the columns
df["tf"] = df['old_enh_id'].apply(lambda x: x.split("_")[0])

df.to_csv(bedout, sep = '\t', index = False)
