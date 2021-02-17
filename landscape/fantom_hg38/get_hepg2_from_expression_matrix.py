import pandas as pd
import os

PATH = "/dors/capra_lab/data/fantom/fantom5/hg38/"
FILE = "F5.hg38.enhancers.expression.matrix"

FULL = os.path.join(PATH, FILE)

OUTPATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/hg38"
OUTFILE = "HEPG2_FANTOM5_hg38.bed"
FULL_OUT = os.path.join(OUTPATH, OUTFILE)


HEPG2 = ["Unnamed: 0", "CNhs12328", "CNhs12329", "CNhs12330"]


#%% load dataframe


df = pd.read_csv(FULL, sep = '\t', usecols = HEPG2)

df = df.loc[~df["Unnamed: 0"].str.contains("chrX")] # remove chrX


#%% make the bed file


df["chr"] = df["Unnamed: 0"].apply(lambda x: x.split(":")[0])
df["start"] = df["Unnamed: 0"].apply(lambda x: (x.split(":")[1]).split("-")[0])
df["end"] = df["Unnamed: 0"].apply(lambda x: (x.split(":")[1]).split("-")[1])


#%% get df with more than 5 reads in each replicate, n = 674


morethan5_df = df[
(df["CNhs12328"]>5)
& (df["CNhs12329"]>5)
& (df["CNhs12330"]>5)
]

# rearrange the columns in bed format
morethan5_df = morethan5_df[[
"chr", "start", "end",
#"CNhs12328", "CNhs12329", "CNhs12330"
]]


morethan5_df.info()

#%% write the file

morethan5_df.to_csv(FULL_OUT, sep = '\t', header = False, index = False)
