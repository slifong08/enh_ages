import pandas as pd
import os

PATH = "/dors/capra_lab/data/fantom/fantom5/fantom5_phase1-2_enhancers"
FILE = "human_permissive_enhancers_phase_1_and_2_expression_count_matrix.txt.gz"

FULL = os.path.join(PATH, FILE)

OUTPATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/"
OUTFILE = "HEPG2_FANTOM5_hg19.bed"
FULL_OUT = os.path.join(OUTPATH, OUTFILE)


HEPG2 = ["Id", "CNhs12328", "CNhs12329", "CNhs12330"]


#%% load dataframe

df = pd.read_csv(FULL, sep = '\t', usecols = HEPG2)

df = df.loc[~df["Id"].str.contains("chrX")] # remove chrX


#%% make the bed file


df["chr"] = df["Id"].apply(lambda x: x.split(":")[0])
df["start"] = df["Id"].apply(lambda x: (x.split(":")[1]).split("-")[0])
df["end"] = df["Id"].apply(lambda x: (x.split(":")[1]).split("-")[1])


#%% get df with more than 5 reads in each replicate, n = 728


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
