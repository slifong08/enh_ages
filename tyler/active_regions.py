import os, sys
import pandas as pd
from scipy import stats

PATH = "/data/hodges_lab/ATAC-STARR_V2/data/ATAC-STARR_cts_matricies/"
FILE = "GM12878inGM12878_counts.tsv"

F = os.path.join(PATH, FILE)
F
#%%
cols = ["peakid", "chr", "start",
"end","strand", "len",
"dna_rep1", "dna_rep2", "rna_rep1", "rna_rep2"]
df = pd.read_csv(F, sep = '\t', skiprows =1, )
df.columns = cols

# 1.14e7 rows

#%%
df.describe()
#%%
df.head(10)
