import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

#%%

path = "/dors/capra_lab/users/fongsl/data/ensembl/breaks/"
f = "%sensGene_hg19_coding_exons-autosome_only_enh_ages_enh_age_arch_summary_matrix.bed" %path

df= pd.read_csv(f, sep ='\t', header = None)

cols = ["chr", "start", "end", 'enh_id', "end_id2", "shuf_id", "core_remodeling",
"arch", "seg_index", "mrca", "len", "taxon", "mrca_2", "taxon2", 'mya', "mya2", "seg_den", "dataset"]
df.columns = cols

df.head()
#%%
sns.countplot(df.arch)
#%%
df.mrca_2 = df.mrca_2.round(3)
sns.countplot(df.mrca_2)
