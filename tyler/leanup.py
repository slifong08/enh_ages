import os, sys
import pandas as pd

PATH = "/dors/capra_lab/users/fongsl/tyler/data/GG-LL_species-specific_OCRs_rank/shuffle/ages/"
FILE = "shuf-GG-LL_species-specific_OCRs_rank-1_ages.bed"
F = os.path.join(PATH, FILE)

age = pd.read_csv(F, sep = '\t', header = None)
#%%
age.head()

#%%
age[0].unique()
#%%
age = age.loc[~age[0].str.contains("random")]

age.to_csv(F, sep ="\t", header = None, index = False)
