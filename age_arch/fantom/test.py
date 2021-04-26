import pandas as pd

PATH = "/dors/capra_lab/projects/enhancer_ages/"
F = f"{PATH}hg19_syn_gen_bkgd.tsv"

df = pd.read_csv(F, sep = '\t')
list(df)

df.head()
#%%
df[["taxon2", "base_count_2"]].drop_duplicates()
#%%
df[["taxon2", "block_count_2"]].drop_duplicates()
