import pandas as pd
import glob

path = "/dors/capra_lab/data/fantom/fantom5/facet_expressed_enhancers/tissues/"
outpath = "/dors/capra_lab/users/fongsl/data/ensembl/"

f = "/dors/capra_lab/users/fongsl/data/ensembl/hg19_blacklist_gap_ensemblexon.bed"
fid = (f.split("/")[-1])


df = pd.read_csv(f, sep = '\t', header = None, usecols = [0,1,2,3,])
df.sort_values(by = [0,1,2,3,])
df = df.fillna(".")
df.head()


#%%
outf = "%sclean_%s" %(outpath, fid)
df.to_csv(outf, sep = '\t', header = False, index = False)
#%%
len(df)

#%%
outf
