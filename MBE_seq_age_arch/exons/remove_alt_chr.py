import os, sys
import pandas as pd

path = "/dors/capra_lab/users/fongsl/data/ensembl/"
f = "%sensGene_hg19_coding_exons.bed" % path

df = pd.read_csv(f, sep = '\t', header = None)
df.head()
#%%
keep_chr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
       'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
       'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
       'chr21', 'chr22']
print(len(df))
#%%
df = df.loc[df[0].isin(keep_chr)] # keep only the chromosomes in list
print(len(df))
df = df[[0,1,2]].drop_duplicates() # drop exons logged twice for different transcript models
print(len(df))
#%%
753667-699667
# lose 54000 non-autosome exons

#260065 unique exon coordinates
#%%
df[0].unique()

outf = "%sensGene_hg19_coding_exons-autosome_only.bed" %path

df.to_csv(outf, sep = '\t', header = None, index = None)
