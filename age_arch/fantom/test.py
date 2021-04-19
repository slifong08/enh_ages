import pandas as pd

PATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/non-genic/"
F = f"{PATH}SHUFFLE_NOEXON_FANTOM_enh_age_arch_summary_matrix.tsv"
F = '/dors/capra_lab/projects/enhancer_ages/fantom/data/non-genic/no-exon_shuffle_fantom_age_arch_summary_matrix.bed'
df = pd.read_csv(F, sep = ',')
list(df)
df = df.drop(['Unnamed: 0'], axis = 1)
df.head()
df.to_csv(F, sep = '\t', index = False),
df = df.rename(columns={'datatype': 'sample_id', 'chr_enh': '#chr_enh'})

cols = ['#chr_enh',
 'start_enh',
 'end_enh',
 'enh_id',
 'sample_id',
 'seg_index',
 'core_remodeling',
 'arch',
 'mrca',
 'taxon',
 'mrca_2',
 'taxon2','enh_len']
df =df[cols]
df = df.loc[df.enh_len>5]
df.to_csv(F)

#%%
F
