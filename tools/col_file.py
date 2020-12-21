import pandas as pd
import numpy as np
import os
#%%
path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/breaks/"
f = "%sROADMAP_E038_enh_and_shuf_age_arch_summary_matrix.tsv"% path

df = pd.read_csv(f, sep = '\t')
df.head()
#%%
headerf = "%ssummary_matrix_header.txt" % (path)


cols = df.columns.to_list() # make the columns into a list
np.savetxt(fname = headerf, X=cols,  delimiter=",",  fmt='%s', header = "ROADMAP_E*_enh_and_shuf_age_arch_summary_matrix.tsv") # save the column headers
#%%
pd.read_csv(headerf)
