import glob
import pandas as pd

fs = glob.glob('/dors/capra_lab/projects/enhancer_ages/fantom/data/architecture_coordinates/all_fantom/*.bed')
fs[0]
#%%
for f in fs:
    print(f)
    df = pd.read_csv(f, sep = '\t', header =None)
    df = df.loc[df[0]!= "chrX"]

    df.to_csv(f, sep ='\t', header = False, index = False)
#%%
SHUF_PATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/non-genic/"
shufFS = glob.glob("%sSHUFFLE_FANTOM_no-exon_shuf-all_fantom_enh_112_tissues-*_age_breaks_summary_matrix.bed" % SHUF_PATH)
for f in shufFS:
    print(f)
    df = pd.read_csv(f, sep = '\t', header =None)
    df = df.loc[df[0]!= "chrX"]

    df.to_csv(f, sep ='\t', header = False, index = False)
    
