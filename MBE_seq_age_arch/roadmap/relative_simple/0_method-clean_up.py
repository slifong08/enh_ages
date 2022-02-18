import glob
import pandas as pd

import os, sys
from scipy import stats


#%% Files


#to_do_list = ["E114"]
to_do_list = []

pre_path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/"
path = "%sbreaks/" % pre_path
path
# file contains both enhancer + 10 shuffled breaks.
enhFs= glob.glob("%sROADMAP_E*_enh_and_shuf_age_arch_summary_matrix.bed" % (path))

len(enhFs)
#%%
enhFs[1]
#%%
for enh in enhFs:

    sid = (enh.split("/")[-1]).split("_")[1]


    print(sid)

    if sid == "E063":
        df = pd.read_csv(enh, sep = '\t', header = None,low_memory=False,
error_bad_lines = False)

        print(len(df))
        df.to_csv(enh, sep = '\t', header =None, index = None)

        del df
#%%
print("hello")
