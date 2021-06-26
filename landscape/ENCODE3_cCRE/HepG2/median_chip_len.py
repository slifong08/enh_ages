import os, sys
import numpy as np
import pandas as pd

PATH = "/dors/capra_lab/data/encode/encode3_hg38/TF/"
F = f"{PATH}encRegTfbsClusteredWithCells.hg38.bed"

df = pd.read_csv(F, sep = '\t', header = None, usecols = [3,4])

np.median(df[4])


PATH_ENCODE = "/dors/capra_lab/projects/enhancer_ages/encode/data/"
F_ENCODE = f"{PATH_ENCODE}ELS_combined_HepG2.bed"

df_ = pd.read_csv(F_ENCODE, sep = '\t', header = None, usecols = [1,2])
df_[3] = df_[2]-df_[1]
np.median(df_[3])
