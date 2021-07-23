import os, sys
import subprocess

# add the path with the file containing all the paths to other files.
PATH = '/dors/capra_lab/users/fongsl/enh_ages/tyler/evolutionary_conservation'
sys.path.append(PATH)

import config # import the config.py file with all the files for tyler's project
import pandas as pd

files = [config.all, config.hu_specific, config.hu_specific_noTE, config.rhe_specific, config.rhe_specific_noTE]
#%%
for file in files:
    df = pd.read_csv(file, sep = '\t', header = None)
    newdf = df.loc[~df[0].str.contains("KI")]
    print(newdf.shape, df.shape)
    newdf.to_csv(file, sep = '\t', header = False, index = False)
    #break


df.shape
