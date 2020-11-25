import glob
import pandas as pd
fs = glob.glob("/dors/capra_lab/projects/enhancer_ages/fantom/data/shuffle/breaks/oldest_ages/*.bed")


#%%

freq_list = []

for f in fs:
    df = pd.read_csv(f, sep = '\t', usecols =[7])
    df.columns = ["cr"]
    complex = len(df.loc[df.cr==1])

    simple = len(df.loc[df.cr!=1])
    freq= simple/(simple + complex)
    freq_list.append(freq)

#%%
import numpy as np
np.mean(freq_list)
