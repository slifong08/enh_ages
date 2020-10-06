import os
import pandas as pd

path = "/dors/capra_lab/data/enhancers/emera16/"
f = "%semera_2016_neocortical_dev_enhancers_hu_ms.csv" % (path)
print(f)


#%% open dataframe


df = pd.read_csv(f)
df.head()


#%% make bedfile


df["chr"] = df["Human enhancer (shared with mouse)"].apply(lambda x: x.split(":")[0])
df["start"] = df["Human enhancer (shared with mouse)"].apply(lambda x: (x.split(":")[1]).split("-")[0])
df["end"] = df["Human enhancer (shared with mouse)"].apply(lambda x: x.split("-")[1])

#%%
df.head()
#%%


# rename columns
df.columns = ['Human_enhancer_shared_w_mouse', 'Stage_in_human',\
'Phylogenetic_age_assignment', 'chr', 'start', 'end']

# reorganize the dataframe
df = df[['chr', 'start', 'end', 'Human_enhancer_shared_w_mouse', \
"Stage_in_human", "Phylogenetic_age_assignment"]]

outf = "%semera_2016_neocortical_dev_enhancers_hu_ms.bed" % path
df.to_csv(outf, sep = '\t', header = False, index = False)
