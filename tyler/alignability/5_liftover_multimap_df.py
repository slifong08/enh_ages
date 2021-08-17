import os, sys
import subprocess
import pandas as pd

def get_liftover_mappings():

    multiple_maps = {
    "Hg38-rheMac8": "/dors/capra_lab/users/fongsl/tyler/data/liftover/hg38_to_RheMac8_multimap.txt",
    "RheMac8-hg38": "/dors/capra_lab/users/fongsl/tyler/data/liftover/rheMac8_to_Hg38_multimap.txt"
    }

    results = {}
    for comp, file in multiple_maps.items():

        # open the file
        df = pd.read_csv(file, header = None)

        # format the annotation column
        count_col = f"{comp}_map_count"
        df = df.rename(columns = {0:"annot"})

        # get multimap counts
        df[count_col] = df["annot"].apply(lambda x: x.split(" Peak_")[0]).map(int)

        # get peak id column
        df["id"] = df["annot"].apply(lambda x: "Peak_" + x.split("Peak_")[1]).map(str)
        df = df.drop(["annot"], axis = 1) # drop annot column
        results[comp] = df

    merged = pd.merge(results["Hg38-rheMac8"], results["RheMac8-hg38"], on = "id")

    return merged
m = get_liftover_mappings()
m.head()
m.groupby(['RheMac8-hg38_map_count', "Hg38-rheMac8_map_count"])["id"].count()
m['Hg38-rheMac8_map_count'].value_counts()
m['Hg38-rheMac8_map_count'].value_counts(normalize = True)
m['RheMac8-hg38_map_count'].value_counts(normalize = True)
