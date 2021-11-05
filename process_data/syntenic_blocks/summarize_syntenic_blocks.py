
import glob
import numpy as np
import os
import pandas as pd
import subprocess
import sys

BUILD = "hg38"
PATH = f"/dors/capra_lab/data/ucsc/{BUILD}/synteny_age_{BUILD}"
CHR = sys.argv[1]

# drop non-important columns and collapse coordinates by MRCA. 

def drop(f, chr_):
    os.chdir(os.path.join(PATH, "summarized"))
    cutf = f"cut_{chr_}_syn_age.bed"
    cmd = f"cut -f 1,2,3,8 {f} > {cutf}"
    subprocess.call(cmd, shell = True)
    
    return cutf

def lump(df):

    # parse through rows, consolidate coordinates of the same age. 
    start, stop = 0, 0
    last_mrca = -1

    results = {}
    
    for i, row in df.iterrows():

        chr_, start_, stop_ = row.iloc[0], row.iloc[1], row.iloc[2]
        mrca = round(row.iloc[3], 3)

        if last_mrca == -1:  # for the first row of the segment
            start = start_
            stop = stop_
            last_mrca = mrca

        elif mrca == last_mrca:  # still in the same MRCA
            stop = stop_  # add to the stop

        elif mrca!= last_mrca:
            new_row = [chr_, start, stop, last_mrca]
            results[i] = new_row

            # reset values
            start = start_
            stop = stop_
            last_mrca = mrca

    return results

# Stack rows, make dataframe. 

def stack_n_save(results, RE):
    
    data = np.array(["#chr", "start", "stop", "mrca"])
    for val in results.values():
        a = np.array(val)
        data = np.vstack((data,a,))


    re = pd.DataFrame(data)
    re.to_csv(RE, sep = '\t', header = None, index = False)



Fzipped = os.path.join(PATH, f"{CHR}_syn_age.bed.gz")

F = os.path.join(PATH, f"{CHR}_syn_age.bed")

RE = os.path.join(PATH, "summarized",  f"{CHR}_syn_age.bed")

if os.path.exists(RE) is False:
    print(RE)
    # unzip
    cmd = f"gunzip {Fzipped}"
    if os.path.exists(Fzipped) is True:
        subprocess.call(cmd, shell = True)

    cut_f = drop(F, CHR)

    df = pd.read_csv(cut_f, sep = '\t', header = None)

    results = lump(df)

    stack_n_save(results, RE)

    # rezip
    cmd = f"gzip {F}"
    subprocess.call(cmd, shell = True)

# zip all the resulting files

os.chdir(os.path.join(PATH, "summarized"))
cmd = "gzip *.bed"
#subprocess.call(cmd, shell = True)