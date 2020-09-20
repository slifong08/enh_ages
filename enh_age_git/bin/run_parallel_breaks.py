import argparse
import os, sys
import glob
import subprocess
import pandas as pd
import numpy as np
from itertools import groupby
from joblib import Parallel, delayed
import multiprocessing

from datetime import datetime

print("begin", datetime.now())

arg_parser = argparse.ArgumentParser(description="Calculate enhancer age.")

arg_parser.add_argument("age_file", help='bed file (enhancers to age) w/ full path')

args = arg_parser.parse_args()

# CONSTANTS
F = args.age_file

agepath = "/".join(F.split("/")[:-1])+"/"
sid = (F.split("/")[-1]).split(".")[0]
print(agepath)

def assemble_breaks(mrca_list, test_id):
    
    ages, arch = np.unique(mrca_list, return_inverse = True) # get the unique ages, architecture order
    arch = [k for k, g in groupby(arch)] # get the architecture w/ itertools groupby

    # create core_remodeling binary
    if len(ages) ==1: 
        core_remodeling = 0
    elif len(ages)>1:
        core_remodeling = 1

    newdf = pd.DataFrame({"test_id": [test_id],
                      "max_seg": [len(arch)], 
                      "core_remodeling": [core_remodeling],
                      "arch":[arch],
                      "max_age": [max(ages).round(3)]})
    return newdf


def format_df(df, agepath, sid):    

    #drop axis
    df = df.drop(["strand","ref","num_species", "patr"], axis = 1)

    # drop syn_blocks > 5 bp in len
    df = df[df.len_syn_overlap > 5]

    # add a test_id
    df["test_id"] = df["chr_enh"] + ":" + df["start_enh"].map(str) + "-" + df["end_enh"].map(str) 

    # figure out difference between enh and syn start
    # positive = syn block inside enh
    # negative = syn block starts before enh - needs to be trimmed!
    df["startdif"] = df.start_enh - df.start_syn

    # trim the start positions (when syn start is before enh start)
    df.loc[df.startdif>0, "start_syn"] = df.start_enh.copy() 

    # figure out difference between enh and syn start
    # positive = syn inside enh
    # negative = syn ends after enh - needs to be trimmed!
    df["enddif"] = df.end_enh - df.end_syn

    # trim the end positions (when syn end is after enh end)
    df.loc[df.enddif<0, "end_syn"] = df.end_enh.copy() 

    collection_dict = {} # collect unique enhancer architectures, max ages, core_remodeling

    # iterate through test_ids
    for test_id in df.test_id.unique():

        if test_id not in collection_dict.keys():

            mrca_list = df.loc[df.test_id==test_id, "mrca"].to_list() # get a list of ages.
            newdf = assemble_breaks(mrca_list, test_id) # assemble the breaks!!!
            collection_dict[test_id] = newdf
    
    archdf = pd.concat(collection_dict.values())
    originaldf = df[["chr_enh", "start_enh", "end_enh", "id", "test_id"]].drop_duplicates()
    results = pd.merge(originaldf, archdf, how = "left")
    chrnum =results.chr_enh.unique().item()
    print('results! %s' % chrnum)
    results.to_csv("%s%s_%s_parallel.txt" % (agepath, chrnum, sid), sep ='\t', header = True, index =False)
    return results

### SPLIT BY CHROMOSOME, TEST PARALLEL ###

chr_dict = {}
df = pd.read_csv(F, sep ='\t', header = None)

df.columns = ["chr_enh", "start_enh", "end_enh", "id",  
                        "chr_syn", "start_syn","end_syn","strand",
                        "ref","num_species", "len_syn","mrca",
                        "patr","len_syn_overlap"]
                        
for chrnum in df.chr_enh.unique():
    chrdf = df.loc[df.chr_enh ==chrnum]
    chr_dict[chrnum] = chrdf
print(chr_dict.keys())

### RUN PARALLEL ###

num_cores = multiprocessing.cpu_count()
num_cores = len(chr_dict.values())

# print start time
start = datetime.now()
print("start parallel", start)

# run parallel jobs
results = Parallel(n_jobs=num_cores, verbose=100, prefer="threads")(delayed(format_df)(i, agepath, sid) for i in chr_dict.values())

# how long did the run take?
end = datetime.now()
print("end parallel", datetime.now())
print("total time", end - start)

catf = "%sshuf-trimmed-310_%s_parallel_breaks.bed" %(agepath, (sid.split("-")[1]).split(".")[0])
cmd = "cat %schr*_%s_parallel.txt > %s" % (agepath, sid, catf)
os.system(cmd)

rm_cmd = "rm %schr*_%s_parallel.txt" % (agepath, sid) # remove all the parallel files
os.system(rm_cmd)

rmF_cmd = "rm %s" % (F) # remove the input file
#os.system(rmF_cmd)


