
from itertools import groupby
from joblib import Parallel, delayed
import multiprocessing
import os, sys
import pandas as pd
import subprocess 

##
# input should be the *_ages.bed file
##

f = sys.argv[1] 

### 
# Functions 
###


def makeDf(f):

    # column names
    columns = ["chr_enh", "start_enh", "end_enh", "chr_syn", "start_syn", "end_syn", "mrca"]

    # read dataframe
    df = pd.read_csv(f, sep = '\t', header = None, usecols= [0, 1, 2, 4, 5, 6, 11])

    # name columns
    df.columns = columns
    
    # blocks that cannot be aged reassigned as human age. 
    df.loc[df.mrca == ".", "mrca"] = 0

    # round mrca
    df.mrca = df.mrca.astype(float).round(3)

    # add enhancer id
    df["enh_id"] = df.chr_enh + ":" + df.start_enh.map(str) + "-" + df.end_enh.map(str)
    
    return df


def assemble_break(df, enh_id, sid, path):

    # base dataframe 
    test = df.loc[df["enh_id"]== enh_id]\
        .drop_duplicates()\
        .sort_values(["chr_syn","start_syn","end_syn"]).reset_index()

    # new dataframe for reassembled breaks
    new_data = test[["chr_enh","start_enh", "end_enh", "enh_id", "chr_syn"]]\
    .drop_duplicates()

    
    collect_df = pd.DataFrame() # collect results here. 

    
    # round the age and sum the number of rows per age.
    age_seg = [(age, sum(1 for i in rows)) for age,rows in groupby(test["mrca"])] 

    
    start, end = test["start_enh"].astype(int), test["end_enh"].astype(int) # change start and end

    df_index = 0 # placeholder - for syntenic segment index
    seg_count = 0 # placeholder - for syntenic segment count
    core_age = max(i for i, k in age_seg) # GET OLDEST/MAX AGE


    for tup in age_seg: # unpack the list of tuples
        age, idx = tup

        # SIMPLE
        
        if len(age_seg)== 1: 

            new_data["start_syn"] = start
            new_data["end_syn"] = end
            new_data["mrca_seg"] = age
            new_data["seg_index"] = 0 # index syntenic segments
            new_data["core"] = 1 # annotate the core
            new_data["core_remodeling"] = 0 # binary. 0 = simple, 1 = complex
            collect_df= collect_df.append(new_data)

        # COMPLEX
        
        else: 

            new_data["seg_index"]= seg_count # index syntenic segments,
            new_data["core_remodeling"] = 1 # binary. 0 = simple, 1 = complex
            
            # get oldest core age
            if age == core_age: 
                new_data["core"] = 1

            else:
                new_data["core"] = 0
            
            # trim first syntenic block to start
            if seg_count == 0: 
                new_data["start_syn"] = start
                new_data["end_syn"] = test.loc[idx-1, "end_syn"]

                new_data["mrca_seg"] = round(age, 3)
                collect_df= collect_df.append(new_data)

            # trim last syntenic block
            elif seg_count == len(age_seg)-1: 
                new_data["mrca_seg"] = age
                new_data["start_syn"] = test.loc[df_index, "start_syn"]
                new_data["end_syn"] = end
                collect_df= collect_df.append(new_data)

            # deal with all the blocks in between first and last syntenic block
            else: 
                new_data["mrca_seg"] = age
                new_data["start_syn"]= test.loc[df_index, "start_syn"]
                new_data["end_syn"]= test.loc[df_index + idx -1, "end_syn"]
                collect_df= collect_df.append(new_data)

        df_index +=idx # go to next index
        seg_count +=1 # count age segments

        
    # re-arrange dataframe 
    collect_df = collect_df[["chr_syn","start_syn", "end_syn", "enh_id",
                             "chr_enh","start_enh", "end_enh",
                             "seg_index", "core_remodeling", "core",
                             "mrca_seg"]]
        
    # write out file 
    outfile = "%s/%s_%s_syn_breaks.bed" % (path, enh_id, sid)
    collect_df.to_csv(outfile, sep ='\t', header = False, index =False)
    
    return collect_df

### Main ###

# get the dataframe
df = makeDf(f)

path = "/".join(f.split("/")[:-1]) + "/"
sid = (f.split("/")[-1]).split(".")[0]
                

### RUN PARALLEL BREAK ASSEMBLY ###


num_cores = multiprocessing.cpu_count()
print("number of cores", num_cores)

# run parallel jobs

Parallel(n_jobs=num_cores, verbose=100, prefer="threads")(delayed(assemble_break)(df, i, sid, path) for i in df.enh_id.unique())


# concatenate all the files

out_cat = "%ssyn_breaks_%s.bed" % (path, sid)
cat_cmd = "%schr*_%s_syn_breaks.bed > %s" % (path, sid, out_cat)

subprocess.call(cat_cmd, shell = True) 
