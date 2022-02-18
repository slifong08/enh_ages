#%%

# 20200928
# sarahfong

# generate 2 .bed files
# first .bed file - all human developmental enhancers, both replicates,
# second .bed file - make a consensus file,replicates from same developmental time
# must show peak overlap > 0 bp.

import os
import pandas as pd
import subprocess

path = "/dors/capra_lab/data/enhancers/reilly15/"
fs = glob.glob("%sGSM*_Mm_*_H3K27ac*.bed" % (path))
fs

#%% collection dictionaries


fdict = {}
Consfs = {}


#%% functions


def get_consensus(devtime, dataset_id, path):
    id_num = int(dataset_id.split("GSM")[1])

    id_num_b = id_num+1
    consensus_f = "%sconsensus/consensus_Mm_H3K27ac_%s.bed" %(path, devtime)

    if devtime == "Hu_8-5": # fancy formatting for this 8.5pcw
        devtime = "8_5pcw"

    a = "%s%s_Mm_%s_H3K27ac_rep1_regions.bed" % (path, dataset_id, devtime)


    b = "%sGSM%s_Mm_%s_H3K27ac_rep2_regions.bed" % (path,id_num_b, devtime)

    cmd = "bedtools intersect -a %s -b %s -wao > %s" % (a, b, consensus_f)
    print(cmd)
    subprocess.call(cmd, shell = True)

    return consensus_f

def format_df(f):

    # get information about the file
    sid = (f.split("/")[-1]).split(".")[0]#"GSM1554668_Mm_8_5pcw_H3K4me2_rep1_regions.bed"

    dataset_id = (sid.split("_")[0])

    rep = sid.split("_")[4]
    devtime = (sid.split("_")[2])

    print( "GSM", dataset_id, devtime, rep)
    # make a dataframe
    df = pd.read_csv(f, sep = '\t', header = None)

    df.columns = ["chr", "start", "end"]
    df["sid"] = sid
    df["devtime"] = devtime
    df["rep"] = rep
    df["enh_id"] = df.chr + ":" + df.start.map(str)+ "-" + df.end.map(str)
    if "rep1" in sid:
        Consf = get_consensus(devtime, dataset_id, path) # run bedtool intersection to get consensus peaks.
    else:
        Consf = ""
    return sid, df, Consf


def format_Consdf(Consf):

    sid = (Consf.split("/")[-1]).split(".")[0]

    devtime = sid.split("_")[-1]

    df = pd.read_csv(Consf, sep = '\t')
    df.columns = ["chr_r1", "start_r1", "end_r1",\
    "chr_r2", "start_r2", "end_r2", "overlap_len"]
    df["rep1_id"] = df.chr_r1 + ":" + df.start_r1.map(str)+ "-" + df.end_r1.map(str)
    df["rep2_id"] = df.chr_r2 + ":" + df.start_r2.map(str)+ "-" + df.end_r2.map(str)
    df["devtime"] = devtime
    df["consensus"] = 0
    df.loc[df.overlap_len>0, "consensus"] = 1


    return df


#%% # run the functions for all the replicates!


for f in fs:
    sid, df, Consf = format_df(f)
    fdict[sid] = df
    print(Consf)
    if Consf != "":
        Consdf = format_Consdf(Consf)
        Consfs[sid] = Consdf
#%%
fdict.keys()


#%%# concatenate all the replicate dataframes together


df = pd.concat(fdict.values())
df.head()


#%%# concatenate all the consensus dataframes together


consdf = pd.concat(Consfs.values())


consdf.head()



#%% mark all the consensus enhancers for replicate 1


df1 = pd.merge(df.loc[df.rep =="rep1"], consdf[["rep1_id", "devtime","consensus"]],
how = "left", left_on = "enh_id", right_on = "rep1_id")

df1.columns = ['chr', 'start','end', 'sid', 'dev_time', 'rep',
 'enh_id', 'rep1_id', 'devtime', 'consensus']

df1.head()


#%% mark all the consensus enhancers for replicate 2


df2 = pd.merge(df[df.rep =="rep2"], consdf[["rep2_id", "devtime", "consensus"]],
how = "left", left_on = "enh_id", right_on = "rep2_id")

# rename the columns
df2.columns = ['chr', 'start','end', 'sid', 'dev_time', 'rep',
'enh_id',  'rep_id2', "devtime", 'consensus']


df2.head()


#%% merge consensus columsn from replicates


new_consdf = pd.concat([df1 , df2], sort=True)
new_consdf = new_consdf.drop(["rep1_id", "rep_id2"], axis = 1)# drop these column

# rearrange the columns
new_consdf = new_consdf[['chr','start', 'end', 'enh_id', 'consensus', 'dev_time','rep', 'sid']]


#%% save the consensus file


outf = "%sconsensus/consensus_Mm_H3K27ac_allpcw.bed" % path

consensus = new_consdf.loc[new_consdf.consensus == 1] # consensus peaks only

# drop replicate coordinates.
consensus = consensus.drop(["rep", "sid"], axis = 1).drop_duplicates()

# write the file
consensus.to_csv(outf, sep = '\t', header = False, index = False)


#%% write the file with all replicates, consensus or not.


outf = "%sall_Mm_H3K27ac.bed" % path
df.to_csv(outf, sep = '\t', header = False, index = False)
#%%
outf
