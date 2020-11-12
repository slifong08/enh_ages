import os
import sys, traceback
import argparse
from datetime import datetime
import numpy as np
import pandas as pd
from functools import partial
from multiprocessing import Pool
from scipy import stats

# TO RUN

# python script input_file sample_id

# python /dors/capra_lab/fongsl/enh_age/bin/age_enhancers.py UBERON_0002372_tonsil_expressed_enhancers.bed UBERON0002372

#%%
"""
###
#   arguments
###
arg_parser = argparse.ArgumentParser(description="Calculate enhancer age.")

arg_parser.add_argument("age_break", help='break file 1 (enhancers to age) w/ full path')

arg_parser.add_argument("-c", "--cutoff", type=int, default=10000,
                        help='threshold enhancer lengths; default=10000')


arg_parser.add_argument("-d", "--dataset", type=str, default="enhancer",
                        help='dataset label')
arg_parser.add_argument("-hd", "--header", type=bool, default="False",
                        help='does file have header?')

args = arg_parser.parse_args()

###
# CONSTANTS
###
#%%

F = args.age_break
SOURCE_F = args.age_break
LEN_CUTOFF = args.cutoff
DATASET = args.dataset
HEADER = args.header
"""
#%%

F = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/non-genic/breaks/no-exon_E038_parallel_breaks_enh_age_arch_summary_matrix.bed"
SOURCE_F = F
LEN_CUTOFF = 10000
SAMPLE_TEN_PERCENT = 0
DATASET = "E038"
HEADER = False

PATH = "/".join(F.split("/")[:-1])
SAMPLE_ID = (F.split("/")[-1]).split(".")[0]

#F = "%s/breaks/%s_age_breaks.bed" % (PATH, SAMPLE_ID)
OUTPATH = "%s/stats/" % PATH

#%%
###
# functions
###

def mkdir(path):
    if os.path.isdir(path) == False:
        cmd = "mkdir %s" % path
        os.system(cmd)
mkdir(OUTPATH)

def get_relative_simple(sid, pdf):
    simple_break_percentile = 0.5
    if sid in pdf.sid2.unique():
        test = pdf.loc[(pdf.sid2 == sid)]

        if simple_break_percentile in test.pct.unique():
            percentile_used = 0.5
            relative_simple = pdf.loc[(pdf.sid2 == sid) &\
             (pdf.pct == simple_break_percentile), "seg_index"].iloc[0].astype(int)

        else: # idk why some datasets do not have a 0.5 percentile
            relative_simple = pdf.loc[(pdf.sid2 == sid) &\
             (pdf.pct == 0.6), "seg_index"].iloc[0].astype(int)
            percentile_used = 0.6
        return relative_simple, percentile_used
    else:
        return "no_relative_simple", "no_percentile_used"

def format_df(df, sample_id, len_cutoff, syn_gen_bkgd, outpath, header, pdf):

    df.columns = ['chr_enh','start_enh',
                      'end_enh', "enh_id", "peak_reads",  "id", 'core_remodeling',"arch",
                      'seg_index','mrca', 'enh_len', "taxon", "mrca_2", "taxons2",
                      "mya", "mya2", "ratio", "dataset"]

    df["sid"]= sample_id # label sample id
    if "shuf" in sample_id:
        sid = (sample_id.split("_")[1]).split("-")[0]
    else:
        sid = sample_id.split("_")[1]
    relative_simple, percentile_used = get_relative_simple(sid, pdf)
    print(sid, relative_simple)

    df = df.loc[df.enh_len <10000] # autosomes only
    df = df.loc[df.chr_enh !="chrX"] # autosomes only
    if relative_simple != "no_relative_simple":
        # code relative simple
        df["relative_arch"] = "rel_simple"
        df.loc[df.seg_index.astype(float) >= float(relative_simple), "relative_arch"] = "rel_complex"
        df.loc[df.relative_arch == "rel_simple", "core_remodeling"] = 0


    # remove complex enhancers from 0.00 mrca. This is no bueno.
    s = df[["enh_id", "core_remodeling","mrca"]].drop_duplicates()

    s = s.enh_id.loc[(s.mrca==0.000)&(s.core_remodeling==1)].tolist()

    df = df[~df.enh_id.isin(s)] # exlude all the data that is mislabeled.
    df.loc[df.dataset.str.contains("age_breaks"), "dataset"] = "roadmap"
    df.loc[df.dataset.str.contains("shuf-"), "dataset"] = "shuf"
    df.loc[df.dataset.str.contains("ucsc"), "dataset"] = "roadmap"

    return df

def get_data_stats(df, sid, dataset, outpath, syn_gen_bkgd):
    sid = sid + "-relative_simple"
    write_dict = {}


    df["mrca"] = df["mrca"].round(3)
    df["mrca_2"] = df["mrca_2"].round(3)

    ## features related to simple and complex
    df_all = df.groupby(["core_remodeling","dataset"]).describe().reset_index() # groupby simple and complex and get basic stats


    all_dict = {}
    all_val = 0
    for cr, new_df in df_all.groupby(level=0): # break apart the levels of the groupby and collect stats

        a = pd.melt(new_df)
        a.columns = ["var", "varStat", "varVal"]

        a["core_remodeling"] = cr
        a["dataset"] = new_df.dataset.iloc[0]

        a["sid"] = sid

        key = str(cr)+ "-"+ str(all_val)
        all_val +=1
        all_dict[key] = a

    df_all_ = pd.concat(all_dict.values())

    write_dict["enh_len_stats"] = df_all_

    df_all_.to_csv("%s%s-df_basic_arch_stats.tsv" % (outpath, sid) , sep = '\t', header= True, index = False)

    ### LENGTHS AND AGES ###

    new_lens = df[["enh_id", "core_remodeling", "enh_len", "mrca_2", "dataset"]].drop_duplicates()
    new_lens.columns = ["id", "core_remodeling", "len", "mrca_2", "dataset"]

    lens = pd.merge(new_lens, syn_gen_bkgd)

    ### ENHANCERS LEN DESC ALL

    desc_all = lens.groupby(["mrca_2", "dataset"])["len"].describe().reset_index()
    desc_all = pd.melt(desc_all, id_vars=['mrca_2', "dataset"], value_vars=['count', 'mean', 'std', 'min', '25%', '50%', '75%', 'max'])

    desc_all.columns = ["mrca_2", "dataset", "len_cat", "len_stat"]

    desc_all["sid"] = sid
    write_dict["enh_len_stats"] = desc_all

    desc_all.to_csv("%s%s-enh_len_stats.tsv" % (outpath, sid) , sep = '\t', header= True, index = False)

    ### ENHANCERS LEN DESC W/ ENHANCER ARCH, AGE
    desc_arch = lens.groupby(["mrca_2", "core_remodeling", "dataset"])["len"].describe().reset_index()
    desc_arch = pd.melt(desc_arch, id_vars=["mrca_2", "core_remodeling", "dataset"], value_vars=['count', 'mean', 'std', 'min', '25%', '50%', '75%', 'max'])
    desc_arch.columns = ["mrca_2","core_remodeling", "dataset", "len_cat", "len_stat", ]

    desc_arch["sid"] = sid
    write_dict["enh_arch_len_stats"] = desc_arch

    desc_arch.to_csv("%s%s-enh_arch_len_stats.tsv" % (outpath, sid) , sep = '\t', header= True, index = False)


    ### REGRESSION MRCA x LEN - calculate slope of the line


    simple_enh = lens.loc[(lens.core_remodeling ==0) ]
    xs = simple_enh.mrca_2
    ys = simple_enh.len

    slope_ss, intercept_ss, r_value_ss, p_value_ss, std_err_ss = stats.linregress(xs,ys)

    complex_enh = lens.loc[(lens.core_remodeling !=0)]
    xc = complex_enh.mrca_2
    yc = complex_enh.len

    slope_cs, intercept_cs, r_value_cs, p_value_cs, std_err_cs = stats.linregress(xc,yc)



    lin_regression = pd.DataFrame({"m": [slope_ss, slope_cs, ],
                                       "b": [intercept_ss, intercept_cs, ],
                                       "r": [r_value_ss,r_value_cs,],
                                       "p": [p_value_ss, p_value_cs,],
                                      "arch":["simple", "complexenh",],
                                      "sid": [sid, sid,],
                                      "dataset" : [sid, sid]})


    write_dict["linear_regression"] = lin_regression
    lin_regression.to_csv("%s%s-linear_regression.tsv" % (outpath, sid) , sep = '\t', header= True, index = False)

    ###  break frequency ###
    breaks = df[["enh_id", "seg_index", "dataset"]].drop_duplicates()
    breaks_freq = breaks.groupby(["dataset","seg_index"])["enh_id"].count().reset_index()
    totals =  breaks.groupby("dataset")["enh_id"].count().reset_index()
    totals.columns = ["dataset", 'totals']
    breaks_freq  = pd.merge(breaks_freq, totals, how = "left")
    breaks_freq["freq"] = breaks_freq["enh_id"].divide(breaks_freq["totals"])
    breaks_freq["sid"] = sid

    write_dict["break_freq"] = breaks_freq

    breaks_freq.to_csv("%s%s-break_freq.tsv" % (outpath, sid) , sep = '\t', header= True, index = False)

    ###  break age ###
    breaks_mrca = df[["enh_id", "seg_index", "mrca_2", "dataset"]].drop_duplicates()
    breaks_mrca = breaks_mrca.groupby(["mrca_2", "dataset"])["seg_index"].describe().reset_index()
    breaks_mrca = pd.melt(breaks_mrca, id_vars=["mrca_2", "dataset"], value_vars=['count', 'mean', 'std', 'min', '25%', '50%', '75%', 'max'])
    breaks_mrca.columns = ["mrca_2", "dataset", "break_cat", "break_stat", ]
    breaks_mrca["sid"] = sid

    write_dict["breaks_mrca"] = breaks_mrca

    breaks_mrca.to_csv("%s%s-breaks_mrca.tsv" % (outpath, sid) , sep = '\t', header= True, index = False)

    ### enh age, arch dataset frequencies ###
    enh_age_freq = lens.groupby(["mrca_2", "core_remodeling", "dataset"])["id"].count().reset_index()

    enh_age_freq.columns = ["mrca_2", "core_remodeling", "dataset", "count_ids"]

    enh_totals = enh_age_freq.groupby("dataset")["count_ids"].sum().reset_index()
    enh_totals.columns= ["dataset", "totals"]
    enh_age_freq = pd.merge(enh_age_freq, enh_totals, how = "left")

    enh_age_freq["frac_of_total"] = enh_age_freq["count_ids"].divide(enh_age_freq.totals).round(3) # get fraction of total

    totals = enh_age_freq.groupby(["core_remodeling","dataset"])["count_ids"].sum().reset_index()
    totals.columns= ["core_remodeling","dataset", "arch_totals"]


    enh_age_freq = pd.merge(enh_age_freq, totals, how = "left")

    enh_age_freq["frac_of_arch"] = enh_age_freq.count_ids.divide(enh_age_freq["arch_totals"]).round(3) # get fraction of simple

    enh_age_freq["total_freq"]= 1
    enh_age_freq = pd.merge(enh_age_freq, syn_gen_bkgd[["mrca_2", "taxon2"]], how = "left", on = "mrca_2").drop_duplicates()

    enh_age_freq["sid"] = sid

    write_dict["enh_age_freq"] = enh_age_freq
    enh_age_freq.to_csv("%s%s-enh_age_freq.tsv" % (outpath, sid), sep = '\t', header= True, index = False)


    return write_dict


#    return desc_all, desc_arch, syn_lens, lin_regression, breaks_freq, breaks_mrca,  enh_age_freq ,syn_arch


def main(argv):
    syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
    syn_gen_bkgd = pd.read_csv(syn_gen_bkgd_file, sep = '\t')
    syn_gen_bkgd["mrca"]=syn_gen_bkgd["mrca"].round(3)
    syn_gen_bkgd["mrca_2"]=syn_gen_bkgd["mrca_2"].round(3)
    syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2"]]

    desc_file = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/roadmap_hg19_sample_id_desc.csv"
    desc_df = pd.read_csv(desc_file, header = None)

    pdf = pd.read_csv("%s/all_noexon_ROADMAP_breaks_percentiles.bed"% PATH, sep = '\t')


    if HEADER == True:
        df = pd.read_csv(F, sep='\t', low_memory=False)
    else:
        df = pd.read_csv(F, sep='\t', header = None, low_memory=False)


    formatted_df = format_df(df, SAMPLE_ID, LEN_CUTOFF, syn_gen_bkgd, OUTPATH, HEADER, pdf)

    #splitup dataframes into enhancer and shuffle
    shuf = formatted_df.loc[formatted_df.dataset =="shuf"]
    df = formatted_df.loc[formatted_df.dataset !="shuf"]

    statsdf = get_data_stats(df, SAMPLE_ID, DATASET, OUTPATH, syn_gen_bkgd)

    shuf_sid = "shuf-"+ SAMPLE_ID
    statsdf = get_data_stats(shuf, shuf_sid, DATASET, OUTPATH, syn_gen_bkgd)


if __name__ == "__main__":
    main(sys.argv[1:])
#%%
syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd = pd.read_csv(syn_gen_bkgd_file, sep = '\t')
syn_gen_bkgd["mrca"]=syn_gen_bkgd["mrca"].round(3)
syn_gen_bkgd["mrca_2"]=syn_gen_bkgd["mrca_2"].round(3)
syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2"]]

desc_file = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/roadmap_hg19_sample_id_desc.csv"
desc_df = pd.read_csv(desc_file, header = None)

pdf = pd.read_csv("%s/all_noexon_ROADMAP_breaks_percentiles.bed"% PATH, sep = '\t')


if HEADER == True:
    df = pd.read_csv(F, sep='\t', low_memory=False)
else:
    df = pd.read_csv(F, sep='\t', header = None, low_memory=False)


formatted_df = format_df(df, SAMPLE_ID, LEN_CUTOFF, syn_gen_bkgd, OUTPATH, HEADER, pdf)
#%%
df = formatted_df.copy()
sid = sid + "-relative_simple"
write_dict = {}


df["mrca"] = df["mrca"].round(3)
df["mrca_2"] = df["mrca_2"].round(3)

## features related to simple and complex
df_all = df.groupby(["core_remodeling","dataset"]).describe().reset_index() # groupby simple and complex and get basic stats
df_all.head()
#%%
all_dict = {}
all_val = 0
for cr, new_df in df_all.groupby(level=0): # break apart the levels of the groupby and collect stats

    a = pd.melt(new_df)
    a.columns = ["var", "varStat", "varVal"]

    a["core_remodeling"] = cr
    a["dataset"] = new_df.dataset.iloc[0]

    a["sid"] = sid

    key = str(cr)+ "-"+ str(all_val)
    all_val +=1
    all_dict[key] = a

df_all_ = pd.concat(all_dict.values())

write_dict["enh_len_stats"] = df_all_
