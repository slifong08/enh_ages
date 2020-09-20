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


###
#   arguments
###
arg_parser = argparse.ArgumentParser(description="Calculate enhancer age.")

arg_parser.add_argument("age_break", help='break file 1 (enhancers to age) w/ full path')

arg_parser.add_argument("-c", "--cutoff", type=int, default=10000,
                        help='threshold enhancer lengths; default=10000')

arg_parser.add_argument("-s", "--sample01", type=int, default=0,
                        help='sample 10% of dataset')

arg_parser.add_argument("-d", "--dataset", type=str, default="enhancer",
                        help='dataset label')
arg_parser.add_argument("-hd", "--header", type=str, default="yes",
                        help='does file have header?')

args = arg_parser.parse_args()

###
# CONSTANTS
###

F = args.age_break
SOURCE_F = args.age_break
LEN_CUTOFF = args.cutoff
SAMPLE_TEN_PERCENT = args.sample01
DATASET = args.dataset
HEADER = args.header

PATH = "/".join(F.split("/")[:-1])
SAMPLE_ID = (F.split("/")[-1]).split(".")[0]
#F = "%s/breaks/%s_age_breaks.bed" % (PATH, SAMPLE_ID) 
OUTPATH = "%s/stats/" % PATH

### 
# functions
###

def mkdir(path):
    if os.path.isdir(path) == False:
        cmd = "mkdir %s" % path
        os.system(cmd)
mkdir(OUTPATH)
        

def format_df(df, sample_id, sample_ten_percent, len_cutoff, syn_gen_bkgd, outpath, header):
    if header == "yes":
        df = df[['chr_enh','start_enh',
                      'end_enh', "test_id",'mrca', 'seg_index',"enh_len"]]
   
    df["enh_id"] = df.test_id
    df["core_remodeling"] = 0 # make a core_remodeling column
    df.loc[df.seg_index.astype(int) >1, "core_remodeling"] = 1
    
    df = df.loc[df.enh_len <10000] # autosomes only
    df = df.loc[df.chr_enh !="chrX"] # autosomes only
    
    # remove complex enhancers from 0.00 mrca. This is no bueno. 
    s = df.groupby(["enh_id", "core_remodeling"])["mrca"].max().reset_index()
    s = s.enh_id.loc[(s.mrca==0.000)&(s.core_remodeling==1)].tolist()

    df = df[~df.enh_id.isin(s)] # exlude all the data that is mislabeled. 
    if sample_ten_percent == 1:
            
        if len_cutoff > 0:
            enh_ids = df.loc[df["enh_len"]<= len_cutoff, "enh_id"].unique() # get a unique list of enhancer ids
        else:
            enh_ids = df["enh_id"].unique() # get a unique list of enhancer ids

        tenpercent = int(round(len(enh_ids)/10, 0)) # int for 10% of ids
            
        print("0.1 of %s enhancers" % sample_id, len(enh_ids), "is", tenpercent)

        sampled = np.random.choice(enh_ids,tenpercent) # select 10% of ids at random

        sampled_df = pd.DataFrame({"enh_id": sampled}) # make a new dataframe

        # get syntenic block information for 10% shuffled enhancer ids

        df = pd.merge(sampled_df, df,  how = "left", on = "enh_id") # write over shuffle_df.

    else:
        if len_cutoff > 0:
            df = df.loc[df.enh_len<= len_cutoff]
    
    ###
    # break pseudocount
    # enh len
    # percent of enhancer

    # enh id
    #
    # code
    #
    # filter synetnic blocks <9bp long
    # merge taxon, mrca_2
    ###

    df["enh_id"] = df["chr_enh"]+ ":" + df["start_enh"].map(str) + "-" + df["end_enh"].map(str)


    df["arch"] = ""
    df.loc[(df["core_remodeling"] == 1),"arch"] = "complexenh"
    df.loc[(df["core_remodeling"] == 0), "arch"] = "simple"

    ###FILTER###

    df["mrca"]=df["mrca"].round(3)

    df = pd.merge(df, syn_gen_bkgd, how = "left", on = "mrca")

    before_shape = df.shape

    df = df.drop_duplicates()
    print("before", before_shape, "after dropping duplicates", df.shape) 
   
    return df

def get_data_stats(df, sid, dataset, outpath, syn_gen_bkgd):
    
    write_dict = {}
    
    ### SHUF len per age ###
    df["mrca"] = df["mrca"].round(3)
    df["mrca_2"] = df["mrca_2"].round(3)

    df_all = df.groupby("core_remodeling").describe().reset_index() # groupby simple and complex and get basic stats

    all_dict = {}
    all_val = 0
    for cr, new_df in df_all.groupby(level=0): # break apart the levels of the groupby and collect stats
        a = pd.melt(new_df)
        a.columns = ["var", "varStat", "varVal"]
        a["core_remodeling"] = cr

        a["sid"] = sid
        a["dataset"] = dataset
        key = str(cr)+ "-"+ str(all_val)
        all_val +=1
        all_dict[key] = a
        
    df_all_ = pd.concat(all_dict.values())
        
    write_dict["enh_len_stats"] = df_all_
    df_all_.to_csv("%s%s-df_basic_arch_stats.tsv" % (outpath, sid) , sep = '\t', header= True, index = False)
    
    new_lens = df.groupby(["enh_id", "core_remodeling", "enh_len"])["mrca_2"].max().reset_index()
    new_lens.columns = ["id", "core_remodeling", "len", "mrca_2"]

    lens = pd.merge(new_lens, syn_gen_bkgd)

    ### ENHANCERS LEN DESC ALL 

    desc_all = lens.groupby("mrca_2")["len"].describe().reset_index()
    desc_all = pd.melt(desc_all, id_vars=['mrca_2'], value_vars=['count', 'mean', 'std', 'min', '25%', '50%', '75%', 'max'])

    desc_all.columns = ["mrca_2", "len_cat", "len_stat"]
    desc_all["dataset"] = dataset
    desc_all["sid"] = sid
    write_dict["enh_len_stats"] = desc_all
    desc_all.to_csv("%s%s-enh_len_stats.tsv" % (outpath, sid) , sep = '\t', header= True, index = False)

    ### ENHANCERS LEN DESC W/ ENHANCER ARCH, AGE
    desc_arch = lens.groupby(["mrca_2", "core_remodeling"])["len"].describe().reset_index()
    desc_arch = pd.melt(desc_arch, id_vars=["mrca_2", "core_remodeling"], value_vars=['count', 'mean', 'std', 'min', '25%', '50%', '75%', 'max'])
    desc_arch.columns = ["mrca_2","core_remodeling", "len_cat", "len_stat", ]
    desc_arch["dataset"] = dataset
    desc_arch["sid"] = sid
    write_dict["enh_arch_len_stats"] = desc_arch
    desc_arch.to_csv("%s%s-enh_arch_len_stats.tsv" % (outpath, sid) , sep = '\t', header= True, index = False)
    
    ### REGRESSION MRCA x LEN - calculate slope of the line
    xs = lens.mrca_2.loc[lens.core_remodeling ==0]
    ys = lens.len.loc[lens.core_remodeling ==0]
    print("XS,YS", xs, ys)
    slope_ss, intercept_ss, r_value_ss, p_value_ss, std_err_ss = stats.linregress(xs,ys)

    xc = lens.mrca_2.loc[lens.core_remodeling ==1]
    yc = lens.len.loc[lens.core_remodeling ==1]
    print(lens.head())
    slope_cs, intercept_cs, r_value_cs, p_value_cs, std_err_cs = stats.linregress(xc,yc)
    
    lin_regression = pd.DataFrame({"m": [slope_ss, slope_cs, ],
                                       "b": [intercept_ss, intercept_cs,],
                                       "r": [r_value_ss,r_value_cs,],
                                       "p": [p_value_ss, p_value_cs, ],
                                      "arch":["simple", "complexenh", ],
                                      "sid": [sid, sid,], 
                                      "dataset" : [dataset, dataset,]})
    
    write_dict["linear_regression"] = lin_regression
    lin_regression.to_csv("%s%s-linear_regression.tsv" % (outpath, sid) , sep = '\t', header= True, index = False)

    ###  break frequency ###
    breaks = df.groupby(["enh_id"])["seg_index"].max().reset_index()
    breaks_freq = breaks.groupby("seg_index")["enh_id"].count().reset_index()
    totals = len(breaks)
    breaks_freq["freq"] = breaks_freq["enh_id"].divide(totals)
    breaks_freq["sid"] = sid
    breaks_freq["dataset"] = dataset

    write_dict["break_freq"] = breaks_freq
    breaks_freq.to_csv("%s%s-break_freq.tsv" % (outpath, sid) , sep = '\t', header= True, index = False)
    
    ###  break age ###
    breaks_mrca = df.groupby(["enh_id"])["seg_index", "mrca_2"].max().reset_index()
    breaks_mrca = breaks_mrca.groupby(["mrca_2"])["seg_index"].describe().reset_index()
    breaks_mrca = pd.melt(breaks_mrca, id_vars=["mrca_2",], value_vars=['count', 'mean', 'std', 'min', '25%', '50%', '75%', 'max'])
    breaks_mrca.columns = ["mrca_2", "break_cat", "break_stat"]
    breaks_mrca["sid"] = sid
    breaks_mrca["dataset"] = dataset
    
    write_dict["breaks_mrca"] = breaks_mrca
    breaks_mrca.to_csv("%s%s-breaks_mrca.tsv" % (outpath, sid) , sep = '\t', header= True, index = False)
    
    ### enh age, arch dataset frequencies ###
    enh_age_freq = lens.groupby(["mrca_2", "core_remodeling"])["id"].count().reset_index()

    enh_age_freq.columns = ["mrca_2", "core_remodeling", "count_ids"]
    
    enh_age_freq["totals"] = enh_age_freq["count_ids"].sum()
    enh_age_freq["frac_of_total"] = enh_age_freq["count_ids"].divide(enh_age_freq.totals).round(3) # get fraction of total

    simple_total = enh_age_freq["count_ids"].loc[enh_age_freq.core_remodeling==0].sum()

    enh_age_freq["frac_of_arch"] = enh_age_freq["count_ids"].divide(simple_total).round(3) # get fraction of simple

    complex_total = enh_age_freq["count_ids"].loc[enh_age_freq.core_remodeling==1].sum()
    enh_age_freq["frac_of_arch"].loc[enh_age_freq.core_remodeling == 1] = enh_age_freq["count_ids"].divide(complex_total).round(3) # get fraction of complex enhancers

    enh_age_freq["total_freq"]= 1
    enh_age_freq = pd.merge(enh_age_freq, syn_gen_bkgd[["mrca_2", "taxon2"]], how = "left", on = "mrca_2").drop_duplicates()

    enh_age_freq["dataset"] = dataset
    enh_age_freq["sid"] = sid

    write_dict["enh_age_freq"] = enh_age_freq
    enh_age_freq.to_csv("%s%s-enh_age_freq.tsv" % (outpath, sid), sep = '\t', header= True, index = False)
    
    return write_dict

def main(argv):
    syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
    syn_gen_bkgd = pd.read_csv(syn_gen_bkgd_file, sep = '\t')
    syn_gen_bkgd["mrca"]=syn_gen_bkgd["mrca"].round(3)
    syn_gen_bkgd["mrca_2"]=syn_gen_bkgd["mrca_2"].round(3)
    syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2"]]

    desc_file = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/roadmap_hg19_sample_id_desc.csv"
    desc_df = pd.read_csv(desc_file, header = None)

    if HEADER == "yes":
        df = pd.read_csv(F, sep='\t', low_memory=False)
    else:
        df = pd.read_csv(F, sep='\t', header = None, low_memory=False)

    formatted_df = format_df(df, SAMPLE_ID, SAMPLE_TEN_PERCENT, LEN_CUTOFF, syn_gen_bkgd, OUTPATH, HEADER)

    statsdf = get_data_stats(formatted_df, SAMPLE_ID, DATASET, OUTPATH, syn_gen_bkgd)

if __name__ == "__main__":
    main(sys.argv[1:])